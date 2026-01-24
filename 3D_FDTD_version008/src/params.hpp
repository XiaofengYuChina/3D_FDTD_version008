// params.hpp — Simulation configuration and parameter file
//
// This file contains the core simulation parameter structures and context
// management. User-configurable values are read from user_config.hpp.

#pragma once

#include <string>
#include <numbers>
#include <functional>
#include <vector>
#include <memory>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <limits>

#include "global_function.hpp"
#include "user_config.hpp"
#include "boundary.hpp"
#include "structure_material.hpp"
#include "detectors.hpp"
#include "auto_mesh_generator.hpp"

// Modular sources
#include "sources/isource.hpp"
#include "sources/dipole_source.hpp"
#include "sources/plane_wave_source.hpp"

// Modular detectors
#include "detectors/idetector.hpp"
#include "detectors/mesh_detector.hpp"
#include "detectors/field_movie_2d.hpp"
#include "detectors/point_field_detector.hpp"

namespace Config {

    // -------------------------- Basic Parameters (from PhysConst) --------------------------

    inline constexpr real PI   = PhysConst::PI;
    inline constexpr real c0   = PhysConst::C0;     // Speed of light (m/s)
    inline constexpr real eps0 = PhysConst::EPS0;   // Vacuum permittivity (F/m)
    inline constexpr real mu0  = PhysConst::MU0;    // Vacuum permeability (H/m)

    // -------------------------- Mesh mode enumeration --------------------------
    enum class MeshMode {
        UNIFORM,            // Uniform mesh (constant dx, dy, dz)
        AUTO_NONUNIFORM     // Automatic non-uniform mesh (ppw-based, interface snapping)
    };

    // -------------------------- Unified simulation parameters --------------------------
    struct SimulationParams {

        // ===== Mesh mode selection =====
        MeshMode mesh_mode = MeshMode::AUTO_NONUNIFORM;

        // Uniform mesh defaults (for UNIFORM mode)
        real dx_uniform = 10e-9;
        real dy_uniform = 10e-9;
        real dz_uniform = 10e-9;

        // Auto mesh configuration (for AUTO_NONUNIFORM mode)
        AutoMeshConfig auto_mesh_config;

        // Grid spacing storage
        GridSpacing grid_spacing;

        // Physical domain boundaries (absolute coordinates)
        real domain_x_min = 0;
        real domain_x_max = 0;
        real domain_y_min = 0;
        real domain_y_max = 0;
        real domain_z_min = 0;
        real domain_z_max = 0;

        // Derived domain sizes (computed from boundaries)
        real domain_size_x() const { return domain_x_max - domain_x_min; }
        real domain_size_y() const { return domain_y_max - domain_y_min; }
        real domain_size_z() const { return domain_z_max - domain_z_min; }

        // ===== Grid and time step (CFL) =====
        // Core dimensions (excluding PML)
        size_t Nx = 100, Ny = 100, Nz = 100;

        // CFL safety factor S (<1), dt = S / c0 / sqrt(1/dx^2 + 1/dy^2 + 1/dz^2)
        real S = 0.9;

        // ===== Time stepping =====
        size_t nSteps = 1000;   // Number of time steps
        size_t saveEvery = 1;   // Save at least one frame/one point every N steps

        // Derived quantities: time step and coefficients used in update formulas
        real dt = 0.0;      // s
        real Ce = 0.0;      // dt/eps0
        real Ch = 0.0;      // dt/mu0

        // ===== Boundary conditions (PEC and CPML) =====
        BoundaryParams bp = [] {
            BoundaryParams p;
            p.type = BcType::CPML_RC;
            p.cpml.npml = 16;
            p.cpml.m = 3.0;
            p.cpml.Rerr = 1e-8;
            p.cpml.alpha0 = 0.2;
            p.cpml.kappa_max = 10.0;
            p.cpml.alpha_linear = true;
            return p;
        }();

        // Output file path: frames/<run_tag>/...        
        std::string run_tag = "3D_FDTD_v008_output";

        // ===== Source type: Ricker/Gaussian pulse =====
        // Peak Current (A)
        real Iz_peak = 1e-5;

        // Frequency setting: f0 (priority) or lambda0
        real f0 = 0.0;                  // Hz
        real lambda0 = 500e-9;          // m (default 500nm)

        // Envelope width parameter (tau or df_FWHM)
        real tau_src = 5e-15;           // s
        real df_fwhm = 2.4e13;          // Hz

        // Time delay t0 = t0_factor * tau_eff, if t0_factor<=0 then 3/f0        
        real t0_factor = 0.4;

        // Source pulse parameters
        real tau_eff = 0.0;     // s   effective pulse width      
        real omega0 = 0.0;      // rad/s  central angular frequency   
        real t0 = 0.0;          // s   time delay

        // 2D slice detector: z-slice position in PHYSICAL coordinates (meters)
        // Must be within [domain_z_min, domain_z_max]
        real z_slice_z = 0;

        // 2D frame file pattern: whether to write float64
        std::string framePattern = "ez_%04d.raw";
        bool writeFloat64 = true;

        // 1D point probe: position in PHYSICAL coordinates (meters)
        // Must be within domain boundaries
        real probe_x = 0;
        real probe_y = 0;
        real probe_z = 0;

        // Source position in PHYSICAL coordinates (meters)
        // Must be within domain boundaries
        real source_x = 0;
        real source_y = 0;
        real source_z = 0;

        // Function to build boundaries and objects in the simulation scene
        std::function<void(StructureScene&)> build_scene;

        // ====== Initialization/derived quantity calculation ======
        void finalize() {
            // Simple initial dt estimate (will be recalculated in make_context)
            real dx_est = dx_uniform;
            real dy_est = dy_uniform;
            real dz_est = dz_uniform;
            
            const real sum_inv2 = (1.0/(dx_est*dx_est)) + 
                                 (1.0/(dy_est*dy_est)) + 
                                 (1.0/(dz_est*dz_est));
            dt = S / (c0 * std::sqrt(sum_inv2));
            Ce = dt / eps0;
            Ch = dt / mu0;
            
            // 2) Frequency setting: priority f0 (directly given) > lambda0 (calculate) > default
            if (f0 <= 0.0) {
                if (lambda0 > 0.0) {
                    f0 = c0 / lambda0;
                } else {
                    const real lam_default = 500e-9;
                    f0 = c0 / lam_default;
                }
            }
            omega0 = 2.0 * std::numbers::pi_v<double> * f0;

            // 3) Envelope width: tau_src priority, then df_FWHM, if both are 0 default to 20*dt
            if (tau_src > 0.0) {
                tau_eff = tau_src;
            } else if (df_fwhm > 0.0) {
                tau_eff = tau_from_df_fwhm(df_fwhm);
            } else {
                tau_eff = 20.0 * dt;
            }

            // 4) Source center time
            if (t0_factor > 0.0) {
                t0 = t0_factor * tau_eff;
            } else {
                t0 = 3.0 / f0;
            }
        }
    };

    // -------------------------- Simulation context --------------------------
    struct SimContext {
        // Grid dimensions (x/y/z) including core and PML regions
        size_t NxT{}, NyT{}, NzT{};
        size_t Nx_core{}, Ny_core{}, Nz_core{};
        size_t npml{};

        // Grid spacing and time step
        GridSpacing grid_spacing;
        real dt{};

        // Boundary conditions and system objects        
        std::unique_ptr<IBoundary> bc;
        StructureScene scene;
        MaterialGrids  mats;

        // Output directory for simulation results
        std::filesystem::path out_root;
    };

    // -------------------------- Runtime --------------------------
    struct Runtime {
        std::vector<std::unique_ptr<Sources::ISource>> sources;

        // Unified detector callback functions, called once per loop iteration
        std::function<void(size_t,
            const std::vector<real>&, const std::vector<real>&, const std::vector<real>&,
            const std::vector<real>&, const std::vector<real>&, const std::vector<real>&)>
            after_H =
            [](auto, const auto&, const auto&, const auto&, const auto&, const auto&, const auto&) {};

        std::function<void(size_t,
            const std::vector<real>&, const std::vector<real>&, const std::vector<real>&,
            const std::vector<real>&, const std::vector<real>&, const std::vector<real>&)>
            after_E =
            [](auto, const auto&, const auto&, const auto&, const auto&, const auto&, const auto&) {};
    };

    // Type aliases for structure, source, and detector setup functions
    using StructurePackage = std::function<void(StructureScene& scene)>;
    using SourcePackage = std::function<void(const SimContext& ctx, Runtime& rt)>;
    using DetectorPackage = std::function<void(const SimContext& ctx, Runtime& rt)>;

    // C++17 inline vectors for storing global structure, source, and detector packages
    inline std::vector<StructurePackage> g_structure_pkgs{};
    inline std::vector<SourcePackage>    g_source_pkgs{};
    inline std::vector<DetectorPackage>  g_detector_pkgs{};

    // Structure bounds for auto mesh generation (pre-registered before make_context)
    inline std::vector<StructureBounds> g_structure_bounds{};

    // Helper function to register structure bounds for auto mesh
    inline void register_structure_bounds(real x_min, real x_max,
                                          real y_min, real y_max,
                                          real z_min, real z_max,
                                          real n_material) {
        g_structure_bounds.push_back({x_min, x_max, y_min, y_max, z_min, z_max, n_material});
    }

    // Clear all registered structure bounds
    inline void clear_structure_bounds() {
        g_structure_bounds.clear();
    }

    // --- Main functions to apply structure, source, and detector packages ---
    inline void apply_structure_packages(SimContext& ctx) {
        for (auto& p : g_structure_pkgs) p(ctx.scene);
    }
    inline void apply_source_packages(const SimContext& ctx, Runtime& rt) {
        for (auto& p : g_source_pkgs) p(ctx, rt);
    }
    inline void apply_detector_packages(const SimContext& ctx, Runtime& rt) {
        for (auto& p : g_detector_pkgs) p(ctx, rt);
    }

    // --- Build simulation context -> initialize scene information ---
    inline SimContext make_context(SimulationParams& P) {
        SimContext ctx;
        std::cout << "\n========================================\n";
        std::cout << "Building Simulation Context\n";
        std::cout << "========================================\n";

        const size_t npml = P.bp.cpml.npml;

        // ===== MESH GENERATION BASED ON MODE =====
        switch (P.mesh_mode) {
            case MeshMode::AUTO_NONUNIFORM: {
                std::cout << "[Grid] Mode: AUTO_NONUNIFORM (two-step: physical mesh + PML extension)\n";

                // Configure auto mesh from simulation parameters
                P.auto_mesh_config.lambda_min = P.lambda0;
                if (P.auto_mesh_config.dx_min <= 0) {
                    P.auto_mesh_config.dx_min = P.lambda0 / 100.0;
                }
                if (P.auto_mesh_config.dx_max <= 0) {
                    P.auto_mesh_config.dx_max = P.lambda0 / 4.0;
                }

                // Determine physical domain size from boundaries
                real domain_x = P.domain_size_x();
                real domain_y = P.domain_size_y();
                real domain_z = P.domain_size_z();

                // If domain sizes not specified (boundaries are 0), estimate from Nx/dx_uniform
                if (domain_x <= 0) domain_x = P.Nx * P.dx_uniform;
                if (domain_y <= 0) domain_y = P.Ny * P.dy_uniform;
                if (domain_z <= 0) domain_z = P.Nz * P.dz_uniform;

                std::cout << "[Grid] Physical domain: " << domain_x*1e9 << " x "
                          << domain_y*1e9 << " x " << domain_z*1e9 << " nm\n";
                std::cout << "[Grid] npml: " << npml << ", dx_max: " << P.auto_mesh_config.dx_max*1e9 << " nm\n";
                std::cout << "[Grid] Registered structure bounds: " << g_structure_bounds.size() << "\n";

                // Convert structure bounds to auto mesh format
                auto auto_mesh_structures = convert_to_auto_mesh_bounds(g_structure_bounds);

                // Generate auto mesh using two-step approach:
                // Step 1: Generate physical domain mesh (0 to domain_size)
                // Step 2: Extend with PML on both sides
                AutoMeshGenerator generator(P.auto_mesh_config);
                P.grid_spacing = generator.generate(
                    domain_x, domain_y, domain_z,  // Physical domain sizes (NOT total domain)
                    npml, auto_mesh_structures, 1.0       // npml cells extended on each side
                );

                // Update grid dimensions based on generated mesh
                P.Nx = P.grid_spacing.dx.size() - 2 * npml;
                P.Ny = P.grid_spacing.dy.size() - 2 * npml;
                P.Nz = P.grid_spacing.dz.size() - 2 * npml;
                break;
            }

            case MeshMode::UNIFORM:
            default: {
                std::cout << "[Grid] Mode: UNIFORM\n";
                const size_t NxT = P.Nx + 2 * npml;
                const size_t NyT = P.Ny + 2 * npml;
                const size_t NzT = P.Nz + 2 * npml;

                std::cout << "  dx = " << P.dx_uniform * 1e9 << " nm\n";
                std::cout << "  dy = " << P.dy_uniform * 1e9 << " nm\n";
                std::cout << "  dz = " << P.dz_uniform * 1e9 << " nm\n";

                P.grid_spacing.dx.assign(NxT, P.dx_uniform);
                P.grid_spacing.dy.assign(NyT, P.dy_uniform);
                P.grid_spacing.dz.assign(NzT, P.dz_uniform);

                build_cumulative_bounds(P.grid_spacing.dx, P.grid_spacing.x_bounds);
                build_cumulative_bounds(P.grid_spacing.dy, P.grid_spacing.y_bounds);
                build_cumulative_bounds(P.grid_spacing.dz, P.grid_spacing.z_bounds);
                break;
            }
        }

        // Pre-compute inverse spacing arrays (OPTIMIZATION: eliminate divisions in inner loops)
        P.grid_spacing.compute_inverses();

        // ===== COMPUTE TIME STEP (CFL condition using minimum spacing) =====
        const size_t NxT = P.grid_spacing.dx.size();
        const size_t NyT = P.grid_spacing.dy.size();
        const size_t NzT = P.grid_spacing.dz.size();

        real dx_min = *std::min_element(P.grid_spacing.dx.begin(), P.grid_spacing.dx.end());
        real dy_min = *std::min_element(P.grid_spacing.dy.begin(), P.grid_spacing.dy.end());
        real dz_min = *std::min_element(P.grid_spacing.dz.begin(), P.grid_spacing.dz.end());
        real dx_max = *std::max_element(P.grid_spacing.dx.begin(), P.grid_spacing.dx.end());
        real dy_max = *std::max_element(P.grid_spacing.dy.begin(), P.grid_spacing.dy.end());
        real dz_max = *std::max_element(P.grid_spacing.dz.begin(), P.grid_spacing.dz.end());

        std::cout << "\n[Grid] Final grid dimensions: " << NxT << " x " << NyT << " x " << NzT << "\n";
        std::cout << "[Grid] Core grid: " << (NxT - 2*npml) << " x " << (NyT - 2*npml) << " x " << (NzT - 2*npml) << "\n";
        std::cout << "[Grid] Spacing range:\n";
        std::cout << "  dx: " << dx_min * 1e9 << " - " << dx_max * 1e9 << " nm\n";
        std::cout << "  dy: " << dy_min * 1e9 << " - " << dy_max * 1e9 << " nm\n";
        std::cout << "  dz: " << dz_min * 1e9 << " - " << dz_max * 1e9 << " nm\n";

        // Check grading ratios
        auto check_grading = [](const std::vector<real>& spacing, const char* name) {
            real max_ratio = 1.0;
            for (size_t i = 1; i < spacing.size(); ++i) {
                real ratio = std::max(spacing[i] / spacing[i-1], spacing[i-1] / spacing[i]);
                max_ratio = std::max(max_ratio, ratio);
            }
            std::cout << "  " << name << " max grading ratio: " << max_ratio;
            if (max_ratio > 1.2) {
                std::cout << " (WARNING: >1.2 may cause reflections)";
            }
            std::cout << "\n";
        };
        check_grading(P.grid_spacing.dx, "dx");
        check_grading(P.grid_spacing.dy, "dy");
        check_grading(P.grid_spacing.dz, "dz");

        const real sum_inv2 = (1.0/(dx_min*dx_min)) + (1.0/(dy_min*dy_min)) + (1.0/(dz_min*dz_min));
        P.dt = P.S / (c0 * std::sqrt(sum_inv2));
        P.Ce = P.dt / eps0;
        P.Ch = P.dt / mu0;

        std::cout << "\n[Time] dt = " << P.dt * 1e18 << " as (" << P.dt * 1e15 << " fs)\n";
        std::cout << "[Time] CFL factor S = " << P.S << "\n";

        // ===== CREATE BOUNDARY CONDITIONS =====
        P.Nx = NxT - 2 * npml;
        P.Ny = NyT - 2 * npml;
        P.Nz = NzT - 2 * npml;

        auto bc = make_boundary(P.bp, P.Nx, P.Ny, P.Nz,
                               P.grid_spacing, P.dt, eps0, mu0, c0);

        ctx.bc = std::move(bc);
        ctx.NxT = ctx.bc->NxT();
        ctx.NyT = ctx.bc->NyT();
        ctx.NzT = ctx.bc->NzT();
        ctx.Nx_core = P.Nx;
        ctx.Ny_core = P.Ny;
        ctx.Nz_core = P.Nz;
        ctx.npml = ctx.bc->npml();
        ctx.dt = P.dt;
        ctx.grid_spacing = P.grid_spacing;
        ctx.grid_spacing.npml = ctx.npml;  // Set PML thickness for physical coordinate conversion

        // ===== INITIALIZE SCENE =====
        ctx.scene.NxT = ctx.NxT;
        ctx.scene.NyT = ctx.NyT;
        ctx.scene.NzT = ctx.NzT;
        ctx.scene.Nx_core = ctx.Nx_core;
        ctx.scene.Ny_core = ctx.Ny_core;
        ctx.scene.Nz_core = ctx.Nz_core;
        ctx.scene.npml_cells = ctx.npml;
        ctx.scene.grid_spacing = ctx.grid_spacing;

        // ===== SAVE GRID SPACING FOR VISUALIZATION =====
        {
            namespace fs = std::filesystem;
            fs::create_directories("frames/" + P.run_tag);

            std::ofstream f("frames/" + P.run_tag + "/grid_spacing.txt");
            f << "# Auto-generated grid spacing file\n";
            f << "# Mode: " << (P.mesh_mode == MeshMode::AUTO_NONUNIFORM ? "AUTO_NONUNIFORM" : "UNIFORM") << "\n";
            f << "# dx (one per line, in meters)\n";
            for (auto dx : ctx.grid_spacing.dx) f << dx << "\n";
            f << "# dy\n";
            for (auto dy : ctx.grid_spacing.dy) f << dy << "\n";
            f << "# dz\n";
            for (auto dz : ctx.grid_spacing.dz) f << dz << "\n";
            f << "# x_bounds\n";
            for (auto x : ctx.grid_spacing.x_bounds) f << x << "\n";
            f << "# y_bounds\n";
            for (auto y : ctx.grid_spacing.y_bounds) f << y << "\n";
            f << "# z_bounds\n";
            for (auto z : ctx.grid_spacing.z_bounds) f << z << "\n";
        }

        // ===== COORDINATE SYSTEM DIAGNOSTIC PRINTS (Step A) =====
        {
            const real pml_offset_x = ctx.grid_spacing.pml_offset_x();
            const real pml_offset_y = ctx.grid_spacing.pml_offset_y();
            const real pml_offset_z = ctx.grid_spacing.pml_offset_z();

            const real total_x_min = ctx.grid_spacing.x_bounds.front();
            const real total_x_max = ctx.grid_spacing.x_bounds.back();
            const real total_y_min = ctx.grid_spacing.y_bounds.front();
            const real total_y_max = ctx.grid_spacing.y_bounds.back();
            const real total_z_min = ctx.grid_spacing.z_bounds.front();
            const real total_z_max = ctx.grid_spacing.z_bounds.back();

            const real core_x_min = ctx.grid_spacing.x_bounds[ctx.npml];
            const real core_x_max = ctx.grid_spacing.x_bounds[ctx.NxT - ctx.npml];
            const real core_y_min = ctx.grid_spacing.y_bounds[ctx.npml];
            const real core_y_max = ctx.grid_spacing.y_bounds[ctx.NyT - ctx.npml];
            const real core_z_min = ctx.grid_spacing.z_bounds[ctx.npml];
            const real core_z_max = ctx.grid_spacing.z_bounds[ctx.NzT - ctx.npml];

            std::cout << "\n========================================\n";
            std::cout << "COORDINATE SYSTEM DIAGNOSTICS (Step A)\n";
            std::cout << "========================================\n";

            std::cout << "\n[Total Domain] (includes PML):\n";
            std::cout << "  x: [" << total_x_min * 1e9 << ", " << total_x_max * 1e9 << "] nm (indices 0 to " << ctx.NxT << ")\n";
            std::cout << "  y: [" << total_y_min * 1e9 << ", " << total_y_max * 1e9 << "] nm (indices 0 to " << ctx.NyT << ")\n";
            std::cout << "  z: [" << total_z_min * 1e9 << ", " << total_z_max * 1e9 << "] nm (indices 0 to " << ctx.NzT << ")\n";

            std::cout << "\n[Core Region] (physical domain, no PML):\n";
            std::cout << "  x: indices [" << ctx.npml << ", " << (ctx.NxT - ctx.npml) << ")\n";
            std::cout << "     total coords: [" << core_x_min * 1e9 << ", " << core_x_max * 1e9 << "] nm\n";
            std::cout << "     physical coords: [0, " << (core_x_max - core_x_min) * 1e9 << "] nm\n";
            std::cout << "  y: indices [" << ctx.npml << ", " << (ctx.NyT - ctx.npml) << ")\n";
            std::cout << "     total coords: [" << core_y_min * 1e9 << ", " << core_y_max * 1e9 << "] nm\n";
            std::cout << "     physical coords: [0, " << (core_y_max - core_y_min) * 1e9 << "] nm\n";
            std::cout << "  z: indices [" << ctx.npml << ", " << (ctx.NzT - ctx.npml) << ")\n";
            std::cout << "     total coords: [" << core_z_min * 1e9 << ", " << core_z_max * 1e9 << "] nm\n";
            std::cout << "     physical coords: [0, " << (core_z_max - core_z_min) * 1e9 << "] nm\n";

            std::cout << "\n[PML Offsets] (physical_coord = total_coord - offset):\n";
            std::cout << "  pml_offset_x = " << pml_offset_x * 1e9 << " nm\n";
            std::cout << "  pml_offset_y = " << pml_offset_y * 1e9 << " nm\n";
            std::cout << "  pml_offset_z = " << pml_offset_z * 1e9 << " nm\n";

            std::cout << "\n[Coordinate Conversion Formula]:\n";
            std::cout << "  x_total = x_physical + pml_offset_x\n";
            std::cout << "  x_physical = x_total - pml_offset_x\n";
            std::cout << "========================================\n";

            // Save coordinate diagnostics to file
            std::ofstream coord_file("frames/" + P.run_tag + "/coordinate_mapping.txt");
            if (coord_file) {
                coord_file << "# Coordinate System Mapping Diagnostics\n";
                coord_file << "# =====================================\n\n";
                coord_file << "# Total Domain (includes PML)\n";
                coord_file << "total_x_min_nm = " << total_x_min * 1e9 << "\n";
                coord_file << "total_x_max_nm = " << total_x_max * 1e9 << "\n";
                coord_file << "total_y_min_nm = " << total_y_min * 1e9 << "\n";
                coord_file << "total_y_max_nm = " << total_y_max * 1e9 << "\n";
                coord_file << "total_z_min_nm = " << total_z_min * 1e9 << "\n";
                coord_file << "total_z_max_nm = " << total_z_max * 1e9 << "\n\n";
                coord_file << "# Core Region indices\n";
                coord_file << "i_core_start = " << ctx.npml << "\n";
                coord_file << "i_core_end = " << (ctx.NxT - ctx.npml) << "\n";
                coord_file << "j_core_start = " << ctx.npml << "\n";
                coord_file << "j_core_end = " << (ctx.NyT - ctx.npml) << "\n";
                coord_file << "k_core_start = " << ctx.npml << "\n";
                coord_file << "k_core_end = " << (ctx.NzT - ctx.npml) << "\n\n";
                coord_file << "# Core Region in total coordinates (nm)\n";
                coord_file << "core_x_min_total_nm = " << core_x_min * 1e9 << "\n";
                coord_file << "core_x_max_total_nm = " << core_x_max * 1e9 << "\n";
                coord_file << "core_y_min_total_nm = " << core_y_min * 1e9 << "\n";
                coord_file << "core_y_max_total_nm = " << core_y_max * 1e9 << "\n";
                coord_file << "core_z_min_total_nm = " << core_z_min * 1e9 << "\n";
                coord_file << "core_z_max_total_nm = " << core_z_max * 1e9 << "\n\n";
                coord_file << "# PML Offsets (nm)\n";
                coord_file << "pml_offset_x_nm = " << pml_offset_x * 1e9 << "\n";
                coord_file << "pml_offset_y_nm = " << pml_offset_y * 1e9 << "\n";
                coord_file << "pml_offset_z_nm = " << pml_offset_z * 1e9 << "\n\n";
                coord_file << "# Physical domain size (nm)\n";
                coord_file << "phys_size_x_nm = " << (core_x_max - core_x_min) * 1e9 << "\n";
                coord_file << "phys_size_y_nm = " << (core_y_max - core_y_min) * 1e9 << "\n";
                coord_file << "phys_size_z_nm = " << (core_z_max - core_z_min) * 1e9 << "\n";
            }
        }

        // Save mesh metadata as JSON (with source position for Python visualization)
        {
            // Calculate source position - target physical coordinates
            const real phys_x_size = ctx.grid_spacing.x_bounds[ctx.NxT - ctx.npml] - ctx.grid_spacing.x_bounds[ctx.npml];
            const real phys_y_size = ctx.grid_spacing.y_bounds[ctx.NyT - ctx.npml] - ctx.grid_spacing.y_bounds[ctx.npml];
            const real phys_z_size = ctx.grid_spacing.z_bounds[ctx.NzT - ctx.npml] - ctx.grid_spacing.z_bounds[ctx.npml];

            // Use explicitly specified source positions (relative to physical domain origin)
            real src_x_target = P.source_x - P.domain_x_min;
            real src_y_target = P.source_y - P.domain_y_min;
            real src_z_target = P.source_z - P.domain_z_min;

            // Convert target physical coordinates to grid indices
            const size_t src_i = ctx.grid_spacing.physical_to_index_x(src_x_target);
            const size_t src_j = ctx.grid_spacing.physical_to_index_y(src_y_target);
            const size_t src_k = ctx.grid_spacing.physical_to_index_z(src_z_target);

            // Get ACTUAL physical coordinates at the mapped indices (Python should use this for marker)
            const real src_x_actual = ctx.grid_spacing.cell_center_physical_x(src_i);
            const real src_y_actual = ctx.grid_spacing.cell_center_physical_y(src_j);
            const real src_z_actual = ctx.grid_spacing.cell_center_physical_z(src_k);

            std::ofstream f("frames/" + P.run_tag + "/mesh_info.json");
            if (!f) {
                std::cerr << "[ERR] Failed to open mesh_info.json for writing\n";
            } else {
                f << std::setprecision(15);
                f << "{\n";
                f << "  \"mode\": \"" << (P.mesh_mode == MeshMode::AUTO_NONUNIFORM ? "AUTO_NONUNIFORM" : "UNIFORM") << "\",\n";
                f << "  \"NxT\": " << ctx.NxT << ",\n";
                f << "  \"NyT\": " << ctx.NyT << ",\n";
                f << "  \"NzT\": " << ctx.NzT << ",\n";
                f << "  \"npml\": " << ctx.npml << ",\n";
                f << "  \"dx_min_nm\": " << dx_min * 1e9 << ",\n";
                f << "  \"dx_max_nm\": " << dx_max * 1e9 << ",\n";
                f << "  \"dy_min_nm\": " << dy_min * 1e9 << ",\n";
                f << "  \"dy_max_nm\": " << dy_max * 1e9 << ",\n";
                f << "  \"dz_min_nm\": " << dz_min * 1e9 << ",\n";
                f << "  \"dz_max_nm\": " << dz_max * 1e9 << ",\n";
                f << "  \"domain_x_nm\": " << P.grid_spacing.x_bounds.back() * 1e9 << ",\n";
                f << "  \"domain_y_nm\": " << P.grid_spacing.y_bounds.back() * 1e9 << ",\n";
                f << "  \"domain_z_nm\": " << P.grid_spacing.z_bounds.back() * 1e9 << ",\n";
                f << "  \"dt_fs\": " << P.dt * 1e15 << ",\n";
                f << "  \"lambda0_nm\": " << P.lambda0 * 1e9 << ",\n";

                // Source position: TARGET physical coordinates (what user requested)
                f << "  \"source_x_target_nm\": " << src_x_target * 1e9 << ",\n";
                f << "  \"source_y_target_nm\": " << src_y_target * 1e9 << ",\n";
                f << "  \"source_z_target_nm\": " << src_z_target * 1e9 << ",\n";

                // Source position: MAPPED grid indices (in total domain)
                f << "  \"source_i\": " << src_i << ",\n";
                f << "  \"source_j\": " << src_j << ",\n";
                f << "  \"source_k\": " << src_k << ",\n";

                // Source position: ACTUAL physical coordinates at mapped indices (USE THIS FOR MARKER)
                f << "  \"source_x_actual_nm\": " << src_x_actual * 1e9 << ",\n";
                f << "  \"source_y_actual_nm\": " << src_y_actual * 1e9 << ",\n";
                f << "  \"source_z_actual_nm\": " << src_z_actual * 1e9 << ",\n";

                // Physical domain boundaries (AUTHORITATIVE - use these, not array indices)
                f << "  \"phys_domain_x_min_nm\": 0.0,\n";
                f << "  \"phys_domain_x_max_nm\": " << phys_x_size * 1e9 << ",\n";
                f << "  \"phys_domain_y_min_nm\": 0.0,\n";
                f << "  \"phys_domain_y_max_nm\": " << phys_y_size * 1e9 << ",\n";
                f << "  \"phys_domain_z_min_nm\": 0.0,\n";
                f << "  \"phys_domain_z_max_nm\": " << phys_z_size * 1e9 << ",\n";

                // Grid type info for Python extent calculation
                f << "  \"grid_type\": \"bounds\",\n";
                f << "  \"grid_x_count\": " << ctx.grid_spacing.x_bounds.size() << ",\n";
                f << "  \"grid_y_count\": " << ctx.grid_spacing.y_bounds.size() << ",\n";
                f << "  \"grid_z_count\": " << ctx.grid_spacing.z_bounds.size() << ",\n";

                f << "  \"num_structures\": " << UserConfig::STRUCTURES.size();

                // Export structure definitions for Python visualization
                if (!UserConfig::STRUCTURES.empty()) {
                    f << ",\n  \"structures\": [\n";
                    for (size_t i = 0; i < UserConfig::STRUCTURES.size(); ++i) {
                        const auto& s = UserConfig::STRUCTURES[i];
                        f << "    {\"type\": \"" << s.type << "\", \"n\": " << s.n << ", \"params_nm\": [";
                        for (int j = 0; j < 6; ++j) {
                            f << s.params[j] * 1e9;
                            if (j < 5) f << ", ";
                        }
                        f << "]}";
                        if (i < UserConfig::STRUCTURES.size() - 1) f << ",";
                        f << "\n";
                    }
                    f << "  ]\n";
                } else {
                    f << ",\n  \"structures\": []\n";
                }
                f << "}\n";
                f.close();
            }
        }

        std::cout << "========================================\n\n";

        return ctx;
    }

    // --- Bake materials and bind to boundary (applies structure packages) ---
    inline void bake_and_bind(SimContext& ctx, SimulationParams& P, int subpixel_samples = 8) {
        ctx.scene.bake(ctx.mats, P.dt, eps0, mu0, subpixel_samples);
        ctx.bc->bind_material_coeffs(&ctx.mats.bEx, &ctx.mats.bEy, &ctx.mats.bEz,
                                      &ctx.mats.bHx, &ctx.mats.bHy, &ctx.mats.bHz);
    }

    // ===== REGISTER MESH OVERRIDE REGIONS FOR AUTO MESH =====
    // Call this BEFORE make_context() when using AUTO_NONUNIFORM mode
    // Regions are read from UserConfig::MESH_OVERRIDE_REGIONS
    // NOTE: Override regions are in PHYSICAL domain coordinates (no conversion needed)
    //       The two-step mesh generator works directly in physical coordinates
    inline void register_mesh_override_regions(SimulationParams& P) {
        P.auto_mesh_config.override_regions.clear();

        if (!UserConfig::MESH_OVERRIDE_ENABLED) {
            P.auto_mesh_config.mesh_override_enabled = false;
            return;
        }

        P.auto_mesh_config.mesh_override_enabled = true;

        // Register each override region from user config
        // Coordinates are in PHYSICAL domain - no conversion needed
        for (const auto& ur : UserConfig::MESH_OVERRIDE_REGIONS) {
            if (!ur.enabled) continue;

            MeshOverrideRegion region;
            region.enabled = true;
            // Direct copy - coordinates are already in physical domain
            region.x0 = ur.x0;
            region.x1 = ur.x1;
            region.y0 = ur.y0;
            region.y1 = ur.y1;
            region.z0 = ur.z0;
            region.z1 = ur.z1;
            region.dx = ur.dx;
            region.dy = ur.dy;
            region.dz = ur.dz;
            region.priority = 0;  // Default priority

            P.auto_mesh_config.override_regions.push_back(region);
        }

        std::cout << "[AutoMesh] Pre-registered " << P.auto_mesh_config.override_regions.size()
                  << " mesh override regions from UserConfig (physical coordinates)\n";
    }

    // ===== PRE-REGISTER STRUCTURE BOUNDS FOR AUTO MESH =====
    // Call this BEFORE make_context() when using AUTO_NONUNIFORM mode
    // Structures are read from UserConfig::STRUCTURES
    // NOTE: Structure bounds are in PHYSICAL domain coordinates (no conversion needed)
    //       The two-step mesh generator works directly in physical coordinates
    inline void register_structure_bounds_for_auto_mesh(SimulationParams& P) {
        clear_structure_bounds();

        // Register each structure from user config
        // Coordinates are in PHYSICAL domain - no conversion needed
        for (const auto& s : UserConfig::STRUCTURES) {
            if (s.type == "box") {
                // params = {x0, x1, y0, y1, z0, z1}
                register_structure_bounds(
                    s.params[0], s.params[1],
                    s.params[2], s.params[3],
                    s.params[4], s.params[5],
                    s.n
                );
            }
            else if (s.type == "sphere") {
                // params = {cx, cy, cz, radius, 0, 0}
                real cx = s.params[0];
                real cy = s.params[1];
                real cz = s.params[2];
                real r  = s.params[3];
                register_structure_bounds(
                    cx - r, cx + r,
                    cy - r, cy + r,
                    cz - r, cz + r,
                    s.n
                );
            }
            else if (s.type == "cylinder") {
                // params = {cx, cy, radius, z0, z1, 0}
                real cx = s.params[0];
                real cy = s.params[1];
                real r  = s.params[2];
                real z0 = s.params[3];
                real z1 = s.params[4];
                register_structure_bounds(
                    cx - r, cx + r,
                    cy - r, cy + r,
                    z0, z1,
                    s.n
                );
            }
        }

        std::cout << "[AutoMesh] Pre-registered " << g_structure_bounds.size()
                  << " structure bounds from UserConfig (physical coordinates)\n";

        // Register mesh override regions
        register_mesh_override_regions(P);
    }

    // ===== ADD STRUCTURES TO SCENE (after mesh is generated) =====
    // Structures are read from UserConfig::STRUCTURES
    inline void register_default_structures() {
        g_structure_pkgs.push_back([](StructureScene& scene) {

            // Background material from UserConfig
            scene.bg = make_nk(UserConfig::BACKGROUND_N);

            // PML offset (structures start after PML)
            real pml_x = scene.grid_spacing.x_bounds[scene.npml_cells];
            real pml_y = scene.grid_spacing.y_bounds[scene.npml_cells];
            real pml_z = scene.grid_spacing.z_bounds[scene.npml_cells];

            // Add each structure from user config
            for (const auto& s : UserConfig::STRUCTURES) {
                Material mat = make_nk(s.n);

                if (s.type == "box") {
                    // params = {x0, x1, y0, y1, z0, z1}
                    scene.add_box(AABB{
                        pml_x + s.params[0], pml_x + s.params[1],
                        pml_y + s.params[2], pml_y + s.params[3],
                        pml_z + s.params[4], pml_z + s.params[5]
                    }, mat);
                }
                else if (s.type == "sphere") {
                    // params = {cx, cy, cz, radius, 0, 0}
                    scene.add_sphere(
                        pml_x + s.params[0],
                        pml_y + s.params[1],
                        pml_z + s.params[2],
                        s.params[3],
                        mat
                    );
                }
                else if (s.type == "cylinder") {
                    // params = {cx, cy, radius, z0, z1, 0}
                    scene.add_cylinder_z(
                        pml_x + s.params[0],
                        pml_y + s.params[1],
                        s.params[2],
                        pml_z + s.params[3],
                        pml_z + s.params[4],
                        mat
                    );
                }
            }

            std::cout << "[Scene] Added " << scene.items.size() << " structures from UserConfig\n";
        });
    }

    // Source package: registers all sources from UserConfig
    // Uses the new modular source system with DIPOLE_SOURCES and PLANE_WAVE_SOURCES
    inline void register_default_sources(const SimulationParams& P) {
        g_source_pkgs.push_back([&P](const SimContext& ctx, Runtime& rt) {
            std::cout << "\n========================================\n";
            std::cout << "REGISTERING SOURCES\n";
            std::cout << "========================================\n";

            // Helper to convert waveform int to enum
            auto to_waveform = [](int wf) -> Sources::Waveform {
                switch (wf) {
                    case 0: return Sources::Waveform::Ricker;
                    case 1: return Sources::Waveform::GaussianModulatedSine;
                    case 2: return Sources::Waveform::RickerLikeGaussian2nd;
                    case 3: return Sources::Waveform::ContinuousWave;
                    default: return Sources::Waveform::Ricker;
                }
            };

            auto to_polarization = [](int pol) -> Sources::Polarization {
                switch (pol) {
                    case 0: return Sources::Polarization::Ex;
                    case 1: return Sources::Polarization::Ey;
                    case 2: return Sources::Polarization::Ez;
                    case 3: return Sources::Polarization::Hx;
                    case 4: return Sources::Polarization::Hy;
                    case 5: return Sources::Polarization::Hz;
                    default: return Sources::Polarization::Ez;
                }
            };

            // Register dipole sources from UserConfig
            for (const auto& def : UserConfig::DIPOLE_SOURCES) {
                if (!def.enabled) continue;

                Sources::SourceConfig cfg;
                cfg.x = def.x;
                cfg.y = def.y;
                cfg.z = def.z;
                cfg.amplitude = def.amplitude;
                cfg.frequency = def.frequency;
                cfg.wavelength = def.wavelength;
                cfg.tau = def.tau;
                cfg.df_fwhm = def.df_fwhm;
                cfg.t0_factor = def.t0_factor;
                cfg.waveform = to_waveform(def.waveform);
                cfg.polarization = to_polarization(def.polarization);

                auto src = Sources::make_dipole_source(
                    cfg, ctx.grid_spacing, ctx.NyT, ctx.NzT, ctx.dt,
                    P.domain_x_min, P.domain_y_min, P.domain_z_min
                );
                if (src) rt.sources.push_back(std::move(src));
            }

            // Register plane wave sources from UserConfig
            for (const auto& def : UserConfig::PLANE_WAVE_SOURCES) {
                if (!def.enabled) continue;

                Sources::PlaneWaveConfig cfg;
                cfg.injection_position = def.injection_position;
                cfg.amplitude = def.amplitude;
                cfg.frequency = def.frequency;
                cfg.wavelength = def.wavelength;
                cfg.tau = def.tau;
                cfg.df_fwhm = def.df_fwhm;
                cfg.t0_factor = def.t0_factor;
                cfg.waveform = to_waveform(def.waveform);

                // Direction
                switch (def.direction) {
                    case 0: cfg.direction = Sources::PlaneWaveDirection::PlusX; break;
                    case 1: cfg.direction = Sources::PlaneWaveDirection::MinusX; break;
                    case 2: cfg.direction = Sources::PlaneWaveDirection::PlusY; break;
                    case 3: cfg.direction = Sources::PlaneWaveDirection::MinusY; break;
                    case 4: cfg.direction = Sources::PlaneWaveDirection::PlusZ; break;
                    case 5: cfg.direction = Sources::PlaneWaveDirection::MinusZ; break;
                    default: cfg.direction = Sources::PlaneWaveDirection::PlusZ; break;
                }
                cfg.polarization = to_polarization(def.polarization);
                cfg.x_min = def.x_min; cfg.x_max = def.x_max;
                cfg.y_min = def.y_min; cfg.y_max = def.y_max;
                cfg.z_min = def.z_min; cfg.z_max = def.z_max;

                auto src = Sources::make_plane_wave_source(
                    cfg, ctx.grid_spacing, ctx.NxT, ctx.NyT, ctx.NzT, ctx.npml, ctx.dt,
                    P.domain_x_min, P.domain_y_min, P.domain_z_min
                );
                if (src) rt.sources.push_back(std::move(src));
            }

            std::cout << "[Sources] Registered " << rt.sources.size() << " sources\n";
            std::cout << "========================================\n\n";
        });
    }

    // New modular detector bundle
    struct NewDetBundle {
        std::vector<std::unique_ptr<Detectors::MeshDetector>> mesh_detectors;
        std::vector<std::unique_ptr<Detectors::FieldMovie2D>> field_movies;
        std::vector<std::unique_ptr<Detectors::PointFieldDetector>> point_probes;
        // Legacy detectors for backward compatibility
        std::unique_ptr<BoxPoyntingDetector> boxFlux;
        std::unique_ptr<BoxEnergyDetector>   boxEM;
    };

    // Detector package: registers all detectors from UserConfig
    // Uses the new modular detector system
    inline void register_default_detectors(const SimulationParams& P) {
        g_detector_pkgs.push_back([&P](const SimContext& ctx, Runtime& rt) {
            using std::filesystem::path;
            auto bundle = std::make_shared<NewDetBundle>();

            std::cout << "\n========================================\n";
            std::cout << "REGISTERING DETECTORS\n";
            std::cout << "========================================\n";

            // Output directory
            const path out_root = path("frames") / P.run_tag;
            std::error_code ec;
            std::filesystem::create_directories(out_root, ec);

            // Get refractive index grid for mesh detectors
            std::vector<real> n_grid;
            make_n_grid_from_scene(ctx.scene, n_grid);

            // Helper to convert field component int to enum
            auto to_field_component = [](int fc) -> Detectors::FieldComponent {
                switch (fc) {
                    case 0: return Detectors::FieldComponent::Ex;
                    case 1: return Detectors::FieldComponent::Ey;
                    case 2: return Detectors::FieldComponent::Ez;
                    case 3: return Detectors::FieldComponent::Hx;
                    case 4: return Detectors::FieldComponent::Hy;
                    case 5: return Detectors::FieldComponent::Hz;
                    case 6: return Detectors::FieldComponent::E_magnitude;
                    case 7: return Detectors::FieldComponent::H_magnitude;
                    case 8: return Detectors::FieldComponent::Sx;
                    case 9: return Detectors::FieldComponent::Sy;
                    case 10: return Detectors::FieldComponent::Sz;
                    case 11: return Detectors::FieldComponent::S_magnitude;
                    default: return Detectors::FieldComponent::Ez;
                }
            };

            // Register mesh detectors from UserConfig
            for (const auto& def : UserConfig::MESH_DETECTORS) {
                if (!def.enabled) continue;

                Detectors::MeshDetectorConfig cfg;
                cfg.name = def.name;
                cfg.slice_position = def.slice_position;
                cfg.slice_plane = def.slice_plane;
                cfg.export_mesh_lines = def.export_mesh_lines;
                cfg.export_spacing_arrays = def.export_spacing;

                auto det = Detectors::make_mesh_detector(
                    out_root, cfg,
                    ctx.NxT, ctx.NyT, ctx.NzT,
                    ctx.Nx_core, ctx.Ny_core, ctx.Nz_core,
                    ctx.npml, ctx.grid_spacing,
                    P.domain_size_x(), P.domain_size_y(), P.domain_size_z(),
                    P.domain_x_min, P.domain_y_min, P.domain_z_min
                );
                if (det) {
                    det->save_index_slice(n_grid);
                    det->save_mesh_info();
                    bundle->mesh_detectors.push_back(std::move(det));
                }
            }

            // Register field movie 2D detectors from UserConfig
            for (const auto& def : UserConfig::FIELD_MOVIE_2D_DETECTORS) {
                if (!def.enabled) continue;

                Detectors::FieldMovie2DConfig cfg;
                cfg.name = def.name;
                cfg.component = to_field_component(def.field_component);
                cfg.slice_plane = def.slice_plane;
                cfg.slice_position = def.slice_position;
                cfg.save_every = P.saveEvery;
                cfg.n_steps = P.nSteps;
                cfg.write_float64 = P.writeFloat64;
                cfg.frame_pattern = def.frame_pattern;

                auto det = Detectors::make_field_movie_2d(
                    out_root, cfg,
                    ctx.NxT, ctx.NyT, ctx.NzT,
                    ctx.Nx_core, ctx.Ny_core, ctx.Nz_core,
                    ctx.npml, ctx.grid_spacing,
                    P.domain_size_x(), P.domain_size_y(), P.domain_size_z(),
                    P.domain_x_min, P.domain_y_min, P.domain_z_min
                );
                if (det) bundle->field_movies.push_back(std::move(det));
            }

            // Register point field detectors from UserConfig
            for (const auto& def : UserConfig::POINT_FIELD_DETECTORS) {
                if (!def.enabled) continue;

                Detectors::PointFieldDetectorConfig cfg;
                cfg.name = def.name;
                cfg.x = def.x;
                cfg.y = def.y;
                cfg.z = def.z;
                cfg.save_every = P.saveEvery;
                cfg.n_steps = P.nSteps;
                cfg.write_float64 = P.writeFloat64;
                for (int c : def.components) {
                    cfg.components.push_back(to_field_component(c));
                }

                auto det = Detectors::make_point_field_detector(
                    out_root, cfg,
                    ctx.NxT, ctx.NyT, ctx.NzT,
                    ctx.npml, ctx.grid_spacing, ctx.dt,
                    P.domain_x_min, P.domain_y_min, P.domain_z_min
                );
                if (det) bundle->point_probes.push_back(std::move(det));
            }

            // Legacy detectors: Box power and energy (always enabled)
            const size_t half = 49;
            const size_t ci = ctx.NxT / 2, cj = ctx.NyT / 2, ck = ctx.NzT / 2;
            size_t i0 = (ci > half ? ci - half : 1);
            size_t j0 = (cj > half ? cj - half : 1);
            size_t k0 = (ck > half ? ck - half : 1);
            size_t i1 = std::min(ctx.NxT - 2, ci + half - 1);
            size_t j1 = std::min(ctx.NyT - 2, cj + half - 1);
            size_t k1 = std::min(ctx.NzT - 2, ck + half - 1);

            bundle->boxFlux = std::make_unique<BoxPoyntingDetector>(
                out_root, "box_flux80",
                ctx.NxT, ctx.NyT, ctx.NzT,
                i0, i1, j0, j1, k0, k1,
                P.saveEvery, P.nSteps,
                ctx.grid_spacing, ctx.dt
            );

            bundle->boxEM = std::make_unique<BoxEnergyDetector>(
                out_root,
                ctx.NxT, ctx.NyT, ctx.NzT,
                i0, i1, j0, j1, k0, k1,
                P.saveEvery, P.nSteps,
                ctx.grid_spacing, ctx.dt,
                ctx.mats.bEx, ctx.mats.bEy, ctx.mats.bEz,
                ctx.mats.bHx, ctx.mats.bHy, ctx.mats.bHz
            );

            std::cout << "[Detectors] Registered:\n";
            std::cout << "  - " << bundle->mesh_detectors.size() << " mesh detectors\n";
            std::cout << "  - " << bundle->field_movies.size() << " field movie 2D detectors\n";
            std::cout << "  - " << bundle->point_probes.size() << " point field detectors\n";
            std::cout << "  - 2 legacy detectors (box flux/energy)\n";
            std::cout << "========================================\n\n";

            // Unified callback chaining
            auto prev_after_E = rt.after_E;
            rt.after_E = [bundle, dt = ctx.dt, prev_after_E](
                size_t n,
                const std::vector<real>& Ex,
                const std::vector<real>& Ey,
                const std::vector<real>& Ez,
                const std::vector<real>& Hx,
                const std::vector<real>& Hy,
                const std::vector<real>& Hz) {
                    prev_after_E(n, Ex, Ey, Ez, Hx, Hy, Hz);

                    // Legacy detectors
                    bundle->boxFlux->try_record(n, Ex, Ey, Ez, Hx, Hy, Hz);
                    bundle->boxEM->try_record(n, Ex, Ey, Ez, Hx, Hy, Hz);

                    // New modular detectors
                    for (auto& det : bundle->field_movies) {
                        det->record_after_E(n, dt, Ex, Ey, Ez, Hx, Hy, Hz);
                    }
                    for (auto& det : bundle->point_probes) {
                        det->record_after_E(n, dt, Ex, Ey, Ez, Hx, Hy, Hz);
                    }
                };
        });
    }

    // Register all default packages: structures, sources, and detectors
    inline void register_all_default_packages(SimulationParams& P) {
        // For AUTO_NONUNIFORM mode, pre-register structure bounds first
        // This must be called BEFORE make_context()
        if (P.mesh_mode == MeshMode::AUTO_NONUNIFORM) {
            register_structure_bounds_for_auto_mesh(P);
        }

        register_default_structures();
        register_default_sources(P);
        register_default_detectors(P);
    }

    // Convenience function for complete setup with AUTO_NONUNIFORM mesh
    // Usage:
    //   Config::sim.mesh_mode = MeshMode::AUTO_NONUNIFORM;
    //   Config::register_all_default_packages(Config::sim);
    //   auto ctx = Config::make_context(Config::sim);
    //   Config::apply_structure_packages(ctx);
    //   Config::bake_and_bind(ctx, Config::sim);

    // -------------------------- Global unique simulation instance --------------------------
    // Initialized from user_config.hpp values
    inline SimulationParams sim = [] {
        SimulationParams p;

        // ===== Mesh mode selection (from UserConfig) =====
        // 0 = UNIFORM, 1 = AUTO_NONUNIFORM (ppw-based, interface snapping)
        switch (UserConfig::MESH_MODE) {
            case 0: p.mesh_mode = MeshMode::UNIFORM; break;
            case 1: default: p.mesh_mode = MeshMode::AUTO_NONUNIFORM; break;
        }

        // Uniform mesh settings
        p.dx_uniform = UserConfig::DX_UNIFORM;
        p.dy_uniform = UserConfig::DY_UNIFORM;
        p.dz_uniform = UserConfig::DZ_UNIFORM;

        // Grid dimensions
        p.Nx = UserConfig::NX;
        p.Ny = UserConfig::NY;
        p.Nz = UserConfig::NZ;

        // Wavelength and domain settings
        p.lambda0 = UserConfig::LAMBDA0;
        p.domain_x_min = UserConfig::DOMAIN_X_MIN;
        p.domain_x_max = UserConfig::DOMAIN_X_MAX;
        p.domain_y_min = UserConfig::DOMAIN_Y_MIN;
        p.domain_y_max = UserConfig::DOMAIN_Y_MAX;
        p.domain_z_min = UserConfig::DOMAIN_Z_MIN;
        p.domain_z_max = UserConfig::DOMAIN_Z_MAX;

        // Time stepping
        p.S = UserConfig::CFL_FACTOR;
        p.nSteps = UserConfig::N_STEPS;
        p.saveEvery = UserConfig::SAVE_EVERY;

        // Source parameters
        p.Iz_peak = UserConfig::SOURCE_IZ_PEAK;
        p.f0 = UserConfig::SOURCE_F0;
        p.tau_src = UserConfig::SOURCE_TAU;
        p.df_fwhm = UserConfig::SOURCE_DF_FWHM;
        p.t0_factor = UserConfig::SOURCE_T0_FACTOR;
        // Source physical coordinates (absolute)
        p.source_x = UserConfig::SOURCE_X;
        p.source_y = UserConfig::SOURCE_Y;
        p.source_z = UserConfig::SOURCE_Z;

        // Detector settings
        p.run_tag = UserConfig::RUN_TAG;
        p.z_slice_z = UserConfig::Z_SLICE_Z;
        p.framePattern = UserConfig::FRAME_PATTERN;
        p.writeFloat64 = UserConfig::WRITE_FLOAT64;
        p.probe_x = UserConfig::PROBE_X;
        p.probe_y = UserConfig::PROBE_Y;
        p.probe_z = UserConfig::PROBE_Z;

        // Boundary conditions
        p.bp.type = (UserConfig::BOUNDARY_TYPE == 0) ? BcType::PEC : BcType::CPML_RC;
        p.bp.cpml.npml = UserConfig::CPML_NPML;
        p.bp.cpml.m = UserConfig::CPML_M;
        p.bp.cpml.Rerr = UserConfig::CPML_RERR;
        p.bp.cpml.alpha0 = UserConfig::CPML_ALPHA0;
        p.bp.cpml.kappa_max = UserConfig::CPML_KAPPA_MAX;
        p.bp.cpml.alpha_linear = UserConfig::CPML_ALPHA_LINEAR;

        // Auto mesh configuration (ppw-based, interface snapping)
        p.auto_mesh_config.lambda_min = p.lambda0;
        p.auto_mesh_config.mesh_accuracy = UserConfig::AUTO_MESH_ACCURACY;
        p.auto_mesh_config.ppw_override = UserConfig::AUTO_MESH_PPW_OVERRIDE;
        p.auto_mesh_config.dx_min = UserConfig::AUTO_MESH_DX_MIN;
        p.auto_mesh_config.dx_max = UserConfig::AUTO_MESH_DX_MAX;
        p.auto_mesh_config.max_grading_ratio = UserConfig::AUTO_MESH_MAX_GRADING_RATIO;
        p.auto_mesh_config.mesh_override_enabled = UserConfig::MESH_OVERRIDE_ENABLED;
        // Note: Override regions will be populated with correct PML offsets in register_mesh_override_regions()

        return p;
    }();

} // namespace Config