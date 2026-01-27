// params.hpp — Simulation configuration and parameter management

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
#include "auto_mesh_generator.hpp"

#include "sources/isource.hpp"
#include "sources/dipole_source.hpp"
#include "sources/plane_wave_source.hpp"

#include "detectors/idetector.hpp"
#include "detectors/mesh_detector.hpp"
#include "detectors/field_movie_2d.hpp"
#include "detectors/point_field_detector.hpp"

namespace Config {

    inline constexpr real PI   = PhysConst::PI;
    inline constexpr real c0   = PhysConst::C0;
    inline constexpr real eps0 = PhysConst::EPS0;
    inline constexpr real mu0  = PhysConst::MU0;

    enum class MeshMode { UNIFORM, AUTO_NONUNIFORM };

    struct SimulationParams {
        MeshMode mesh_mode = MeshMode::AUTO_NONUNIFORM;

        real dx_uniform = 10e-9;
        real dy_uniform = 10e-9;
        real dz_uniform = 10e-9;

        AutoMeshConfig auto_mesh_config;
        GridSpacing grid_spacing;

        real domain_x_min = 0, domain_x_max = 0;
        real domain_y_min = 0, domain_y_max = 0;
        real domain_z_min = 0, domain_z_max = 0;

        real domain_size_x() const { return domain_x_max - domain_x_min; }
        real domain_size_y() const { return domain_y_max - domain_y_min; }
        real domain_size_z() const { return domain_z_max - domain_z_min; }

        size_t Nx = 100, Ny = 100, Nz = 100;
        real S = 0.9;  // CFL factor

        size_t nSteps = 1000;
        size_t saveEvery = 1;

        real dt = 0.0;
        real Ce = 0.0;
        real Ch = 0.0;

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

        std::string run_tag = "3D_FDTD_v008_output";

        real Iz_peak = 1e-5;
        real f0 = 0.0;
        real lambda0 = 500e-9;
        real tau_src = 5e-15;
        real df_fwhm = 2.4e13;
        real t0_factor = 0.4;

        real tau_eff = 0.0;
        real omega0 = 0.0;
        real t0 = 0.0;

        real z_slice_z = 0;
        std::string framePattern = "ez_%04d.raw";
        bool writeFloat64 = true;

        real probe_x = 0, probe_y = 0, probe_z = 0;
        real source_x = 0, source_y = 0, source_z = 0;

        std::function<void(StructureScene&)> build_scene;

        void finalize() {
            real dx_est = dx_uniform;
            real dy_est = dy_uniform;
            real dz_est = dz_uniform;
            
            const real sum_inv2 = 1.0/(dx_est*dx_est) + 1.0/(dy_est*dy_est) + 1.0/(dz_est*dz_est);
            dt = S / (c0 * std::sqrt(sum_inv2));
            Ce = dt / eps0;
            Ch = dt / mu0;

            if (f0 <= 0.0) {
                f0 = (lambda0 > 0.0) ? c0 / lambda0 : c0 / 500e-9;
            }
            omega0 = 2.0 * std::numbers::pi_v<double> * f0;

            if (tau_src > 0.0) tau_eff = tau_src;
            else if (df_fwhm > 0.0) tau_eff = tau_from_df_fwhm(df_fwhm);
            else tau_eff = 20.0 * dt;

            t0 = (t0_factor > 0.0) ? t0_factor * tau_eff : 3.0 / f0;
        }
    };

    struct SimContext {
        size_t NxT{}, NyT{}, NzT{};
        size_t Nx_core{}, Ny_core{}, Nz_core{};
        size_t npml{};
        GridSpacing grid_spacing;
        real dt{};
        std::unique_ptr<IBoundary> bc;
        StructureScene scene;
        MaterialGrids mats;
        std::filesystem::path out_root;
    };

    struct Runtime {
        std::vector<std::unique_ptr<Sources::ISource>> sources;

        std::function<void(size_t,
            const std::vector<real>&, const std::vector<real>&, const std::vector<real>&,
            const std::vector<real>&, const std::vector<real>&, const std::vector<real>&)>
            after_H = [](auto, const auto&, const auto&, const auto&, const auto&, const auto&, const auto&) {};

        std::function<void(size_t,
            const std::vector<real>&, const std::vector<real>&, const std::vector<real>&,
            const std::vector<real>&, const std::vector<real>&, const std::vector<real>&)>
            after_E = [](auto, const auto&, const auto&, const auto&, const auto&, const auto&, const auto&) {};
    };

    using StructurePackage = std::function<void(StructureScene& scene)>;
    using SourcePackage = std::function<void(const SimContext& ctx, Runtime& rt)>;
    using DetectorPackage = std::function<void(const SimContext& ctx, Runtime& rt)>;

    inline std::vector<StructurePackage> g_structure_pkgs{};
    inline std::vector<SourcePackage>    g_source_pkgs{};
    inline std::vector<DetectorPackage>  g_detector_pkgs{};
    inline std::vector<StructureBounds> g_structure_bounds{};

    inline void register_structure_bounds(real x_min, real x_max,
                                          real y_min, real y_max,
                                          real z_min, real z_max,
                                          real n_material) {
        g_structure_bounds.push_back({x_min, x_max, y_min, y_max, z_min, z_max, n_material});
    }

    inline void clear_structure_bounds() { g_structure_bounds.clear(); }

    inline void apply_structure_packages(SimContext& ctx) {
        for (auto& p : g_structure_pkgs) p(ctx.scene);
    }
    inline void apply_source_packages(const SimContext& ctx, Runtime& rt) {
        for (auto& p : g_source_pkgs) p(ctx, rt);
    }
    inline void apply_detector_packages(const SimContext& ctx, Runtime& rt) {
        for (auto& p : g_detector_pkgs) p(ctx, rt);
    }

    inline SimContext make_context(SimulationParams& P) {
        SimContext ctx;
        std::cout << "\n========================================\n";
        std::cout << "Building Simulation Context\n";
        std::cout << "========================================\n";

        const size_t npml = P.bp.cpml.npml;

        switch (P.mesh_mode) {
            case MeshMode::AUTO_NONUNIFORM: {
                std::cout << "[Grid] Mode: AUTO_NONUNIFORM\n";

                P.auto_mesh_config.lambda_min = P.lambda0;
                if (P.auto_mesh_config.dx_min <= 0) P.auto_mesh_config.dx_min = P.lambda0 / 100.0;
                if (P.auto_mesh_config.dx_max <= 0) P.auto_mesh_config.dx_max = P.lambda0 / 4.0;

                real domain_x = P.domain_size_x();
                real domain_y = P.domain_size_y();
                real domain_z = P.domain_size_z();

                if (domain_x <= 0) domain_x = P.Nx * P.dx_uniform;
                if (domain_y <= 0) domain_y = P.Ny * P.dy_uniform;
                if (domain_z <= 0) domain_z = P.Nz * P.dz_uniform;

                std::cout << "[Grid] Physical domain: " << domain_x*1e9 << " x "
                          << domain_y*1e9 << " x " << domain_z*1e9 << " nm\n";
                std::cout << "[Grid] npml: " << npml << ", dx_max: " << P.auto_mesh_config.dx_max*1e9 << " nm\n";
                std::cout << "[Grid] Registered structure bounds: " << g_structure_bounds.size() << "\n";

                auto auto_mesh_structures = convert_to_auto_mesh_bounds(g_structure_bounds);
                AutoMeshGenerator generator(P.auto_mesh_config);
                P.grid_spacing = generator.generate(domain_x, domain_y, domain_z, npml, auto_mesh_structures, 1.0);

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
                std::cout << "  dx=" << P.dx_uniform*1e9 << " dy=" << P.dy_uniform*1e9 << " dz=" << P.dz_uniform*1e9 << " nm\n";

                P.grid_spacing.dx.assign(NxT, P.dx_uniform);
                P.grid_spacing.dy.assign(NyT, P.dy_uniform);
                P.grid_spacing.dz.assign(NzT, P.dz_uniform);
                build_cumulative_bounds(P.grid_spacing.dx, P.grid_spacing.x_bounds);
                build_cumulative_bounds(P.grid_spacing.dy, P.grid_spacing.y_bounds);
                build_cumulative_bounds(P.grid_spacing.dz, P.grid_spacing.z_bounds);
                break;
            }
        }

        P.grid_spacing.compute_inverses();

        const size_t NxT = P.grid_spacing.dx.size();
        const size_t NyT = P.grid_spacing.dy.size();
        const size_t NzT = P.grid_spacing.dz.size();

        real dx_min = *std::min_element(P.grid_spacing.dx.begin(), P.grid_spacing.dx.end());
        real dy_min = *std::min_element(P.grid_spacing.dy.begin(), P.grid_spacing.dy.end());
        real dz_min = *std::min_element(P.grid_spacing.dz.begin(), P.grid_spacing.dz.end());
        real dx_max = *std::max_element(P.grid_spacing.dx.begin(), P.grid_spacing.dx.end());
        real dy_max = *std::max_element(P.grid_spacing.dy.begin(), P.grid_spacing.dy.end());
        real dz_max = *std::max_element(P.grid_spacing.dz.begin(), P.grid_spacing.dz.end());

        std::cout << "\n[Grid] " << NxT << "x" << NyT << "x" << NzT << " (core " << (NxT-2*npml) << "x" << (NyT-2*npml) << "x" << (NzT-2*npml) << ")\n";
        std::cout << "[Grid] dx:" << dx_min*1e9 << "-" << dx_max*1e9 << " dy:" << dy_min*1e9 << "-" << dy_max*1e9 << " dz:" << dz_min*1e9 << "-" << dz_max*1e9 << " nm\n";

        auto check_grading = [](const std::vector<real>& spacing, const char* name) {
            real max_ratio = 1.0;
            for (size_t i = 1; i < spacing.size(); ++i) {
                real ratio = std::max(spacing[i]/spacing[i-1], spacing[i-1]/spacing[i]);
                max_ratio = std::max(max_ratio, ratio);
            }
            std::cout << "  " << name << " grading: " << max_ratio << (max_ratio > 1.2 ? " (WARNING)" : "") << "\n";
        };
        check_grading(P.grid_spacing.dx, "dx");
        check_grading(P.grid_spacing.dy, "dy");
        check_grading(P.grid_spacing.dz, "dz");

        const real sum_inv2 = (1.0/(dx_min*dx_min)) + (1.0/(dy_min*dy_min)) + (1.0/(dz_min*dz_min));
        P.dt = P.S / (c0 * std::sqrt(sum_inv2));
        P.Ce = P.dt / eps0;
        P.Ch = P.dt / mu0;
        std::cout << "[Time] dt=" << P.dt*1e15 << " fs, S=" << P.S << "\n";

        P.Nx = NxT - 2 * npml;
        P.Ny = NyT - 2 * npml;
        P.Nz = NzT - 2 * npml;

        auto bc = make_boundary(P.bp, P.Nx, P.Ny, P.Nz, P.grid_spacing, P.dt, eps0, mu0, c0);
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
        ctx.grid_spacing.npml = ctx.npml;

        ctx.scene.NxT = ctx.NxT;
        ctx.scene.NyT = ctx.NyT;
        ctx.scene.NzT = ctx.NzT;
        ctx.scene.Nx_core = ctx.Nx_core;
        ctx.scene.Ny_core = ctx.Ny_core;
        ctx.scene.Nz_core = ctx.Nz_core;
        ctx.scene.npml_cells = ctx.npml;
        ctx.scene.grid_spacing = ctx.grid_spacing;

        // Save grid spacing
        {
            namespace fs = std::filesystem;
            fs::create_directories("frames/" + P.run_tag);
            std::ofstream f("frames/" + P.run_tag + "/grid_spacing.txt");
            f << "# Mode: " << (P.mesh_mode == MeshMode::AUTO_NONUNIFORM ? "AUTO_NONUNIFORM" : "UNIFORM") << "\n";
            f << "# dx\n"; for (auto dx : ctx.grid_spacing.dx) f << dx << "\n";
            f << "# dy\n"; for (auto dy : ctx.grid_spacing.dy) f << dy << "\n";
            f << "# dz\n"; for (auto dz : ctx.grid_spacing.dz) f << dz << "\n";
            f << "# x_bounds\n"; for (auto x : ctx.grid_spacing.x_bounds) f << x << "\n";
            f << "# y_bounds\n"; for (auto y : ctx.grid_spacing.y_bounds) f << y << "\n";
            f << "# z_bounds\n"; for (auto z : ctx.grid_spacing.z_bounds) f << z << "\n";
        }

        // Coordinate diagnostics
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
            // Calculate physical domain size
            const real phys_x_size = ctx.grid_spacing.x_bounds[ctx.NxT - ctx.npml] - ctx.grid_spacing.x_bounds[ctx.npml];
            const real phys_y_size = ctx.grid_spacing.y_bounds[ctx.NyT - ctx.npml] - ctx.grid_spacing.y_bounds[ctx.npml];
            const real phys_z_size = ctx.grid_spacing.z_bounds[ctx.NzT - ctx.npml] - ctx.grid_spacing.z_bounds[ctx.npml];

            // Get source position from DIPOLE_SOURCES (use first enabled source, or domain center if none)
            real src_x_target = phys_x_size / 2.0;  // Default: domain center
            real src_y_target = phys_y_size / 2.0;
            real src_z_target = phys_z_size / 2.0;
            bool has_source = false;
            for (const auto& def : UserConfig::DIPOLE_SOURCES) {
                if (def.enabled) {
                    src_x_target = def.x - P.domain_x_min;
                    src_y_target = def.y - P.domain_y_min;
                    src_z_target = def.z - P.domain_z_min;
                    has_source = true;
                    break;  // Use first enabled source
                }
            }

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

    inline void bake_and_bind(SimContext& ctx, SimulationParams& P, int subpixel_samples = 8) {
        ctx.scene.bake(ctx.mats, P.dt, eps0, mu0, subpixel_samples);
        ctx.bc->bind_material_coeffs(&ctx.mats.bEx, &ctx.mats.bEy, &ctx.mats.bEz,
                                      &ctx.mats.bHx, &ctx.mats.bHy, &ctx.mats.bHz);
    }

    inline void register_mesh_override_regions(SimulationParams& P) {
        P.auto_mesh_config.override_regions.clear();
        if (!UserConfig::MESH_OVERRIDE_ENABLED) {
            P.auto_mesh_config.mesh_override_enabled = false;
            return;
        }
        P.auto_mesh_config.mesh_override_enabled = true;

        for (const auto& ur : UserConfig::MESH_OVERRIDE_REGIONS) {
            if (!ur.enabled) continue;
            MeshOverrideRegion region;
            region.enabled = true;
            region.x0 = ur.x0; region.x1 = ur.x1;
            region.y0 = ur.y0; region.y1 = ur.y1;
            region.z0 = ur.z0; region.z1 = ur.z1;
            region.dx = ur.dx; region.dy = ur.dy; region.dz = ur.dz;
            region.priority = 0;
            P.auto_mesh_config.override_regions.push_back(region);
        }
        std::cout << "[AutoMesh] Registered " << P.auto_mesh_config.override_regions.size() << " override regions\n";
    }

    inline void register_structure_bounds_for_auto_mesh(SimulationParams& P) {
        clear_structure_bounds();
        for (const auto& s : UserConfig::STRUCTURES) {
            if (s.type == "box") {
                register_structure_bounds(s.params[0], s.params[1], s.params[2], s.params[3], s.params[4], s.params[5], s.n);
            } else if (s.type == "sphere") {
                real cx = s.params[0], cy = s.params[1], cz = s.params[2], r = s.params[3];
                register_structure_bounds(cx-r, cx+r, cy-r, cy+r, cz-r, cz+r, s.n);
            } else if (s.type == "cylinder") {
                real cx = s.params[0], cy = s.params[1], r = s.params[2], z0 = s.params[3], z1 = s.params[4];
                register_structure_bounds(cx-r, cx+r, cy-r, cy+r, z0, z1, s.n);
            }
        }
        std::cout << "[AutoMesh] Registered " << g_structure_bounds.size() << " structure bounds\n";
        register_mesh_override_regions(P);
    }

    inline void register_default_structures() {
        g_structure_pkgs.push_back([](StructureScene& scene) {
            scene.bg = make_nk(UserConfig::BACKGROUND_N);
            real pml_x = scene.grid_spacing.x_bounds[scene.npml_cells];
            real pml_y = scene.grid_spacing.y_bounds[scene.npml_cells];
            real pml_z = scene.grid_spacing.z_bounds[scene.npml_cells];
            size_t tls_count = 0;

            for (const auto& s : UserConfig::STRUCTURES) {
                Material mat = make_nk(s.n);
                StructureTLSConfig tls_cfg;
                const auto* tls_mat = TLSMaterials::find(s.tls_material);
                if (tls_mat && UserConfig::TLS_ENABLED) {
                    tls_cfg.enabled = true;
                    tls_cfg.lambda0 = tls_mat->lambda0;
                    tls_cfg.gamma = tls_mat->gamma;
                    tls_cfg.tau = tls_mat->tau;
                    tls_cfg.N0 = tls_mat->N0;
                    tls_cfg.inversion_fraction = tls_mat->inversion_fraction;
                    tls_count++;
                    std::cout << "[Structure] Using TLS: " << tls_mat->name << "\n";
                } else if (!s.tls_material.empty() && !tls_mat) {
                    std::cerr << "[WARNING] TLS material \"" << s.tls_material << "\" not found!\n";
                }

                if (s.type == "box") {
                    scene.add_box(AABB{pml_x+s.params[0], pml_x+s.params[1], pml_y+s.params[2], pml_y+s.params[3], pml_z+s.params[4], pml_z+s.params[5]}, mat, tls_cfg);
                } else if (s.type == "sphere") {
                    scene.add_sphere(pml_x+s.params[0], pml_y+s.params[1], pml_z+s.params[2], s.params[3], mat, tls_cfg);
                } else if (s.type == "cylinder") {
                    scene.add_cylinder_z(pml_x+s.params[0], pml_y+s.params[1], s.params[2], pml_z+s.params[3], pml_z+s.params[4], mat, tls_cfg);
                }
            }
            std::cout << "[Scene] Added " << scene.items.size() << " structures" << (tls_count > 0 ? " (" + std::to_string(tls_count) + " with TLS)" : "") << "\n";
        });
    }

    inline void register_default_sources(const SimulationParams& P) {
        g_source_pkgs.push_back([&P](const SimContext& ctx, Runtime& rt) {
            std::cout << "\n[SOURCES]\n";

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
        });
    }

    struct DetectorBundle {
        std::vector<std::unique_ptr<Detectors::MeshDetector>> mesh_detectors;
        std::vector<std::unique_ptr<Detectors::FieldMovie2D>> field_movies;
        std::vector<std::unique_ptr<Detectors::PointFieldDetector>> point_probes;
    };

    inline void register_default_detectors(const SimulationParams& P) {
        g_detector_pkgs.push_back([&P](const SimContext& ctx, Runtime& rt) {
            using std::filesystem::path;
            auto bundle = std::make_shared<DetectorBundle>();
            std::cout << "\n[DETECTORS]\n";

            const path out_root = path("frames") / P.run_tag;
            std::error_code ec;
            std::filesystem::create_directories(out_root, ec);

            std::vector<real> n_grid;
            make_n_grid_from_scene(ctx.scene, n_grid);

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

            std::cout << "[Detectors] " << bundle->mesh_detectors.size() << " mesh, "
                      << bundle->field_movies.size() << " movie, " << bundle->point_probes.size() << " probe\n";

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

                    for (auto& det : bundle->field_movies) {
                        det->record_after_E(n, dt, Ex, Ey, Ez, Hx, Hy, Hz);
                    }
                    for (auto& det : bundle->point_probes) {
                        det->record_after_E(n, dt, Ex, Ey, Ez, Hx, Hy, Hz);
                    }
                };
        });
    }

    inline void register_all_default_packages(SimulationParams& P) {
        if (P.mesh_mode == MeshMode::AUTO_NONUNIFORM) {
            register_structure_bounds_for_auto_mesh(P);
        }
        register_default_structures();
        register_default_sources(P);
        register_default_detectors(P);
    }

    inline SimulationParams sim = [] {
        SimulationParams p;
        switch (UserConfig::MESH_MODE) {
            case 0: p.mesh_mode = MeshMode::UNIFORM; break;
            case 1: default: p.mesh_mode = MeshMode::AUTO_NONUNIFORM; break;
        }

        p.dx_uniform = UserConfig::DX_UNIFORM;
        p.dy_uniform = UserConfig::DY_UNIFORM;
        p.dz_uniform = UserConfig::DZ_UNIFORM;
        p.Nx = UserConfig::NX;
        p.Ny = UserConfig::NY;
        p.Nz = UserConfig::NZ;
        p.lambda0 = UserConfig::LAMBDA0;
        p.domain_x_min = UserConfig::DOMAIN_X_MIN;
        p.domain_x_max = UserConfig::DOMAIN_X_MAX;
        p.domain_y_min = UserConfig::DOMAIN_Y_MIN;
        p.domain_y_max = UserConfig::DOMAIN_Y_MAX;
        p.domain_z_min = UserConfig::DOMAIN_Z_MIN;
        p.domain_z_max = UserConfig::DOMAIN_Z_MAX;
        p.S = UserConfig::CFL_FACTOR;
        p.nSteps = UserConfig::N_STEPS;
        p.saveEvery = UserConfig::SAVE_EVERY;
        p.run_tag = UserConfig::RUN_TAG;
        p.z_slice_z = UserConfig::Z_SLICE_Z;
        p.framePattern = UserConfig::FRAME_PATTERN;
        p.writeFloat64 = UserConfig::WRITE_FLOAT64;
        p.probe_x = UserConfig::PROBE_X;
        p.probe_y = UserConfig::PROBE_Y;
        p.probe_z = UserConfig::PROBE_Z;
        p.bp.type = (UserConfig::BOUNDARY_TYPE == 0) ? BcType::PEC : BcType::CPML_RC;
        p.bp.cpml.npml = UserConfig::CPML_NPML;
        p.bp.cpml.m = UserConfig::CPML_M;
        p.bp.cpml.Rerr = UserConfig::CPML_RERR;
        p.bp.cpml.alpha0 = UserConfig::CPML_ALPHA0;
        p.bp.cpml.kappa_max = UserConfig::CPML_KAPPA_MAX;
        p.bp.cpml.alpha_linear = UserConfig::CPML_ALPHA_LINEAR;
        p.auto_mesh_config.lambda_min = p.lambda0;
        p.auto_mesh_config.mesh_accuracy = UserConfig::AUTO_MESH_ACCURACY;
        p.auto_mesh_config.ppw_override = UserConfig::AUTO_MESH_PPW_OVERRIDE;
        p.auto_mesh_config.dx_min = UserConfig::AUTO_MESH_DX_MIN;
        p.auto_mesh_config.dx_max = UserConfig::AUTO_MESH_DX_MAX;
        p.auto_mesh_config.max_grading_ratio = UserConfig::AUTO_MESH_MAX_GRADING_RATIO;
        p.auto_mesh_config.mesh_override_enabled = UserConfig::MESH_OVERRIDE_ENABLED;
        return p;
    }();

} // namespace Config