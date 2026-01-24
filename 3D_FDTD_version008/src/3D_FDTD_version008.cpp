// 3D_FDTD_version008.cpp — Main program
//
// This program simulates 3D electromagnetic wave propagation using FDTD
// with optional two-level atomic gain medium coupling.
//
// Configuration: Modify user_config.hpp to change simulation parameters.

#include "user_config.hpp"          // User-configurable parameters (MODIFY THIS FILE)
#include "params.hpp"               // Simulation parameter management (includes sources/detectors)
#include "fdtd_stepper.hpp"         // FDTD update equations
#include "boundary.hpp"             // Boundary conditions (PML)
#include "structure_material.hpp"   // Materials and geometry
#include "global_function.hpp"      // Utility functions and constants
#include "omp_config.hpp"           // OpenMP configuration
#include "two_level_system.hpp"     // Two-level gain medium
#include "two_level_detectors.hpp"  // TLS diagnostics
#include "stability_monitor.hpp"    // Numerical stability monitoring
#include "parallel_domain_v2.hpp"   // Parallel domain decomposition (V2 - true domain decomposition)

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <filesystem>
#include <chrono>

namespace fs = std::filesystem;

int main() {
    std::cout << "========================================\n";
    std::cout << "3D FDTD Electromagnetic Simulation\n";
    std::cout << "========================================\n\n";

#if FDTD_OMP_ENABLED
    std::cout << "OpenMP enabled, max threads = " << omp_get_max_threads() << "\n";
#else
    std::cout << "OpenMP NOT enabled\n";
#endif

    // ===== 1. Initialize simulation parameters from UserConfig =====
    auto& P = Config::sim;
    P.finalize();

    std::cout << "\nSimulation parameters (from user_config.hpp):\n";
    std::cout << "  Domain X: [" << UnitConv::m_to_nm(UserConfig::DOMAIN_X_MIN) << ", "
              << UnitConv::m_to_nm(UserConfig::DOMAIN_X_MAX) << "] nm\n";
    std::cout << "  Domain Y: [" << UnitConv::m_to_nm(UserConfig::DOMAIN_Y_MIN) << ", "
              << UnitConv::m_to_nm(UserConfig::DOMAIN_Y_MAX) << "] nm\n";
    std::cout << "  Domain Z: [" << UnitConv::m_to_nm(UserConfig::DOMAIN_Z_MIN) << ", "
              << UnitConv::m_to_nm(UserConfig::DOMAIN_Z_MAX) << "] nm\n";
    std::cout << "  Wavelength: " << UnitConv::m_to_nm(P.lambda0) << " nm\n";
    std::cout << "  CFL factor: " << P.S << "\n";
    std::cout << "  Time steps: " << P.nSteps << "\n";
    std::cout << "  Mesh mode: " << (P.mesh_mode == Config::MeshMode::AUTO_NONUNIFORM ? "AUTO_NONUNIFORM" : "UNIFORM") << "\n\n";

    // ===== 2. Register packages and build context =====
    Config::register_all_default_packages(P);
    Config::SimContext ctx = Config::make_context(P);
    Config::apply_structure_packages(ctx);
    Config::bake_and_bind(ctx, P, /*subpixel_samples=*/8);

    // ===== 3. Initialize sources and detectors =====
    Config::Runtime rt;
    Config::apply_source_packages(ctx, rt);
    Config::apply_detector_packages(ctx, rt);

    // ===== 4. Allocate field arrays =====
    const size_t NxT = ctx.NxT, NyT = ctx.NyT, NzT = ctx.NzT;
    const size_t Ntot = NxT * NyT * NzT;
    std::vector<real> Ex(Ntot, 0), Ey(Ntot, 0), Ez(Ntot, 0);
    std::vector<real> Hx(Ntot, 0), Hy(Ntot, 0), Hz(Ntot, 0);
    std::vector<real> Jx(Ntot, 0), Jy(Ntot, 0), Jz(Ntot, 0);

    // ===== 4.5. Initialize Parallel Domain Decomposition V2 (if enabled) =====
    // V2 uses TRUE domain decomposition: each subdomain has its own PML,
    // only halo regions are exchanged (O(N²) per step, not O(N³))
    ParallelV2::DomainDecomposition domain_decomp;
    ParallelV2::ParallelConfig parallel_config;
    parallel_config.enabled = UserConfig::PARALLEL_ENABLED;
    parallel_config.num_domains = UserConfig::PARALLEL_NUM_DOMAINS;
    parallel_config.axis = static_cast<ParallelV2::DecompAxis>(UserConfig::PARALLEL_DECOMP_AXIS);
    parallel_config.halo_width = UserConfig::PARALLEL_HALO_WIDTH;
    parallel_config.npml = ctx.bc->npml();

    if (parallel_config.enabled) {
        std::cout << "\n========================================\n";
        std::cout << "Initializing Parallel Domain Decomposition V2\n";
        std::cout << "========================================\n";

        domain_decomp.initialize(
            parallel_config,
            NxT, NyT, NzT,
            ctx.bc->npml(),
            ctx.grid_spacing,
            ctx.mats,
            P.dt, Config::eps0, Config::mu0, Config::c0
        );

        if (domain_decomp.is_parallel()) {
            std::cout << "Parallel mode: ACTIVE with " << domain_decomp.num_domains() << " domains\n";
            std::cout << "  - Each subdomain has independent PML\n";
            std::cout << "  - Only halo exchange per step (no scatter/gather)\n";
        } else {
            std::cout << "Parallel mode: INACTIVE (falling back to serial)\n";
        }
    }

    // Flag to track if we're using parallel mode
    const bool use_parallel = domain_decomp.is_parallel();

    // ===== 5. Initialize Two-Level System (if enabled) =====
    TwoLevelParams tls_params;  // Uses defaults from UserConfig
    TwoLevelState tls_state;

    if (UserConfig::TLS_ENABLED) {
        // Finalize TLS parameters
        tls_params.finalize();
        tls_params.finalize_with_dt(P.dt);

        // Allocate TLS state
        tls_state.allocate(NxT, NyT, NzT);

        // Calculate gain region bounds from ABSOLUTE PHYSICAL COORDINATES
        // Convert from absolute coordinates to coordinates relative to physical domain origin
        const real gain_x_min_rel = UserConfig::GAIN_X_MIN - P.domain_x_min;
        const real gain_x_max_rel = UserConfig::GAIN_X_MAX - P.domain_x_min;
        const real gain_y_min_rel = UserConfig::GAIN_Y_MIN - P.domain_y_min;
        const real gain_y_max_rel = UserConfig::GAIN_Y_MAX - P.domain_y_min;
        const real gain_z_min_rel = UserConfig::GAIN_Z_MIN - P.domain_z_min;
        const real gain_z_max_rel = UserConfig::GAIN_Z_MAX - P.domain_z_min;

        // Convert physical coordinates to grid indices
        const size_t i0 = ctx.grid_spacing.physical_to_index_x(gain_x_min_rel);
        const size_t i1 = ctx.grid_spacing.physical_to_index_x(gain_x_max_rel) + 1;
        const size_t j0 = ctx.grid_spacing.physical_to_index_y(gain_y_min_rel);
        const size_t j1 = ctx.grid_spacing.physical_to_index_y(gain_y_max_rel) + 1;
        const size_t k0 = ctx.grid_spacing.physical_to_index_z(gain_z_min_rel);
        const size_t k1 = ctx.grid_spacing.physical_to_index_z(gain_z_max_rel) + 1;

        std::cout << "Gain region physical coords: ["
                  << UserConfig::GAIN_X_MIN * 1e9 << ", " << UserConfig::GAIN_X_MAX * 1e9 << "] x ["
                  << UserConfig::GAIN_Y_MIN * 1e9 << ", " << UserConfig::GAIN_Y_MAX * 1e9 << "] x ["
                  << UserConfig::GAIN_Z_MIN * 1e9 << ", " << UserConfig::GAIN_Z_MAX * 1e9 << "] nm\n";
        std::cout << "Gain region grid indices: [" << i0 << "," << i1 << ") x ["
                  << j0 << "," << j1 << ") x [" << k0 << "," << k1 << ")\n";
        std::cout << "Gain region size: " << (i1 - i0) << " x " << (j1 - j0) << " x " << (k1 - k0) << " cells\n";

        tls_state.initialize_gain_region(
            i0, i1, j0, j1, k0, k1,
            tls_params.N0_total,
            ctx.grid_spacing,  // Pass grid spacing for cell volume calculation
            UserConfig::GAIN_INVERSION_FRACTION
        );

        std::cout << "Two-level gain region initialized\n";
    } else {
        std::cout << "Two-level system: DISABLED\n";
    }

    // ===== 6. Initialize TLS detectors (if enabled) =====
    fs::path out_root = fs::path("frames") / P.run_tag;
    // Calculate probe position indices from physical coordinates
    const real probe_x_rel = P.probe_x - P.domain_x_min;
    const real probe_y_rel = P.probe_y - P.domain_y_min;
    const real probe_z_rel = P.probe_z - P.domain_z_min;
    const size_t ic = ctx.grid_spacing.physical_to_index_x(probe_x_rel);
    const size_t jc = ctx.grid_spacing.physical_to_index_y(probe_y_rel);
    const size_t kc = ctx.grid_spacing.physical_to_index_z(probe_z_rel);

    std::unique_ptr<PopulationTimeSeriesDetector> pop_detector;
    std::unique_ptr<GainLossMonitor> gain_monitor;
    std::unique_ptr<PolarizationDiagnostics> pol_diagnostics;
    std::unique_ptr<TLSGlobalHistory> tls_global_history;

    if (UserConfig::TLS_ENABLED) {
        // VOLUME DENSITY FORMULATION: Pass GridSpacing for proper integration
        pop_detector = std::make_unique<PopulationTimeSeriesDetector>(
            out_root, "populations",
            NxT, NyT, NzT,
            P.saveEvery, P.nSteps, P.dt,
            ctx.grid_spacing  // For density × dV integration
        );

        gain_monitor = std::make_unique<GainLossMonitor>(
            out_root, "gain_monitor",
            NxT, NyT, NzT,
            P.saveEvery, P.dt,
            ctx.grid_spacing  // For density × dV integration
        );

        pol_diagnostics = std::make_unique<PolarizationDiagnostics>(
            out_root, "polarization",
            NxT, NyT, NzT,
            ic, jc, kc,
            P.saveEvery, P.dt
        );

        // TLS Global History - uses gain bounds from tls_state (unified source)
        tls_global_history = std::make_unique<TLSGlobalHistory>(
            out_root,
            tls_state,           // Pass state to get gain bounds
            ctx.grid_spacing,    // For density × dV integration
            P.saveEvery,
            P.dt,
            false                // auto_plot disabled by default (portable)
        );

        std::cout << "TLS detectors initialized\n";
    }

    // ===== 7. Initialize stability monitor =====
    StabilityMonitor stability_monitor;
    stability_monitor.print_config();

    // ===== 8. Main Time-Stepping Loop =====
    std::cout << "\n========================================\n";
    std::cout << "Starting FDTD Time-Stepping\n";
    if (use_parallel) {
        std::cout << "Mode: PARALLEL V2 (" << domain_decomp.num_domains() << " domains)\n";
        std::cout << "  - TRUE domain decomposition (each subdomain has own PML)\n";
        std::cout << "  - Only O(N^2) halo exchange per step (no O(N^3) scatter/gather)\n";
    } else {
#if FDTD_OMP_ENABLED
        std::cout << "Mode: SERIAL (OpenMP parallelized, " << omp_get_max_threads() << " threads)\n";
#else
        std::cout << "Mode: SERIAL\n";
#endif
    }
    std::cout << "========================================\n\n";

    // Timer for monitoring computation time per interval
    auto last_monitor_time = std::chrono::high_resolution_clock::now();

    // Initial scatter for parallel mode (only done ONCE at start)
    if (use_parallel) {
        domain_decomp.initial_scatter(Ex, Ey, Ez, Hx, Hy, Hz);

        // Initialize TLS in subdomains (if TLS enabled)
        if (UserConfig::TLS_ENABLED) {
            domain_decomp.initialize_tls(tls_state);
            std::cout << "[ParallelV2] TLS will be computed locally in each subdomain\n";
            std::cout << "  - NO gather/scatter needed for TLS during computation\n";
        }
    }

    for (size_t n = 0; n < P.nSteps; ++n) {

        // ==================== PARALLEL MODE V2 (True Domain Decomposition) ====================
        // Each subdomain has its own PML, only halo regions exchanged per step
        if (use_parallel) {
            // Step 1: Clear and inject source current into subdomains
            ParallelV2::parallel_clear_J(domain_decomp);
            for (auto& src : rt.sources) {
                // Inject into global arrays first, then scatter to appropriate subdomain
                src->inject_half_step(n, P.dt, Jx, Jy, Jz);
            }
            // Copy source currents to subdomains
            for (auto& dom : domain_decomp.subdomains) {
                for (size_t li = 0; li < dom.local_Nx; ++li) {
                    for (size_t lj = 0; lj < dom.local_Ny; ++lj) {
                        for (size_t lk = 0; lk < dom.local_Nz; ++lk) {
                            size_t gi, gj, gk;
                            dom.local_to_global(li, lj, lk, gi, gj, gk);
                            size_t lid = dom.local_idx(li, lj, lk);
                            size_t gid = idx3(gi, gj, gk, NyT, NzT);
                            dom.Jz[lid] = Jz[gid];
                        }
                    }
                }
            }

            // Step 2: Update H-field in all subdomains (includes PML correction)
            // Each subdomain handles its own PML internally
            ParallelV2::parallel_update_H(domain_decomp);

            // Step 3: Update E-field in all subdomains (includes PML and TLS)
            if (UserConfig::TLS_ENABLED) {
                // TLS is computed locally in each subdomain - NO gather/scatter needed!
                ParallelV2::parallel_update_E_with_tls(domain_decomp, P.dt, tls_params);
            } else {
                ParallelV2::parallel_update_E(domain_decomp);
            }

            // Step 4: Gather to global arrays only when needed for diagnostics/output
            // This is ONLY for monitoring/output, NOT for computation!
            bool need_global_data = (n % UserConfig::STABILITY_MONITOR_INTERVAL == 0) ||
                                    (n % P.saveEvery == 0);
            if (need_global_data) {
                domain_decomp.gather_fields_for_output(Ex, Ey, Ez, Hx, Hy, Hz);
                // Also gather TLS state for diagnostics
                if (UserConfig::TLS_ENABLED) {
                    domain_decomp.gather_tls_for_output(tls_state);
                }
            }

            // Callbacks (operate on global arrays)
            if (need_global_data) {
                rt.after_H(n, Ex, Ey, Ez, Hx, Hy, Hz);
                rt.after_E(n, Ex, Ey, Ez, Hx, Hy, Hz);
            }

        // ==================== SERIAL MODE ====================
        } else {
            // Step 1: Inject source current
            std::fill(Jx.begin(), Jx.end(), 0.0);
            std::fill(Jy.begin(), Jy.end(), 0.0);
            std::fill(Jz.begin(), Jz.end(), 0.0);
            for (auto& src : rt.sources) {
                src->inject_half_step(n, P.dt, Jx, Jy, Jz);
            }

            // Step 2: Update polarization (if TLS enabled)
            if (UserConfig::TLS_ENABLED) {
                update_polarization(NxT, NyT, NzT, tls_params, Ez, tls_state);
            }

            // Step 3: Update H-field
            fdtd_update_H(
                NxT, NyT, NzT,
                ctx.grid_spacing.dx.data(),
                ctx.grid_spacing.dy.data(),
                ctx.grid_spacing.dz.data(),
                ctx.grid_spacing.inv_dx.data(),
                ctx.grid_spacing.inv_dy.data(),
                ctx.grid_spacing.inv_dz.data(),
                ctx.mats.aHx.data(), ctx.mats.bHx.data(),
                ctx.mats.aHy.data(), ctx.mats.bHy.data(),
                ctx.mats.aHz.data(), ctx.mats.bHz.data(),
                Ex.data(), Ey.data(), Ez.data(),
                Hx.data(), Hy.data(), Hz.data()
            );

            ctx.bc->apply_after_H(Ex, Ey, Ez, Hx, Hy, Hz);
            rt.after_H(n, Ex, Ey, Ez, Hx, Hy, Hz);

            // Step 4: Update E-field
            if (UserConfig::TLS_ENABLED) {
                // With TLS: includes polarization source term -dP/dt
                fdtd_update_E_with_gain(
                    NxT, NyT, NzT,
                    ctx.grid_spacing.dx.data(),
                    ctx.grid_spacing.dy.data(),
                    ctx.grid_spacing.dz.data(),
                    ctx.grid_spacing.inv_dx.data(),
                    ctx.grid_spacing.inv_dy.data(),
                    ctx.grid_spacing.inv_dz.data(),
                    ctx.mats.aEx.data(), ctx.mats.bEx.data(),
                    ctx.mats.aEy.data(), ctx.mats.bEy.data(),
                    ctx.mats.aEz.data(), ctx.mats.bEz.data(),
                    Ex.data(), Ey.data(), Ez.data(),
                    Hx.data(), Hy.data(), Hz.data(),
                    Jx.data(), Jy.data(), Jz.data(),
                    tls_state
                );
            } else {
                // Without TLS: standard E-field update
                fdtd_update_E(
                    NxT, NyT, NzT,
                    ctx.grid_spacing.dx.data(),
                    ctx.grid_spacing.dy.data(),
                    ctx.grid_spacing.dz.data(),
                    ctx.grid_spacing.inv_dx.data(),
                    ctx.grid_spacing.inv_dy.data(),
                    ctx.grid_spacing.inv_dz.data(),
                    ctx.mats.aEx.data(), ctx.mats.bEx.data(),
                    ctx.mats.aEy.data(), ctx.mats.bEy.data(),
                    ctx.mats.aEz.data(), ctx.mats.bEz.data(),
                    Ex.data(), Ey.data(), Ez.data(),
                    Hx.data(), Hy.data(), Hz.data(),
                    Jx.data(), Jy.data(), Jz.data()
                );
            }

            ctx.bc->apply_after_E(Ex, Ey, Ez, Hx, Hy, Hz);
            rt.after_E(n, Ex, Ey, Ez, Hx, Hy, Hz);
        }

        // Step 5: Update populations (if TLS enabled)
        // For parallel mode: populations already updated in parallel_update_E_with_tls
        // For serial mode: update populations here
        if (UserConfig::TLS_ENABLED && !use_parallel) {
            update_populations(NxT, NyT, NzT, P.dt, tls_params, Ez, tls_state);
        }

        // Record TLS diagnostics (need global data which was gathered above)
        if (UserConfig::TLS_ENABLED) {
            // For parallel mode: only record when we have gathered data
            bool can_record = !use_parallel || (n % UserConfig::STABILITY_MONITOR_INTERVAL == 0) || (n % P.saveEvery == 0);
            if (can_record) {
                if (pop_detector) pop_detector->record(n, tls_state);
                if (gain_monitor) gain_monitor->record(n, tls_state);
                if (pol_diagnostics) pol_diagnostics->record(n, tls_state, ctx.grid_spacing);
                if (tls_global_history) tls_global_history->record(n, tls_state);
            }
        }

        // Stability check every STABILITY_MONITOR_INTERVAL steps
        if (n % UserConfig::STABILITY_MONITOR_INTERVAL == 0) {
            real conservation_error = stability_monitor.calculate_conservation_error(tls_state);

            bool stable = stability_monitor.check_stability(
                n, Ex, Ey, Ez, Hx, Hy, Hz, tls_state, conservation_error
            );

            if (!stable) {
                std::cerr << "[FATAL] Simulation terminated at step " << n << " due to instability.\n";
                stability_monitor.print_summary();
                return 1;
            }

            // Calculate elapsed time since last monitor
            auto current_time = std::chrono::high_resolution_clock::now();
            double elapsed_seconds = std::chrono::duration<double>(current_time - last_monitor_time).count();
            last_monitor_time = current_time;

            // Progress output
            real max_E = stability_monitor.get_max_E(Ex, Ey, Ez);

            if (UserConfig::TLS_ENABLED) {
                real max_P = NumericUtils::max_abs(tls_state.Pz);
                // VOLUME DENSITY FORMULATION: Use integrated functions
                real total_Nu = compute_integrated_Nu(tls_state, ctx.grid_spacing);
                real total_Ng = compute_integrated_Ng(tls_state, ctx.grid_spacing);
                real inversion = total_Nu - total_Ng;
                std::cout << "step " << n << " | max|E| = " << std::scientific << std::setprecision(3) << max_E
                          << " | max|P| = " << max_P
                          << " | Nu = " << std::fixed << std::setprecision(2) << total_Nu
                          << " | Ng = " << total_Ng
                          << " | inv = " << std::setprecision(2) << inversion
                          << " | err = " << std::scientific << conservation_error
                          << " | time = " << std::fixed << std::setprecision(3) << elapsed_seconds << " s"
                          << "\n";
            } else {
                std::cout << "step " << n << " | max|E| = " << std::scientific << std::setprecision(3) << max_E
                          << " | time = " << std::fixed << std::setprecision(3) << elapsed_seconds << " s"
                          << "\n";
            }
        }

        // Save progress indicator (using probe position)
        if (n % P.saveEvery == 0 && n % UserConfig::STABILITY_MONITOR_INTERVAL != 0) {
            // Use probe position for progress output
            std::cout << "step " << n << " | Ez(probe) = "
                      << std::scientific << Ez[idx3(ic, jc, kc, NyT, NzT)] << std::fixed << "\n";
        }
    }

    // ===== 9. Finalization =====
    stability_monitor.print_summary();

    // Finalize TLS global history (generate plots)
    if (tls_global_history) {
        tls_global_history->finalize();
    }

    std::cout << "\n========================================\n";
    std::cout << "Simulation Complete\n";
    std::cout << "========================================\n";
    std::cout << "  dt = " << UnitConv::s_to_fs(P.dt) << " fs\n";
    std::cout << "  t_total = " << UnitConv::s_to_fs(P.dt * P.nSteps) << " fs\n";
    std::cout << "  Output: frames/" << P.run_tag << "/\n";

    return 0;
}