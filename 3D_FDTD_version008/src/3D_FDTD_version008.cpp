// 3D_FDTD_version008.cpp - Main program for 3D FDTD electromagnetic simulation
// with optional two-level atomic gain medium coupling.
// Configuration: Modify user_config.hpp to change simulation parameters.

#include "user_config.hpp"
#include "params.hpp"
#include "fdtd_stepper.hpp"
#include "boundary.hpp"
#include "structure_material.hpp"
#include "global_function.hpp"
#include "omp_config.hpp"
#include "two_level_system.hpp"
#include "two_level_detectors.hpp"
#include "stability_monitor.hpp"
#include "parallel_domain.hpp"

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

    Config::register_all_default_packages(P);
    Config::SimContext ctx = Config::make_context(P);
    Config::apply_structure_packages(ctx);
    Config::bake_and_bind(ctx, P, /*subpixel_samples=*/8);

    Config::Runtime rt;
    Config::apply_source_packages(ctx, rt);
    Config::apply_detector_packages(ctx, rt);

    const size_t NxT = ctx.NxT, NyT = ctx.NyT, NzT = ctx.NzT;
    const size_t Ntot = NxT * NyT * NzT;
    std::vector<real> Ex(Ntot, 0), Ey(Ntot, 0), Ez(Ntot, 0);
    std::vector<real> Hx(Ntot, 0), Hy(Ntot, 0), Hz(Ntot, 0);
    std::vector<real> Jx(Ntot, 0), Jy(Ntot, 0), Jz(Ntot, 0);

    Parallel::DomainDecomposition domain_decomp;
    Parallel::ParallelConfig parallel_config;
    parallel_config.enabled = UserConfig::PARALLEL_ENABLED;
    parallel_config.num_domains = UserConfig::PARALLEL_NUM_DOMAINS;
    parallel_config.axis = static_cast<Parallel::DecompAxis>(UserConfig::PARALLEL_DECOMP_AXIS);
    parallel_config.halo_width = UserConfig::PARALLEL_HALO_WIDTH;
    parallel_config.npml = ctx.bc->npml();

    if (parallel_config.enabled) {
        std::cout << "\n========================================\n";
        std::cout << "Initializing Parallel Domain Decomposition\n";
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

    const bool use_parallel = domain_decomp.is_parallel();

    TwoLevelParams tls_params;
    TwoLevelState tls_state;
    TLSRegionManager tls_manager;
    bool any_tls_enabled = UserConfig::TLS_ENABLED && ctx.scene.has_any_tls();

    if (any_tls_enabled) {
        std::cout << "\n========================================\n";
        std::cout << "Initializing Structure-Bound TLS\n";
        std::cout << "========================================\n";

        tls_state.allocate(NxT, NyT, NzT);
        auto tls_structures = ctx.scene.get_tls_structures();
        std::cout << "[TLS] Found " << tls_structures.size() << " structure(s) with TLS\n";

        for (size_t idx = 0; idx < tls_structures.size(); ++idx) {
            const auto* item = tls_structures[idx];
            const auto& tls_cfg = item->tls_config;
            AABB bb = item->shape->bounding_box();

            real pml_x = ctx.grid_spacing.x_bounds[ctx.npml];
            real pml_y = ctx.grid_spacing.y_bounds[ctx.npml];
            real pml_z = ctx.grid_spacing.z_bounds[ctx.npml];

            TwoLevelParams region_params;
            region_params.lambda0 = tls_cfg.lambda0;
            region_params.gamma = tls_cfg.gamma;
            region_params.tau = tls_cfg.tau;
            region_params.N0_total = tls_cfg.N0;
            region_params.enable_population_clamp = UserConfig::TLS_ENABLE_CLAMP;

            std::string region_name = "structure_" + std::to_string(idx);
            tls_manager.add_region_from_structure(
                ctx.grid_spacing,
                bb.x0 - pml_x, bb.x1 - pml_x,
                bb.y0 - pml_y, bb.y1 - pml_y,
                bb.z0 - pml_z, bb.z1 - pml_z,
                region_params,
                tls_cfg.inversion_fraction,
                item->shape.get(),
                region_name
            );
        }

        tls_manager.initialize(tls_state, ctx.grid_spacing);
        tls_manager.finalize_with_dt(P.dt);

        if (tls_manager.num_regions() > 0) {
            tls_params = tls_manager.regions[0].params;
            tls_state.gain_i0 = tls_manager.regions[0].i0;
            tls_state.gain_i1 = tls_manager.regions[0].i1;
            tls_state.gain_j0 = tls_manager.regions[0].j0;
            tls_state.gain_j1 = tls_manager.regions[0].j1;
            tls_state.gain_k0 = tls_manager.regions[0].k0;
            tls_state.gain_k1 = tls_manager.regions[0].k1;

            for (size_t r = 1; r < tls_manager.num_regions(); ++r) {
                const auto& region = tls_manager.regions[r];
                tls_state.gain_i0 = std::min(tls_state.gain_i0, region.i0);
                tls_state.gain_i1 = std::max(tls_state.gain_i1, region.i1);
                tls_state.gain_j0 = std::min(tls_state.gain_j0, region.j0);
                tls_state.gain_j1 = std::max(tls_state.gain_j1, region.j1);
                tls_state.gain_k0 = std::min(tls_state.gain_k0, region.k0);
                tls_state.gain_k1 = std::max(tls_state.gain_k1, region.k1);
            }
        }
        std::cout << "========================================\n\n";
    } else {
        if (!UserConfig::TLS_ENABLED) {
            std::cout << "Two-level system: DISABLED (TLS_ENABLED = false)\n";
        } else {
            std::cout << "Two-level system: DISABLED (no structures with has_tls = true)\n";
        }
    }

    fs::path out_root = fs::path("frames") / P.run_tag;
    const real probe_x_rel = P.probe_x - P.domain_x_min;
    const real probe_y_rel = P.probe_y - P.domain_y_min;
    const real probe_z_rel = P.probe_z - P.domain_z_min;
    const size_t ic = ctx.grid_spacing.physical_to_index_x(probe_x_rel);
    const size_t jc = ctx.grid_spacing.physical_to_index_y(probe_y_rel);
    const size_t kc = ctx.grid_spacing.physical_to_index_z(probe_z_rel);

    std::unique_ptr<PopulationTimeSeriesDetector> pop_detector;
    std::unique_ptr<PolarizationDiagnostics> pol_diagnostics;
    std::unique_ptr<TLSGlobalHistory> tls_global_history;

    if (any_tls_enabled) {
        pop_detector = std::make_unique<PopulationTimeSeriesDetector>(
            out_root, "populations", NxT, NyT, NzT,
            P.saveEvery, P.nSteps, P.dt, ctx.grid_spacing);

        pol_diagnostics = std::make_unique<PolarizationDiagnostics>(
            out_root, "polarization", NxT, NyT, NzT,
            ic, jc, kc, P.saveEvery, P.dt);

        tls_global_history = std::make_unique<TLSGlobalHistory>(
            out_root, tls_state, ctx.grid_spacing, P.saveEvery, P.dt, false);

        std::cout << "TLS detectors initialized (" << tls_manager.num_regions() << " region(s))\n";
    }

    StabilityMonitor stability_monitor;
    stability_monitor.print_config();

    std::cout << "\n========================================\n";
    std::cout << "Starting FDTD Time-Stepping\n";
    if (use_parallel) {
        std::cout << "Mode: PARALLEL (" << domain_decomp.num_domains() << " domains)\n";
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

    auto last_monitor_time = std::chrono::high_resolution_clock::now();

    if (use_parallel) {
        domain_decomp.initial_scatter(Ex, Ey, Ez, Hx, Hy, Hz);
        if (any_tls_enabled) {
            domain_decomp.initialize_tls(tls_state);
            std::cout << "[Parallel] TLS will be computed locally in each subdomain\n";
            std::cout << "  - " << tls_manager.num_regions() << " TLS region(s) bound to structures\n";
            std::cout << "  - NO gather/scatter needed for TLS during computation\n";
        }
    }

    for (size_t n = 0; n < P.nSteps; ++n) {
        if (use_parallel) {
            Parallel::parallel_clear_J(domain_decomp);
            for (auto& src : rt.sources) {
                src->inject_half_step(n, P.dt, Jx, Jy, Jz);
            }
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

            Parallel::parallel_update_H(domain_decomp);

            if (any_tls_enabled) {
                Parallel::parallel_update_E_with_tls(domain_decomp, P.dt, tls_params);
            } else {
                Parallel::parallel_update_E(domain_decomp);
            }

            bool need_global_data = (n % UserConfig::STABILITY_MONITOR_INTERVAL == 0) ||
                                    (n % P.saveEvery == 0);
            if (need_global_data) {
                domain_decomp.gather_fields_for_output(Ex, Ey, Ez, Hx, Hy, Hz);
                if (any_tls_enabled) {
                    domain_decomp.gather_tls_for_output(tls_state);
                }
                rt.after_H(n, Ex, Ey, Ez, Hx, Hy, Hz);
                rt.after_E(n, Ex, Ey, Ez, Hx, Hy, Hz);
            }
        } else {
            std::fill(Jx.begin(), Jx.end(), 0.0);
            std::fill(Jy.begin(), Jy.end(), 0.0);
            std::fill(Jz.begin(), Jz.end(), 0.0);
            for (auto& src : rt.sources) {
                src->inject_half_step(n, P.dt, Jx, Jy, Jz);
            }

            if (any_tls_enabled) {
                update_polarization_multi_region(NxT, NyT, NzT, tls_manager, Ez, tls_state);
            }

            fdtd_update_H(
                NxT, NyT, NzT,
                ctx.grid_spacing.dx.data(), ctx.grid_spacing.dy.data(), ctx.grid_spacing.dz.data(),
                ctx.grid_spacing.inv_dx.data(), ctx.grid_spacing.inv_dy.data(), ctx.grid_spacing.inv_dz.data(),
                ctx.mats.aHx.data(), ctx.mats.bHx.data(),
                ctx.mats.aHy.data(), ctx.mats.bHy.data(),
                ctx.mats.aHz.data(), ctx.mats.bHz.data(),
                Ex.data(), Ey.data(), Ez.data(),
                Hx.data(), Hy.data(), Hz.data()
            );

            ctx.bc->apply_after_H(Ex, Ey, Ez, Hx, Hy, Hz);
            rt.after_H(n, Ex, Ey, Ez, Hx, Hy, Hz);

            if (any_tls_enabled) {
                fdtd_update_E_with_gain(
                    NxT, NyT, NzT,
                    ctx.grid_spacing.dx.data(), ctx.grid_spacing.dy.data(), ctx.grid_spacing.dz.data(),
                    ctx.grid_spacing.inv_dx.data(), ctx.grid_spacing.inv_dy.data(), ctx.grid_spacing.inv_dz.data(),
                    ctx.mats.aEx.data(), ctx.mats.bEx.data(),
                    ctx.mats.aEy.data(), ctx.mats.bEy.data(),
                    ctx.mats.aEz.data(), ctx.mats.bEz.data(),
                    Ex.data(), Ey.data(), Ez.data(),
                    Hx.data(), Hy.data(), Hz.data(),
                    Jx.data(), Jy.data(), Jz.data(),
                    tls_state
                );
            } else {
                fdtd_update_E(
                    NxT, NyT, NzT,
                    ctx.grid_spacing.dx.data(), ctx.grid_spacing.dy.data(), ctx.grid_spacing.dz.data(),
                    ctx.grid_spacing.inv_dx.data(), ctx.grid_spacing.inv_dy.data(), ctx.grid_spacing.inv_dz.data(),
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

        if (any_tls_enabled && !use_parallel) {
            update_populations_multi_region(NxT, NyT, NzT, P.dt, tls_manager, Ez, tls_state);
        }

        if (any_tls_enabled) {
            bool can_record = !use_parallel || (n % UserConfig::STABILITY_MONITOR_INTERVAL == 0) || (n % P.saveEvery == 0);
            if (can_record) {
                if (pop_detector) pop_detector->record(n, tls_state);
                if (pol_diagnostics) pol_diagnostics->record(n, tls_state, Ez, ctx.grid_spacing);
                if (tls_global_history) tls_global_history->record(n, tls_state);
            }
        }

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

            auto current_time = std::chrono::high_resolution_clock::now();
            double elapsed_seconds = std::chrono::duration<double>(current_time - last_monitor_time).count();
            last_monitor_time = current_time;

            real max_E = stability_monitor.get_max_E(Ex, Ey, Ez);

            if (any_tls_enabled) {
                real max_P = NumericUtils::max_abs(tls_state.Pz);
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

        if (n % P.saveEvery == 0 && n % UserConfig::STABILITY_MONITOR_INTERVAL != 0) {
            std::cout << "step " << n << " | Ez(probe) = "
                      << std::scientific << Ez[idx3(ic, jc, kc, NyT, NzT)] << std::fixed << "\n";
        }
    }

    stability_monitor.print_summary();
    if (pop_detector) {
        pop_detector->finalize(tls_state);
    }
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
