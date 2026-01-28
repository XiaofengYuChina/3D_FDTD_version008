// two_level_system.hpp - Two-level atomic system for active gain media
// Implements laser physics: population inversion, stimulated emission, etc.

#pragma once

#include <vector>
#include <cmath>
#include <numbers>
#include <iostream>
#include "global_function.hpp"
#include "user_config.hpp"
#include "omp_config.hpp"
#include "structure_material.hpp"

struct StructureTLSConfig;

struct TwoLevelParams {
    real lambda0 = 1500e-9;     // Transition wavelength (m)
    real omega_a = 0.0;         // Angular frequency (rad/s)
    real gamma = 7e12;          // Polarization damping rate (s^-1)
    real tau = 1e-12;           // Upper level lifetime (s)
    real mu_z = 0.0;            // Transition dipole moment
    real N0_total = 1e25;       // Total dipole density (m^-3)
    real hbar = PhysConst::HBAR;
    bool enable_population_clamp = UserConfig::TLS_ENABLE_CLAMP;

    // Three-point recursion coefficients for damped oscillator:
    // d^2P/dt^2 + gamma*dP/dt + omega^2*P = F
    // Discretized: P^{n+1} = kapa*F^n + pa2*P^n + pa3*P^{n-1}
    real kapa_coeff = 0.0;
    real pa2 = 0.0;
    real pa3 = 0.0;
    real inv_dt = 0.0;
    bool coefficients_initialized = false;

    real Az_value() const { return 0.0; }  // Low-intensity approximation

    void finalize(real c0 = 299792458.0, real eps0 = 8.854187817e-12, real hbar_in = 1.054571817e-34) {
        hbar = hbar_in;
        omega_a = 2.0 * std::numbers::pi * c0 / lambda0;

        // |mu_z| = sqrt(3*pi*eps0*hbar*c^3 / (omega_a^3 * tau))
        real numerator = 3.0 * std::numbers::pi * eps0 * hbar * c0 * c0 * c0;
        real denominator = omega_a * omega_a * omega_a * tau;
        mu_z = std::sqrt(numerator / denominator);

        std::cout << "\n=== Two-Level System Parameter Check ===\n";
        std::cout << "  tau (upper level lifetime) = " << tau << " s\n";
        std::cout << "  mu_z (dipole moment)       = " << mu_z << " C·m\n";
        std::cout << "  omega_a (angular freq)     = " << omega_a << " rad/s\n";
        std::cout << "  gamma (damping rate)       = " << gamma << " s⁻¹\n";

        real coupling_indicator = mu_z * std::sqrt(N0_total / (eps0 * hbar * omega_a));
        std::cout << "  Coupling strength indicator = " << coupling_indicator << "\n";

        // Small tau leads to large mu_z, causing numerical instability
        if (tau < 1e-12) {
            std::cerr << "\n*** WARNING (Issue D): tau = " << tau << " s is very small! ***\n";
            std::cerr << "    This leads to extremely large mu_z = " << mu_z << " C·m\n";
            std::cerr << "    Strong coupling may cause numerical instability/divergence.\n";
            std::cerr << "    Recommended range: tau = 1e-12 to 1e-9 s\n";
            std::cerr << "    Consider increasing tau or reducing time step dt.\n\n";
        }
        if (coupling_indicator > 1e6) {
            std::cerr << "*** WARNING: Coupling indicator = " << coupling_indicator << " is very large! ***\n";
            std::cerr << "    This may cause rapid energy exchange and numerical issues.\n\n";
        }
        std::cout << "========================================\n\n";
    }

    void finalize_with_dt(real dt) {
        if (omega_a <= 0.0) {
            std::cerr << "ERROR: finalize_with_dt called before finalize()!\n";
            return;
        }

        const real dt2 = dt * dt;
        const real omega_sq = omega_a * omega_a;
        const real beta = 2.0 + gamma * dt;

        const real two_omega_over_hbar = 2.0 * omega_a / hbar;
        const real mu_z_sq = mu_z * mu_z;
        kapa_coeff = (2.0 * dt2 / beta) * two_omega_over_hbar * mu_z_sq;
        pa2 = 2.0 * (2.0 - omega_sq * dt2) / beta;
        pa3 = (gamma * dt - 2.0) / beta;
        inv_dt = 1.0 / dt;
        coefficients_initialized = true;

        std::cout << "=== Three-point recursion coefficients ===\n";
        std::cout << "  dt = " << dt * 1e15 << " fs\n";
        std::cout << "  beta (2+γ·dt) = " << beta << "\n";
        std::cout << "  kapa_coeff = " << kapa_coeff << "\n";
        std::cout << "  pa2 = " << pa2 << "\n";
        std::cout << "  pa3 = " << pa3 << "\n";

        real omega_dt = omega_a * dt;
        std::cout << "  ω·dt = " << omega_dt << " (should be < 2 for stability)\n";

        if (omega_dt >= 2.0) {
            std::cerr << "*** WARNING: ω·dt = " << omega_dt << " >= 2 ***\n";
            std::cerr << "    Polarization dynamics may be unstable!\n";
            std::cerr << "    Reduce time step dt or transition frequency ω_a.\n";
        } else if (omega_dt > 1.5) {
            std::cout << "  (Approaching stability limit, consider smaller dt)\n";
        } else {
            std::cout << "  (Stable regime)\n";
        }
        std::cout << "==========================================\n\n";
    }
};

struct TwoLevelState {
    size_t NxT, NyT, NzT;
    size_t gain_i0 = 0, gain_i1 = 0;
    size_t gain_j0 = 0, gain_j1 = 0;
    size_t gain_k0 = 0, gain_k1 = 0;

    // Population densities [m^-3] (volume density formulation)
    std::vector<real> Ng;       // Ground state density
    std::vector<real> Nu;       // Upper state density
    std::vector<real> Ng0;      // Total density (conserved)

    // Polarization - three-point recursion storage
    std::vector<real> Pz;       // Current polarization
    std::vector<real> Pz_prev;  // Previous time step
    std::vector<real> dPz_dt;   // Time derivative for E-field update

    std::vector<real> Ez_old;   // Previous Ez for stimulated term calculation
    std::vector<real> Ndip;     // Dipole density distribution [m^-3]

    void allocate(size_t Nx, size_t Ny, size_t Nz) {
        NxT = Nx; NyT = Ny; NzT = Nz;
        const size_t N = NxT * NyT * NzT;
        Ng.assign(N, 0.0);
        Nu.assign(N, 0.0);
        Ng0.assign(N, 0.0);
        Pz.assign(N, 0.0);
        Pz_prev.assign(N, 0.0);
        dPz_dt.assign(N, 0.0);
        Ez_old.assign(N, 0.0);
        Ndip.assign(N, 0.0);
    }

    void initialize_gain_region(
        size_t i0, size_t i1, size_t j0, size_t j1, size_t k0, size_t k1,
        real N0_density, const GridSpacing& grid, real inversion_fraction = 0.1
    ) {
        i1 = std::min(i1, NxT);
        j1 = std::min(j1, NyT);
        k1 = std::min(k1, NzT);

        if (i0 >= i1 || j0 >= j1 || k0 >= k1) {
            std::cerr << "Warning: Invalid gain region bounds, skipping initialization\n";
            return;
        }

        gain_i0 = i0; gain_i1 = i1;
        gain_j0 = j0; gain_j1 = j1;
        gain_k0 = k0; gain_k1 = k1;

        real total_volume = 0.0;
        real min_cell_vol = 1e30, max_cell_vol = 0.0;

        for (size_t i = i0; i < i1; ++i) {
            for (size_t j = j0; j < j1; ++j) {
                for (size_t k = k0; k < k1; ++k) {
                    size_t id = idx3(i, j, k, NyT, NzT);
                    Ndip[id] = N0_density;
                    Nu[id] = inversion_fraction * N0_density;
                    Ng[id] = (1.0 - inversion_fraction) * N0_density;
                    Ng0[id] = N0_density;

                    real cell_volume = grid.dx[i] * grid.dy[j] * grid.dz[k];
                    total_volume += cell_volume;
                    min_cell_vol = std::min(min_cell_vol, cell_volume);
                    max_cell_vol = std::max(max_cell_vol, cell_volume);
                }
            }
        }

        real total_atoms = N0_density * total_volume;
        real inversion_density = (2.0 * inversion_fraction - 1.0) * N0_density;
        real integrated_inversion = inversion_density * total_volume;
        size_t n_cells = (i1 - i0) * (j1 - j0) * (k1 - k0);

        std::cout << "\n[TLS Init] Gain region initialized (DENSITY FORMULATION):\n";
        std::cout << "  N0_density = " << N0_density << " atoms/m³\n";
        std::cout << "  Gain region volume V_total = " << total_volume * 1e18 << " μm³\n";
        std::cout << "  Equivalent total atoms N = " << total_atoms << "\n";
        std::cout << "  Cell volume range: " << min_cell_vol * 1e27 << " - "
                  << max_cell_vol * 1e27 << " nm³\n";
        std::cout << "  Number of cells: " << n_cells << "\n";
        std::cout << "  Initial inversion fraction: " << inversion_fraction << "\n";
        std::cout << "  Inversion density (Nu-Ng): " << inversion_density << " m^-3\n";
        std::cout << "  Integrated inversion: " << integrated_inversion << " atoms\n\n";
    }
};

// Update polarization using three-point recursion for damped harmonic oscillator
// d^2Pz/dt^2 + gamma*dPz/dt + omega_a^2*Pz = (2*omega_a/hbar)*|mu_z|^2*Ndip*[(Ng-Nu)/Ng0]*Ez
inline void update_polarization(
    size_t NxT, size_t NyT, size_t NzT,
    const TwoLevelParams& params,
    const std::vector<real>& Ez,
    TwoLevelState& state
) {
    if (!params.coefficients_initialized) {
        std::cerr << "ERROR: TwoLevelParams::finalize_with_dt() not called!\n";
        return;
    }

    const real kapa_coeff = params.kapa_coeff;
    const real pa2 = params.pa2;
    const real pa3 = params.pa3;
    const real inv_dt = params.inv_dt;

#if FDTD_OMP_ENABLED
#pragma omp parallel for
#endif
    for (int i = 0; i < (int)NxT; ++i) {
        for (int j = 0; j < (int)NyT; ++j) {
            for (int k = 0; k < (int)NzT; ++k) {
                size_t id = idx3((size_t)i, (size_t)j, (size_t)k, NyT, NzT);
                if (state.Ndip[id] <= 0) continue;

                // Population inversion: positive = absorption, negative = gain
                real Ng_frac = (state.Ng0[id] > 0) ?
                    (state.Ng[id] - state.Nu[id]) / state.Ng0[id] : 0.0;

                real driving_factor = state.Ndip[id] * Ng_frac * Ez[id];
                real Pz_new = kapa_coeff * driving_factor + pa2 * state.Pz[id] + pa3 * state.Pz_prev[id];
                real Pz_old = state.Pz[id];
                real dPz_dt_half = (Pz_new - Pz_old) * inv_dt;

                state.Pz_prev[id] = Pz_old;
                state.Pz[id] = Pz_new;
                state.dPz_dt[id] = dPz_dt_half;
            }
        }
    }
}

// Update population densities using two-level rate equations
// stim_rate = E_avg*(dP/dt)/(hbar*omega_a)
// dNu/dt = -Nu/tau + stim_rate, dNg/dt = +Nu/tau - stim_rate
inline void update_populations(
    size_t NxT, size_t NyT, size_t NzT, real dt,
    const TwoLevelParams& params,
    const std::vector<real>& Ez,
    TwoLevelState& state
) {
    const real inv_tau = 1.0 / params.tau;
    const real inv_hbar_omega = 1.0 / (params.hbar * params.omega_a);
    const real inv_dt = params.inv_dt;

#if FDTD_OMP_ENABLED
#pragma omp parallel for
#endif
    for (int i = 0; i < (int)NxT; ++i) {
        for (int j = 0; j < (int)NyT; ++j) {
            for (int k = 0; k < (int)NzT; ++k) {
                size_t id = idx3((size_t)i, (size_t)j, (size_t)k, NyT, NzT);
                if (state.Ndip[id] <= 0) continue;

                real Nu_curr = state.Nu[id];
                real Ntotal = state.Ng0[id];
                if (Ntotal <= 0) continue;

                real E_avg = 0.5 * (Ez[id] + state.Ez_old[id]);
                real delta_P = state.Pz[id] - state.Pz_prev[id];
                real dP_dt = delta_P * inv_dt;
                real stim_rate = E_avg * dP_dt * inv_hbar_omega;

                real dNu_dt = -inv_tau * Nu_curr + stim_rate;
                real dNg_dt = +inv_tau * Nu_curr - stim_rate;

                state.Nu[id] += dt * dNu_dt;
                state.Ng[id] += dt * dNg_dt;

                state.Nu[id] = std::max(0.0, state.Nu[id]);
                state.Ng[id] = std::max(0.0, state.Ng[id]);

                if (params.enable_population_clamp) {
                    real total = state.Nu[id] + state.Ng[id];
                    if (total > 1e-30) {
                        real scale = Ntotal / total;
                        state.Nu[id] *= scale;
                        state.Ng[id] *= scale;
                    }
                    state.Nu[id] = std::min(state.Nu[id], Ntotal);
                    state.Ng[id] = Ntotal - state.Nu[id];
                }
            }
        }
    }
}

// E-field update with polarization source: dE/dt = (1/eps0*n^2)[curl(H) - dP/dt - J]
template<typename Real>
inline void fdtd_update_E_with_gain(
    size_t NxT, size_t NyT, size_t NzT,
    const Real* __restrict dx_array, const Real* __restrict dy_array, const Real* __restrict dz_array,
    const Real* __restrict inv_dx_array, const Real* __restrict inv_dy_array, const Real* __restrict inv_dz_array,
    const Real* __restrict aEx, const Real* __restrict bEx,
    const Real* __restrict aEy, const Real* __restrict bEy,
    const Real* __restrict aEz, const Real* __restrict bEz,
    Real* __restrict Ex, Real* __restrict Ey, Real* __restrict Ez,
    Real* __restrict Hx, Real* __restrict Hy, Real* __restrict Hz,
    const Real* __restrict Jx, const Real* __restrict Jy, const Real* __restrict Jz,
    TwoLevelState& state
) {
    using namespace fdtd_math;
    const size_t sI = NyT * NzT;
    const size_t sJ = NzT;
    const size_t sK = 1;

#if FDTD_OMP_ENABLED
#pragma omp parallel for
#endif
    for (int i = 1; i < (int)NxT; ++i) {
        for (int j = 1; j < (int)NyT; ++j) {
            for (int k = 1; k < (int)NzT; ++k) {
                size_t id = idx3((size_t)i, (size_t)j, (size_t)k, NyT, NzT);

                const Real inv_dx = inv_dx_array[i];
                const Real inv_dy = inv_dy_array[j];
                const Real inv_dz = inv_dz_array[k];

                Real curlHx = diff_ym(Hz, id, sJ, inv_dy) - diff_zm(Hy, id, sK, inv_dz);
                Real curlHy = diff_zm(Hx, id, sK, inv_dz) - diff_xm(Hz, id, sI, inv_dx);
                Real curlHz = diff_xm(Hy, id, sI, inv_dx) - diff_ym(Hx, id, sJ, inv_dy);

                Ex[id] = aEx[id] * Ex[id] + bEx[id] * (curlHx - Jx[id]);
                Ey[id] = aEy[id] * Ey[id] + bEy[id] * (curlHy - Jy[id]);

                state.Ez_old[id] = Ez[id];
                Real Pz_source = state.dPz_dt[id];
                Ez[id] = aEz[id] * Ez[id] + bEz[id] * (curlHz - Jz[id] - Pz_source);
            }
        }
    }
}

inline real compute_integrated_Nu(const TwoLevelState& state, const GridSpacing& grid) {
    real total_Nu = 0.0;
    for (size_t i = 0; i < state.NxT; ++i) {
        for (size_t j = 0; j < state.NyT; ++j) {
            for (size_t k = 0; k < state.NzT; ++k) {
                size_t id = idx3(i, j, k, state.NyT, state.NzT);
                if (state.Ndip[id] <= 0) continue;
                total_Nu += state.Nu[id] * grid.dx[i] * grid.dy[j] * grid.dz[k];
            }
        }
    }
    return total_Nu;
}

inline real compute_integrated_Ng(const TwoLevelState& state, const GridSpacing& grid) {
    real total_Ng = 0.0;
    for (size_t i = 0; i < state.NxT; ++i) {
        for (size_t j = 0; j < state.NyT; ++j) {
            for (size_t k = 0; k < state.NzT; ++k) {
                size_t id = idx3(i, j, k, state.NyT, state.NzT);
                if (state.Ndip[id] <= 0) continue;
                total_Ng += state.Ng[id] * grid.dx[i] * grid.dy[j] * grid.dz[k];
            }
        }
    }
    return total_Ng;
}

inline real compute_integrated_inversion(const TwoLevelState& state, const GridSpacing& grid) {
    real total_inversion = 0.0;
    for (size_t i = 0; i < state.NxT; ++i) {
        for (size_t j = 0; j < state.NyT; ++j) {
            for (size_t k = 0; k < state.NzT; ++k) {
                size_t id = idx3(i, j, k, state.NyT, state.NzT);
                if (state.Ndip[id] <= 0) continue;
                total_inversion += (state.Nu[id] - state.Ng[id]) * grid.dx[i] * grid.dy[j] * grid.dz[k];
            }
        }
    }
    return total_inversion;
}

inline real compute_avg_inversion_density(const TwoLevelState& state, const GridSpacing& grid) {
    real total_inversion = 0.0;
    real total_volume = 0.0;
    for (size_t i = 0; i < state.NxT; ++i) {
        for (size_t j = 0; j < state.NyT; ++j) {
            for (size_t k = 0; k < state.NzT; ++k) {
                size_t id = idx3(i, j, k, state.NyT, state.NzT);
                if (state.Ndip[id] <= 0) continue;
                real dV = grid.dx[i] * grid.dy[j] * grid.dz[k];
                total_inversion += (state.Nu[id] - state.Ng[id]) * dV;
                total_volume += dV;
            }
        }
    }
    return (total_volume > 0) ? total_inversion / total_volume : 0.0;
}

inline real compute_polarization_energy(
    const TwoLevelState& state, const std::vector<real>& Ez, const GridSpacing& grid
) {
    real total_energy = 0.0;
    for (size_t i = 0; i < state.NxT; ++i) {
        for (size_t j = 0; j < state.NyT; ++j) {
            for (size_t k = 0; k < state.NzT; ++k) {
                size_t id = idx3(i, j, k, state.NyT, state.NzT);
                if (state.Ndip[id] <= 0) continue;
                real dV = grid.dx[i] * grid.dy[j] * grid.dz[k];
                total_energy += -0.5 * state.Pz[id] * Ez[id] * dV;
            }
        }
    }
    return total_energy;
}

// TLS region bound to a structure
struct TLSRegion {
    size_t region_id = 0;
    std::string name = "region_0";
    size_t i0 = 0, i1 = 0;
    size_t j0 = 0, j1 = 0;
    size_t k0 = 0, k1 = 0;
    TwoLevelParams params;
    real inversion_fraction = 1.0;
    const Shape* shape = nullptr;  // Shape geometry for precise containment check

    bool contains(size_t i, size_t j, size_t k, const GridSpacing& grid) const {
        if (!(i >= i0 && i < i1 && j >= j0 && j < j1 && k >= k0 && k < k1))
            return false;
        if (shape) {
            real x = grid.cell_center_physical_x(i);
            real y = grid.cell_center_physical_y(j);
            real z = grid.cell_center_physical_z(k);
            return shape->contains(x, y, z);
        }
        return true;  // Fallback: bounding box only
    }

    void get_physical_bounds(const GridSpacing& grid,
                             real& x_min, real& x_max,
                             real& y_min, real& y_max,
                             real& z_min, real& z_max) const {
        x_min = grid.x_bounds[i0]; x_max = grid.x_bounds[i1];
        y_min = grid.y_bounds[j0]; y_max = grid.y_bounds[j1];
        z_min = grid.z_bounds[k0]; z_max = grid.z_bounds[k1];
    }

    void initialize(TwoLevelState& state, const GridSpacing& grid) {
        if (i0 >= i1 || j0 >= j1 || k0 >= k1) {
            std::cerr << "[TLSRegion " << name << "] Warning: Invalid bounds, skipping\n";
            return;
        }

        real total_volume = 0.0;
        real N0_density = params.N0_total;
        size_t n_cells = 0;

        for (size_t i = i0; i < i1; ++i) {
            for (size_t j = j0; j < j1; ++j) {
                for (size_t k = k0; k < k1; ++k) {
                    if (!contains(i, j, k, grid)) continue;
                    size_t id = idx3(i, j, k, state.NyT, state.NzT);
                    state.Ndip[id] = N0_density;
                    state.Nu[id] = inversion_fraction * N0_density;
                    state.Ng[id] = (1.0 - inversion_fraction) * N0_density;
                    state.Ng0[id] = N0_density;
                    total_volume += grid.dx[i] * grid.dy[j] * grid.dz[k];
                    ++n_cells;
                }
            }
        }

        std::cout << "[TLSRegion " << name << "] Initialized:\n";
        std::cout << "  Bounds: [" << i0 << "," << i1 << ") x [" << j0 << "," << j1 << ") x [" << k0 << "," << k1 << ")\n";
        std::cout << "  Cells: " << n_cells << ", Volume: " << total_volume * 1e18 << " μm³\n";
        std::cout << "  λ₀ = " << params.lambda0 * 1e9 << " nm, γ = " << params.gamma << " s⁻¹, τ = " << params.tau * 1e12 << " ps\n";
        std::cout << "  N₀ = " << N0_density << " m⁻³, Inversion = " << inversion_fraction << "\n";
    }
};

// Manages multiple TLS regions, each potentially with different parameters
class TLSRegionManager {
public:
    std::vector<TLSRegion> regions;
    size_t NxT = 0, NyT = 0, NzT = 0;
    bool enabled = false;
    std::vector<int> cell_region_map;  // -1 = no TLS, >= 0 = region index

    void add_region_from_structure(
        const GridSpacing& grid,
        real x_min, real x_max, real y_min, real y_max, real z_min, real z_max,
        const TwoLevelParams& params, real inversion_frac,
        const Shape* structure_shape = nullptr, const std::string& name = ""
    ) {
        TLSRegion region;
        region.region_id = regions.size();
        region.name = name.empty() ? "region_" + std::to_string(region.region_id) : name;

        region.i0 = grid.physical_to_index_x(x_min);
        region.i1 = grid.physical_to_index_x(x_max) + 1;
        region.j0 = grid.physical_to_index_y(y_min);
        region.j1 = grid.physical_to_index_y(y_max) + 1;
        region.k0 = grid.physical_to_index_z(z_min);
        region.k1 = grid.physical_to_index_z(z_max) + 1;

        region.params = params;
        region.inversion_fraction = inversion_frac;
        region.shape = structure_shape;

        regions.push_back(region);
        std::cout << "[TLSManager] Added region '" << region.name << "' at physical coords ["
                  << x_min * 1e9 << ", " << x_max * 1e9 << "] x ["
                  << y_min * 1e9 << ", " << y_max * 1e9 << "] x ["
                  << z_min * 1e9 << ", " << z_max * 1e9 << "] nm\n";
    }

    void initialize(TwoLevelState& state, const GridSpacing& grid) {
        if (regions.empty()) {
            enabled = false;
            return;
        }

        enabled = true;
        NxT = state.NxT; NyT = state.NyT; NzT = state.NzT;
        cell_region_map.assign(NxT * NyT * NzT, -1);

        for (size_t r = 0; r < regions.size(); ++r) {
            auto& region = regions[r];
            region.params.finalize();
            region.initialize(state, grid);

            for (size_t i = region.i0; i < region.i1; ++i) {
                for (size_t j = region.j0; j < region.j1; ++j) {
                    for (size_t k = region.k0; k < region.k1; ++k) {
                        if (!region.contains(i, j, k, grid)) continue;
                        cell_region_map[idx3(i, j, k, NyT, NzT)] = static_cast<int>(r);
                    }
                }
            }
        }
        std::cout << "[TLSManager] Initialized " << regions.size() << " TLS region(s)\n";
    }

    void finalize_with_dt(real dt) {
        for (auto& region : regions) {
            region.params.finalize_with_dt(dt);
        }
    }

    const TLSRegion* get_region_for_cell(size_t i, size_t j, size_t k) const {
        if (!enabled) return nullptr;
        size_t id = idx3(i, j, k, NyT, NzT);
        if (id >= cell_region_map.size()) return nullptr;
        int region_idx = cell_region_map[id];
        if (region_idx < 0 || region_idx >= (int)regions.size()) return nullptr;
        return &regions[region_idx];
    }

    TLSRegion* get_region(size_t idx) {
        return (idx < regions.size()) ? &regions[idx] : nullptr;
    }

    size_t num_regions() const { return regions.size(); }
    bool has_regions() const { return !regions.empty(); }
};

// Update polarization for multi-region TLS
inline void update_polarization_multi_region(
    size_t NxT, size_t NyT, size_t NzT,
    const TLSRegionManager& manager,
    const std::vector<real>& Ez,
    TwoLevelState& state
) {
    if (!manager.enabled) return;

#if FDTD_OMP_ENABLED
#pragma omp parallel for
#endif
    for (int i = 0; i < (int)NxT; ++i) {
        for (int j = 0; j < (int)NyT; ++j) {
            for (int k = 0; k < (int)NzT; ++k) {
                size_t id = idx3((size_t)i, (size_t)j, (size_t)k, NyT, NzT);
                if (state.Ndip[id] <= 0) continue;

                const TLSRegion* region = manager.get_region_for_cell(i, j, k);
                if (!region || !region->params.coefficients_initialized) continue;

                const auto& params = region->params;
                real Ng_frac = (state.Ng0[id] > 0) ?
                    (state.Ng[id] - state.Nu[id]) / state.Ng0[id] : 0.0;

                real driving_factor = state.Ndip[id] * Ng_frac * Ez[id];
                real Pz_new = params.kapa_coeff * driving_factor
                            + params.pa2 * state.Pz[id] + params.pa3 * state.Pz_prev[id];

                real Pz_old = state.Pz[id];
                state.Pz_prev[id] = Pz_old;
                state.Pz[id] = Pz_new;
                state.dPz_dt[id] = (Pz_new - Pz_old) * params.inv_dt;
            }
        }
    }
}

// Update populations for multi-region TLS
inline void update_populations_multi_region(
    size_t NxT, size_t NyT, size_t NzT, real dt,
    const TLSRegionManager& manager,
    const std::vector<real>& Ez,
    TwoLevelState& state
) {
    if (!manager.enabled) return;

#if FDTD_OMP_ENABLED
#pragma omp parallel for
#endif
    for (int i = 0; i < (int)NxT; ++i) {
        for (int j = 0; j < (int)NyT; ++j) {
            for (int k = 0; k < (int)NzT; ++k) {
                size_t id = idx3((size_t)i, (size_t)j, (size_t)k, NyT, NzT);
                if (state.Ndip[id] <= 0) continue;

                const TLSRegion* region = manager.get_region_for_cell(i, j, k);
                if (!region) continue;

                const auto& params = region->params;
                real Nu_curr = state.Nu[id];
                real Ntotal = state.Ng0[id];
                if (Ntotal <= 0) continue;

                real inv_tau = 1.0 / params.tau;
                real inv_hbar_omega = 1.0 / (params.hbar * params.omega_a);

                real E_avg = 0.5 * (Ez[id] + state.Ez_old[id]);
                real delta_P = state.Pz[id] - state.Pz_prev[id];
                real stim_rate = E_avg * delta_P * params.inv_dt * inv_hbar_omega;

                real dNu_dt = -inv_tau * Nu_curr + stim_rate;
                real dNg_dt = +inv_tau * Nu_curr - stim_rate;

                state.Nu[id] += dt * dNu_dt;
                state.Ng[id] += dt * dNg_dt;

                state.Nu[id] = std::max(0.0, state.Nu[id]);
                state.Ng[id] = std::max(0.0, state.Ng[id]);

                if (params.enable_population_clamp) {
                    real total = state.Nu[id] + state.Ng[id];
                    if (total > 1e-30) {
                        real scale = Ntotal / total;
                        state.Nu[id] *= scale;
                        state.Ng[id] *= scale;
                    }
                    state.Nu[id] = std::min(state.Nu[id], Ntotal);
                    state.Ng[id] = Ntotal - state.Nu[id];
                }
            }
        }
    }
}
