// two_level_system.hpp — Two-level atomic system for active gain media
// Implements laser physics: population inversion, stimulated emission, etc.

#pragma once

#include <vector>
#include <cmath>
#include <numbers>
#include <iostream>
#include "global_function.hpp"
#include "user_config.hpp"
#include "omp_config.hpp"

// Forward declaration
struct StructureTLSConfig;

// ========== Two-Level System Parameters ==========
struct TwoLevelParams {
    // Atomic transition parameters
    // These will be set from TLS_MATERIALS, defaults are just placeholders
    real lambda0 = 1500e-9;     // Transition wavelength (m)
    real omega_a = 0.0;         // Angular frequency ωa = 2πc/λ₀ (rad/s)

    // Decay and damping
    real gamma = 7e12;          // Polarization damping rate γ (s⁻¹)
    real tau = 1e-12;           // Upper level lifetime τ (s)

    // Transition dipole moment (calculated from Einstein A coefficient)
    real mu_z = 0.0;            // |μz| = √(3πε₀ℏc³/(ωa³τ))

    // Spatial distribution of dipoles
    real N0_total = 1e25;       // Total dipole density (m⁻³)

    // Physical constants (stored for use in stimulated term calculation)
    real hbar = PhysConst::HBAR;            // Reduced Planck constant (J·s)

    // Debug options (default from UserConfig)
    bool enable_population_clamp = UserConfig::TLS_ENABLE_CLAMP;

    // ========== Three-point recursion coefficients (computed in finalize_with_dt) ==========
    // For damped oscillator: d²P/dt² + γ·dP/dt + ω²·P = F
    // Discretized: P^{n+1} = kapa·F^n + pa2·P^n + pa3·P^{n-1}
    // Reference: Shijie_fdtd2dmetal4level_new.f90
    real kapa_coeff = 0.0;          // Driving term coefficient: 2·dt²/(2+γ·dt) × (2ω/ℏ)|μ|²
    real pa2 = 0.0;                 // P^n coefficient: 2·(2-ω²·dt²)/(2+γ·dt)
    real pa3 = 0.0;                 // P^{n-1} coefficient: (γ·dt-2)/(2+γ·dt)
    real inv_dt = 0.0;              // Pre-computed 1/dt (OPTIMIZATION: avoid divisions in loops)
    bool coefficients_initialized = false;

    // FIX #2: Az_value() now returns 0 (low-intensity approximation)
    // In low-intensity regime, the A² term in the polarization equation is negligible.
    // The original code incorrectly used mu_z (dipole moment) instead of vector potential.
    // Reference: oe-14-8-3569.pdf - at low intensity, omega_eff² ≈ omega_a²
    real Az_value() const { return 0.0; }

    // Initialize basic derived quantities (without dt)
    void finalize(real c0 = 299792458.0, real eps0 = 8.854187817e-12, real hbar_in = 1.054571817e-34) {
        hbar = hbar_in;
        omega_a = 2.0 * std::numbers::pi * c0 / lambda0;

        // Calculate transition dipole moment from Einstein A coefficient
        // |μz| = √(3πε₀ℏc³/(ωa³τ))
        real numerator = 3.0 * std::numbers::pi * eps0 * hbar * c0 * c0 * c0;
        real denominator = omega_a * omega_a * omega_a * tau;
        mu_z = std::sqrt(numerator / denominator);

        // CRITICAL FIX (Issue D): Print warnings for potentially dangerous parameters
        // Small tau leads to large mu_z, causing extremely strong coupling → numerical instability
        // mu_z ~ 1/sqrt(tau), so tau=1e-14 gives mu_z ~100x larger than tau=1e-10
        std::cout << "\n=== Two-Level System Parameter Check ===\n";
        std::cout << "  tau (upper level lifetime) = " << tau << " s\n";
        std::cout << "  mu_z (dipole moment)       = " << mu_z << " C·m\n";
        std::cout << "  omega_a (angular freq)     = " << omega_a << " rad/s\n";
        std::cout << "  gamma (damping rate)       = " << gamma << " s⁻¹\n";

        // Dimensionless coupling strength indicator: mu_z * sqrt(N0 / (eps0 * hbar * omega_a))
        real coupling_indicator = mu_z * std::sqrt(N0_total / (eps0 * hbar * omega_a));
        std::cout << "  Coupling strength indicator = " << coupling_indicator << "\n";

        // Warn if tau is dangerously small (typically should be 1e-12 to 1e-9 s)
        if (tau < 1e-12) {
            std::cerr << "\n*** WARNING (Issue D): tau = " << tau << " s is very small! ***\n";
            std::cerr << "    This leads to extremely large mu_z = " << mu_z << " C·m\n";
            std::cerr << "    Strong coupling may cause numerical instability/divergence.\n";
            std::cerr << "    Recommended range: tau = 1e-12 to 1e-9 s\n";
            std::cerr << "    Consider increasing tau or reducing time step dt.\n\n";
        }

        // Warn if coupling is extremely strong
        if (coupling_indicator > 1e6) {
            std::cerr << "*** WARNING: Coupling indicator = " << coupling_indicator << " is very large! ***\n";
            std::cerr << "    This may cause rapid energy exchange and numerical issues.\n\n";
        }

        std::cout << "========================================\n\n";
    }

    // Initialize three-point recursion coefficients (requires dt from simulation)
    // Must be called AFTER finalize() and after dt is determined
    void finalize_with_dt(real dt) {
        if (omega_a <= 0.0) {
            std::cerr << "ERROR: finalize_with_dt called before finalize()!\n";
            return;
        }

        const real dt2 = dt * dt;
        const real omega_sq = omega_a * omega_a;

        // Denominator: β = 2 + γ·dt
        const real beta = 2.0 + gamma * dt;

        // Three-point recursion coefficients for damped harmonic oscillator:
        // d²P/dt² + γ·dP/dt + ω²·P = F
        // Using central difference discretization:
        // (P^{n+1} - 2P^n + P^{n-1})/dt² + γ·(P^{n+1} - P^{n-1})/(2dt) + ω²·P^n = F^n
        // Rearranging:
        // P^{n+1}·(1 + γ·dt/2) = F^n·dt² + P^n·(2 - ω²·dt²) + P^{n-1}·(γ·dt/2 - 1)
        // Multiply through by 2:
        // P^{n+1}·(2 + γ·dt) = 2·F^n·dt² + P^n·2·(2 - ω²·dt²) + P^{n-1}·(γ·dt - 2)

        // Driving coefficient: 2·dt²/β × (2ω/ℏ)|μ|²
        // Note: Full driving term = kapa_coeff × Ndip × (Ng-Nu)/Ng0 × Ez
        const real two_omega_over_hbar = 2.0 * omega_a / hbar;
        const real mu_z_sq = mu_z * mu_z;
        kapa_coeff = (2.0 * dt2 / beta) * two_omega_over_hbar * mu_z_sq;

        // P^n coefficient: 2·(2 - ω²·dt²)/β
        pa2 = 2.0 * (2.0 - omega_sq * dt2) / beta;

        // P^{n-1} coefficient: (γ·dt - 2)/β
        pa3 = (gamma * dt - 2.0) / beta;

        // Pre-compute 1/dt (OPTIMIZATION: avoid divisions in inner loops)
        inv_dt = 1.0 / dt;

        coefficients_initialized = true;

        std::cout << "=== Three-point recursion coefficients ===\n";
        std::cout << "  dt = " << dt * 1e15 << " fs\n";
        std::cout << "  beta (2+γ·dt) = " << beta << "\n";
        std::cout << "  kapa_coeff = " << kapa_coeff << "\n";
        std::cout << "  pa2 = " << pa2 << "\n";
        std::cout << "  pa3 = " << pa3 << "\n";

        // Stability check for damped harmonic oscillator
        // For explicit central difference scheme:
        // - Undamped (γ=0): requires ω·dt < 2
        // - Damped (γ>0): stable if ω·dt < 2 (damping helps stability)
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

// ========== Two-Level System State ==========
struct TwoLevelState {
    size_t NxT, NyT, NzT;

    // Gain region bounds (half-open intervals) - stored for external access
    size_t gain_i0 = 0, gain_i1 = 0;
    size_t gain_j0 = 0, gain_j1 = 0;
    size_t gain_k0 = 0, gain_k1 = 0;

    // Population densities [m^-3] (VOLUME DENSITY FORMULATION)
    // All population quantities are densities, NOT per-cell counts
    // This ensures physics is independent of grid resolution (dx, dy, dz)
    std::vector<real> Ng;           // Ground state density Ng(r⃗,t) [m^-3]
    std::vector<real> Nu;           // Upper state density Nu(r⃗,t) [m^-3]
    std::vector<real> Ng0;          // Total density N_total(r⃗) = Ng + Nu [m^-3] (conserved)

    // Polarization (C/m²) - Three-point recursion storage
    std::vector<real> Pz;           // Current polarization Pz(r⃗,t) = P^n
    std::vector<real> Pz_prev;      // Previous time step Pz(r⃗,t-dt) = P^{n-1}
    std::vector<real> dPz_dt;       // Time derivative dPz/dt (for E-field update)

    // FIX #4: Store previous Ez for stimulated term calculation
    // stim = E_avg * ΔP / (ℏωa), where E_avg = 0.5*(Ez_new + Ez_old)
    std::vector<real> Ez_old;       // Previous time step Ez for stim calculation

    // Dipole density distribution Ndip(r⃗) [m^-3]
    std::vector<real> Ndip;         // Total dipole density, nonzero only in gain region

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
    
    // Initialize: set up gain region with population inversion
    // CRITICAL FIX (Issue C): Use half-open interval [i0, i1) to avoid off-by-one
    // Previous: for (i = i0; i <= i1) gave size = (i1-i0+1), one extra cell
    // Now: for (i = i0; i < i1) gives size = (i1-i0), exactly as specified
    //
    // VOLUME DENSITY FORMULATION:
    // All population quantities (Ndip, Nu, Ng, Ng0) are now VOLUME DENSITIES [m^-3]
    // This ensures grid-independent physics: changing resolution (dx,dy,dz) does NOT
    // change macroscopic behavior (gain, threshold, etc.)
    void initialize_gain_region(
        size_t i0, size_t i1,   // Half-open interval [i0, i1)
        size_t j0, size_t j1,   // Half-open interval [j0, j1)
        size_t k0, size_t k1,   // Half-open interval [k0, k1)
        real N0_density,        // Dipole density [atoms/m³]
        const GridSpacing& grid, // Grid spacing for computing total volume
        real inversion_fraction = 0.1  // Initial population in upper state
    ) {
        // Boundary validation: ensure indices are within grid bounds
        i1 = std::min(i1, NxT);
        j1 = std::min(j1, NyT);
        k1 = std::min(k1, NzT);

        // Ensure i0 <= i1 (no underflow or invalid range)
        if (i0 >= i1 || j0 >= j1 || k0 >= k1) {
            std::cerr << "Warning: Invalid gain region bounds, skipping initialization\n";
            return;
        }

        // Store validated gain region bounds for external access
        gain_i0 = i0; gain_i1 = i1;
        gain_j0 = j0; gain_j1 = j1;
        gain_k0 = k0; gain_k1 = k1;

        // Statistics for debugging
        real total_volume = 0.0;
        real min_cell_vol = 1e30, max_cell_vol = 0.0;

        for (size_t i = i0; i < i1; ++i) {       // Half-open: [i0, i1)
            for (size_t j = j0; j < j1; ++j) {   // Half-open: [j0, j1)
                for (size_t k = k0; k < k1; ++k) { // Half-open: [k0, k1)
                    size_t id = idx3(i, j, k, NyT, NzT);

                    // VOLUME DENSITY FORMULATION:
                    // Ndip, Nu, Ng, Ng0 are all densities [m^-3], NOT atoms per cell
                    // This makes the physics grid-independent
                    Ndip[id] = N0_density;                          // [m^-3]
                    Nu[id]   = inversion_fraction * N0_density;     // [m^-3]
                    Ng[id]   = (1.0 - inversion_fraction) * N0_density;  // [m^-3]
                    Ng0[id]  = N0_density;  // Total density (conserved) [m^-3]

                    // Track volume statistics
                    real cell_volume = grid.dx[i] * grid.dy[j] * grid.dz[k];
                    total_volume += cell_volume;
                    min_cell_vol = std::min(min_cell_vol, cell_volume);
                    max_cell_vol = std::max(max_cell_vol, cell_volume);
                }
            }
        }

        // Compute equivalent total atom count for reference
        real total_atoms = N0_density * total_volume;
        real inversion_density = (2.0 * inversion_fraction - 1.0) * N0_density;  // [m^-3]
        real integrated_inversion = inversion_density * total_volume;

        // Print initialization summary
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

// ========== Update Functions ==========

// Update polarization using THREE-POINT RECURSION (higher accuracy than Verlet)
//
// Physics: Damped harmonic oscillator ODE
// d²Pz/dt² + γ·dPz/dt + ωa²·Pz = (2ωa/ℏ)·|μz|²·Ndip·[(Ng-Nu)/Ng⁰]·Ez
//
// Discretization using central difference (Reference: Shijie_fdtd2dmetal4level_new.f90):
// (P^{n+1} - 2P^n + P^{n-1})/dt² + γ·(P^{n+1} - P^{n-1})/(2dt) + ω²·P^n = F^n
//
// Rearranging to explicit form:
// P^{n+1} = kapa·F^n + pa2·P^n + pa3·P^{n-1}
//
// where:
//   kapa = 2·dt²/(2+γ·dt) × driving_coefficient
//   pa2 = 2·(2-ω²·dt²)/(2+γ·dt)
//   pa3 = (γ·dt-2)/(2+γ·dt)
//
// This three-point recursion is more accurate than Verlet for damped oscillators
// because it exactly captures the damping term without first-order error.

inline void update_polarization(
    size_t NxT, size_t NyT, size_t NzT,
    const TwoLevelParams& params,
    const std::vector<real>& Ez,  // Electric field Ez
    TwoLevelState& state
) {
    // Ensure coefficients are initialized
    if (!params.coefficients_initialized) {
        std::cerr << "ERROR: TwoLevelParams::finalize_with_dt() not called!\n";
        return;
    }

    const real kapa_coeff = params.kapa_coeff;
    const real pa2 = params.pa2;
    const real pa3 = params.pa3;
    const real inv_dt = params.inv_dt;  // Pre-computed 1/dt

#if FDTD_OMP_ENABLED
#pragma omp parallel for
#endif
    for (int i = 0; i < (int)NxT; ++i) {
        for (int j = 0; j < (int)NyT; ++j) {
            for (int k = 0; k < (int)NzT; ++k) {
                size_t id = idx3((size_t)i, (size_t)j, (size_t)k, NyT, NzT);

                // Skip if no dipoles here (Ndip is now density [m^-3])
                if (state.Ndip[id] <= 0) continue;  // Skip cells with no gain medium

                // Calculate population inversion term: (Ng - Nu) / Ng0
                // Positive when more atoms in ground state (absorption)
                // Negative when more atoms in upper state (gain/amplification)
                real Ng_frac = (state.Ng0[id] > 0) ?
                    (state.Ng[id] - state.Nu[id]) / state.Ng0[id] : 0.0;

                // Driving force coefficient: Ndip·[(Ng-Nu)/Ng⁰]·Ez
                // Full driving = kapa_coeff × driving_factor
                real driving_factor = state.Ndip[id] * Ng_frac * Ez[id];

                // Three-point recursion: P^{n+1} = kapa·F + pa2·P^n + pa3·P^{n-1}
                real Pz_new = kapa_coeff * driving_factor
                            + pa2 * state.Pz[id]
                            + pa3 * state.Pz_prev[id];

                // Calculate dPz/dt at half-step n+1/2 for E-field update
                // dPz/dt|_{n+1/2} = (P^{n+1} - P^n) / dt
                real Pz_old = state.Pz[id];
                real dPz_dt_half = (Pz_new - Pz_old) * inv_dt;  // Use pre-computed 1/dt

                // Shift polarization history: P^{n-1} ← P^n, P^n ← P^{n+1}
                state.Pz_prev[id] = Pz_old;
                state.Pz[id] = Pz_new;
                state.dPz_dt[id] = dPz_dt_half;
            }
        }
    }
}

// Update population densities using CORRECTED two-level rate equations:
//
// Rate equations (VOLUME DENSITY FORMULATION, all quantities in m^-3):
//   stim_rate = E_avg · (dP/dt) / (ℏωa)   [m^-3 s^-1]
//   dNu/dt = -Nu/τ + stim_rate
//   dNg/dt = +Nu/τ - stim_rate
//
// where:
//   E_avg = 0.5 * (Ez_new + Ez_old)   [V/m, time-centered field]
//   dP/dt = (Pz_new - Pz_old) / dt    [C/(m²·s)]
//
// Sign convention:
//   stim_rate > 0: field does positive work → ABSORPTION → Nu increases
//   stim_rate < 0: field does negative work → STIMULATED EMISSION → Nu decreases
//
// Conservation properties:
//   1. Nu + Ng = Ntotal is conserved (no external pumping/losses)
//   2. Energy conservation: field energy change = -ℏωa × population change
//   3. Positive inversion (Nu > Ng) causes amplification

inline void update_populations(
    size_t NxT, size_t NyT, size_t NzT,
    real dt,
    const TwoLevelParams& params,
    const std::vector<real>& Ez,  // Current (new) electric field
    TwoLevelState& state
) {
    const real tau = params.tau;
    const real hbar = params.hbar;
    const real omega_a = params.omega_a;
    const real inv_dt = params.inv_dt;  // Pre-computed 1/dt

    const real inv_tau = 1.0 / tau;
    // Coefficient for stimulated term: 1/(ℏωa)
    const real inv_hbar_omega = 1.0 / (hbar * omega_a);

#if FDTD_OMP_ENABLED
#pragma omp parallel for
#endif
    for (int i = 0; i < (int)NxT; ++i) {
        for (int j = 0; j < (int)NyT; ++j) {
            for (int k = 0; k < (int)NzT; ++k) {
                size_t id = idx3((size_t)i, (size_t)j, (size_t)k, NyT, NzT);

                // Skip if no dipoles here (Ndip is now density [m^-3])
                if (state.Ndip[id] <= 0) continue;  // Skip cells with no gain medium

                real Nu_curr = state.Nu[id];
                real Ntotal = state.Ng0[id];  // Total density [m^-3] (conserved)

                if (Ntotal <= 0) continue;  // Skip cells with no atoms

                // Stimulated term using energy-conserving RATE form
                // stim_rate = E_avg · (dP/dt) / (ℏωa)  [units: m^-3 s^-1]
                // E_avg = 0.5 * (Ez_new + Ez_old)   [V/m]
                // dP/dt = ΔP/dt where ΔP = Pz_new - Pz_old  [C/(m²·s)]
                //
                // Dimensional analysis:
                //   E·(dP/dt)/(ℏω) = [V/m]·[C/(m²·s)]/[J·s·s^-1]
                //                  = [J/(m³·s)]/[J] = [m^-3 s^-1] ✓
                //
                // Sign convention:
                //   stim_rate > 0: field does positive work → ABSORPTION → Nu increases
                //   stim_rate < 0: field does negative work → STIMULATED EMISSION → Nu decreases
                real E_avg = 0.5 * (Ez[id] + state.Ez_old[id]);
                real delta_P = state.Pz[id] - state.Pz_prev[id];
                real dP_dt = delta_P * inv_dt;  // [C/(m²·s)]
                real stim_rate = E_avg * dP_dt * inv_hbar_omega;  // [m^-3 s^-1]

                // FIX #3: Simple linear spontaneous decay: -Nu/τ
                // (Original code had -Nu·(1-Ng/Ng0)/τ = -Nu²/(τ·Ng0), which was WRONG)
                //
                // FIX #5: Stimulated term sign correction
                // Physical meaning of stim_rate = E·(dP/dt)/(ℏω):
                //   - stim_rate > 0: field does positive work on medium → ABSORPTION → Nu increases
                //   - stim_rate < 0: field does negative work on medium → STIMULATED EMISSION → Nu decreases
                // Therefore: dNu/dt = -Nu/τ + stim_rate (NOT minus!)
                real dNu_dt = -inv_tau * Nu_curr + stim_rate;
                real dNg_dt = +inv_tau * Nu_curr - stim_rate;

                // Forward Euler update
                state.Nu[id] += dt * dNu_dt;
                state.Ng[id] += dt * dNg_dt;

                // Basic protection: ensure populations stay non-negative
                state.Nu[id] = std::max(0.0, state.Nu[id]);
                state.Ng[id] = std::max(0.0, state.Ng[id]);

                // Optional debug clamp: enforce Nu + Ng = Ntotal
                if (params.enable_population_clamp) {
                    real total = state.Nu[id] + state.Ng[id];
                    if (total > 1e-30) {
                        real scale = Ntotal / total;
                        state.Nu[id] *= scale;
                        state.Ng[id] *= scale;
                    }
                    // Additional clamp to valid range
                    state.Nu[id] = std::min(state.Nu[id], Ntotal);
                    state.Ng[id] = Ntotal - state.Nu[id];
                }
            }
        }
    }
}

// Modified E-field update to include polarization source term:
// ∂E/∂t = (1/ε₀n²)[∇×H - ∂P/∂t - J]
//
// FIX #1: Polarization source term sign correction
// The Maxwell-Ampere equation requires: Ez += bEz * (curlHz - Jz - dPz_dt)
// The "-dPz_dt" term represents energy transfer between field and medium.
//
// Original INCORRECT code:
//   Pz_source = -dPz_dt;
//   Ez = ... - Pz_source;  // This gave +dPz_dt (WRONG)
//
// CORRECT code:
//   Pz_source = dPz_dt;
//   Ez = ... - Pz_source;  // This gives -dPz_dt (CORRECT)
//
// Reference: oe-14-8-3569.pdf, Eq. for Maxwell-Ampere law

template<typename Real>
inline void fdtd_update_E_with_gain(
    size_t NxT, size_t NyT, size_t NzT,
    const Real* __restrict dx_array,
    const Real* __restrict dy_array,
    const Real* __restrict dz_array,
    // Pre-computed inverse arrays (OPTIMIZATION: avoid division in inner loops)
    const Real* __restrict inv_dx_array,
    const Real* __restrict inv_dy_array,
    const Real* __restrict inv_dz_array,
    const Real* __restrict aEx, const Real* __restrict bEx,
    const Real* __restrict aEy, const Real* __restrict bEy,
    const Real* __restrict aEz, const Real* __restrict bEz,
    Real* __restrict Ex, Real* __restrict Ey, Real* __restrict Ez,
    Real* __restrict Hx, Real* __restrict Hy, Real* __restrict Hz,
    const Real* __restrict Jx, const Real* __restrict Jy, const Real* __restrict Jz,
    TwoLevelState& state  // Non-const: we update Ez_old
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

                // Use pre-computed inverses (eliminates O(N³) divisions)
                const Real inv_dx = inv_dx_array[i];
                const Real inv_dy = inv_dy_array[j];
                const Real inv_dz = inv_dz_array[k];

                // curl H calculation (same as before)
                Real curlHx = diff_ym(Hz, id, sJ, inv_dy) - diff_zm(Hy, id, sK, inv_dz);
                Real curlHy = diff_zm(Hx, id, sK, inv_dz) - diff_xm(Hz, id, sI, inv_dx);
                Real curlHz = diff_xm(Hy, id, sI, inv_dx) - diff_ym(Hx, id, sJ, inv_dy);

                // Standard E updates for x and y components
                Ex[id] = aEx[id] * Ex[id] + bEx[id] * (curlHx - Jx[id]);
                Ey[id] = aEy[id] * Ey[id] + bEy[id] * (curlHy - Jy[id]);

                // Store Ez_old BEFORE updating Ez (for next step's stim calculation)
                state.Ez_old[id] = Ez[id];

                // FIX #1: Polarization source term with CORRECT sign
                // The term is -∂P/∂t in Maxwell-Ampere equation
                // dPz_dt holds the positive value of ∂Pz/∂t
                // So we need to SUBTRACT it: Ez += bEz*(curlHz - Jz - dPz_dt)
                Real Pz_source = state.dPz_dt[id];  // This is +dP/dt
                Ez[id] = aEz[id] * Ez[id] + bEz[id] * (curlHz - Jz[id] - Pz_source);
            }
        }
    }
}

// ========== Diagnostic Functions ==========
// VOLUME DENSITY FORMULATION: All population quantities are densities [m^-3]
// To get physically meaningful totals, must integrate with cell volume dV

// Calculate integrated total population in upper state: ∫ Nu dV [atoms]
inline real compute_integrated_Nu(const TwoLevelState& state, const GridSpacing& grid) {
    real total_Nu = 0.0;
    const size_t NyT = state.NyT;
    const size_t NzT = state.NzT;

    for (size_t i = 0; i < state.NxT; ++i) {
        for (size_t j = 0; j < state.NyT; ++j) {
            for (size_t k = 0; k < state.NzT; ++k) {
                size_t id = idx3(i, j, k, NyT, NzT);
                if (state.Ndip[id] <= 0) continue;

                real dV = grid.dx[i] * grid.dy[j] * grid.dz[k];
                total_Nu += state.Nu[id] * dV;
            }
        }
    }
    return total_Nu;
}

// Calculate integrated total population in ground state: ∫ Ng dV [atoms]
inline real compute_integrated_Ng(const TwoLevelState& state, const GridSpacing& grid) {
    real total_Ng = 0.0;
    const size_t NyT = state.NyT;
    const size_t NzT = state.NzT;

    for (size_t i = 0; i < state.NxT; ++i) {
        for (size_t j = 0; j < state.NyT; ++j) {
            for (size_t k = 0; k < state.NzT; ++k) {
                size_t id = idx3(i, j, k, NyT, NzT);
                if (state.Ndip[id] <= 0) continue;

                real dV = grid.dx[i] * grid.dy[j] * grid.dz[k];
                total_Ng += state.Ng[id] * dV;
            }
        }
    }
    return total_Ng;
}

// Calculate integrated population inversion: ∫ (Nu - Ng) dV [atoms]
inline real compute_integrated_inversion(const TwoLevelState& state, const GridSpacing& grid) {
    real total_inversion = 0.0;
    const size_t NyT = state.NyT;
    const size_t NzT = state.NzT;

    for (size_t i = 0; i < state.NxT; ++i) {
        for (size_t j = 0; j < state.NyT; ++j) {
            for (size_t k = 0; k < state.NzT; ++k) {
                size_t id = idx3(i, j, k, NyT, NzT);
                if (state.Ndip[id] <= 0) continue;

                real dV = grid.dx[i] * grid.dy[j] * grid.dz[k];
                total_inversion += (state.Nu[id] - state.Ng[id]) * dV;
            }
        }
    }
    return total_inversion;
}

// Compute average inversion density in gain region: (∫(Nu-Ng)dV) / (∫dV) [m^-3]
inline real compute_avg_inversion_density(const TwoLevelState& state, const GridSpacing& grid) {
    real total_inversion = 0.0;
    real total_volume = 0.0;
    const size_t NyT = state.NyT;
    const size_t NzT = state.NzT;

    for (size_t i = 0; i < state.NxT; ++i) {
        for (size_t j = 0; j < state.NyT; ++j) {
            for (size_t k = 0; k < state.NzT; ++k) {
                size_t id = idx3(i, j, k, NyT, NzT);
                if (state.Ndip[id] <= 0) continue;

                real dV = grid.dx[i] * grid.dy[j] * grid.dz[k];
                total_inversion += (state.Nu[id] - state.Ng[id]) * dV;
                total_volume += dV;
            }
        }
    }

    return (total_volume > 0) ? total_inversion / total_volume : 0.0;
}

// Calculate dipole-field interaction energy: U_int = -∫ P·E dV
// For z-polarized case: U_int = -∫ Pz·Ez dV
// We use 0.5 factor for time-averaged energy: u_int = -0.5 * Pz * Ez
//
// Physical interpretation:
//   u_int < 0 when P and E are aligned → stable configuration
//   u_int > 0 when P and E are anti-aligned → unstable
//
// CRITICAL FIX: The sign must be NEGATIVE (dipole interaction u_int = -P·E)
inline real compute_polarization_energy(
    const TwoLevelState& state,
    const std::vector<real>& Ez,
    const GridSpacing& grid
) {
    real total_energy = 0.0;
    size_t NyT = state.NyT;
    size_t NzT = state.NzT;

    for (size_t i = 0; i < state.NxT; ++i) {
        for (size_t j = 0; j < state.NyT; ++j) {
            for (size_t k = 0; k < state.NzT; ++k) {
                size_t id = idx3(i, j, k, NyT, NzT);

                if (state.Ndip[id] <= 0) continue;  // Skip cells with no gain medium

                real dV = grid.dx[i] * grid.dy[j] * grid.dz[k];

                // Dipole interaction energy density: u_int = -P·E
                // With time-averaging factor: -0.5 * Pz * Ez
                real energy_density = -0.5 * state.Pz[id] * Ez[id];
                total_energy += energy_density * dV;
            }
        }
    }
    return total_energy;
}

// ============================================================================
//                    MULTI-REGION TLS SUPPORT (Structure-Bound TLS)
// ============================================================================
//
// This section provides support for binding TLS to structures.
// Each structure can have its own TLS region with independent parameters.
// The TLS regions are managed together but each maintains its own physics.

// ========== Single TLS Region ==========
// Represents a TLS region bound to a specific structure
struct TLSRegion {
    // Region identifier
    size_t region_id = 0;
    std::string name = "region_0";

    // Bounds in grid indices (half-open intervals [i0, i1))
    size_t i0 = 0, i1 = 0;
    size_t j0 = 0, j1 = 0;
    size_t k0 = 0, k1 = 0;

    // TLS parameters for this region
    TwoLevelParams params;

    // Check if a cell index is within this region
    bool contains(size_t i, size_t j, size_t k) const {
        return (i >= i0 && i < i1 &&
                j >= j0 && j < j1 &&
                k >= k0 && k < k1);
    }

    // Get the bounding box in physical coordinates
    void get_physical_bounds(const GridSpacing& grid,
                             real& x_min, real& x_max,
                             real& y_min, real& y_max,
                             real& z_min, real& z_max) const {
        x_min = grid.x_bounds[i0];
        x_max = grid.x_bounds[i1];
        y_min = grid.y_bounds[j0];
        y_max = grid.y_bounds[j1];
        z_min = grid.z_bounds[k0];
        z_max = grid.z_bounds[k1];
    }

    // Initial population inversion fraction for this region
    // This is set from the structure's TLS config
    real inversion_fraction = 1.0;

    // Initialize gain region for this TLS region
    void initialize(TwoLevelState& state, const GridSpacing& grid) {
        // Validate bounds
        if (i0 >= i1 || j0 >= j1 || k0 >= k1) {
            std::cerr << "[TLSRegion " << name << "] Warning: Invalid bounds, skipping\n";
            return;
        }

        size_t NyT = state.NyT;
        size_t NzT = state.NzT;

        real total_volume = 0.0;
        real N0_density = params.N0_total;

        // Initialize cells in this region
        for (size_t i = i0; i < i1; ++i) {
            for (size_t j = j0; j < j1; ++j) {
                for (size_t k = k0; k < k1; ++k) {
                    size_t id = idx3(i, j, k, NyT, NzT);

                    state.Ndip[id] = N0_density;
                    state.Nu[id] = inversion_fraction * N0_density;
                    state.Ng[id] = (1.0 - inversion_fraction) * N0_density;
                    state.Ng0[id] = N0_density;

                    real cell_volume = grid.dx[i] * grid.dy[j] * grid.dz[k];
                    total_volume += cell_volume;
                }
            }
        }

        size_t n_cells = (i1 - i0) * (j1 - j0) * (k1 - k0);
        std::cout << "[TLSRegion " << name << "] Initialized:\n";
        std::cout << "  Bounds: [" << i0 << "," << i1 << ") x [" << j0 << "," << j1 << ") x [" << k0 << "," << k1 << ")\n";
        std::cout << "  Cells: " << n_cells << ", Volume: " << total_volume * 1e18 << " μm³\n";
        std::cout << "  λ₀ = " << params.lambda0 * 1e9 << " nm, γ = " << params.gamma << " s⁻¹, τ = " << params.tau * 1e12 << " ps\n";
        std::cout << "  N₀ = " << N0_density << " m⁻³, Inversion = " << inversion_fraction << "\n";
    }
};

// ========== Multi-Region TLS Manager ==========
// Manages multiple TLS regions, each potentially with different parameters
class TLSRegionManager {
public:
    // All TLS regions
    std::vector<TLSRegion> regions;

    // Grid dimensions (stored for convenience)
    size_t NxT = 0, NyT = 0, NzT = 0;

    // Global enable flag
    bool enabled = false;

    // Index mapping: for each cell, which region does it belong to?
    // -1 means no TLS region, >= 0 is the region index
    std::vector<int> cell_region_map;

    // Add a TLS region from structure bounds
    // Converts physical coordinates to grid indices
    void add_region_from_structure(
        const GridSpacing& grid,
        real x_min, real x_max,
        real y_min, real y_max,
        real z_min, real z_max,
        const TwoLevelParams& params,
        real inversion_frac,
        const std::string& name = ""
    ) {
        TLSRegion region;
        region.region_id = regions.size();
        region.name = name.empty() ? "region_" + std::to_string(region.region_id) : name;

        // Convert physical coordinates to grid indices
        // Note: grid coordinates are relative to core origin (after PML)
        region.i0 = grid.physical_to_index_x(x_min);
        region.i1 = grid.physical_to_index_x(x_max) + 1;
        region.j0 = grid.physical_to_index_y(y_min);
        region.j1 = grid.physical_to_index_y(y_max) + 1;
        region.k0 = grid.physical_to_index_z(z_min);
        region.k1 = grid.physical_to_index_z(z_max) + 1;

        region.params = params;
        region.inversion_fraction = inversion_frac;

        regions.push_back(region);
        std::cout << "[TLSManager] Added region '" << region.name << "' at physical coords ["
                  << x_min * 1e9 << ", " << x_max * 1e9 << "] x ["
                  << y_min * 1e9 << ", " << y_max * 1e9 << "] x ["
                  << z_min * 1e9 << ", " << z_max * 1e9 << "] nm\n";
    }

    // Initialize all regions
    void initialize(TwoLevelState& state, const GridSpacing& grid) {
        if (regions.empty()) {
            enabled = false;
            return;
        }

        enabled = true;
        NxT = state.NxT;
        NyT = state.NyT;
        NzT = state.NzT;

        // Build cell-to-region map
        cell_region_map.assign(NxT * NyT * NzT, -1);

        for (size_t r = 0; r < regions.size(); ++r) {
            auto& region = regions[r];

            // Finalize parameters
            region.params.finalize();

            // Initialize this region's cells
            region.initialize(state, grid);

            // Update cell-to-region map
            for (size_t i = region.i0; i < region.i1; ++i) {
                for (size_t j = region.j0; j < region.j1; ++j) {
                    for (size_t k = region.k0; k < region.k1; ++k) {
                        size_t id = idx3(i, j, k, NyT, NzT);
                        cell_region_map[id] = static_cast<int>(r);
                    }
                }
            }
        }

        std::cout << "[TLSManager] Initialized " << regions.size() << " TLS region(s)\n";
    }

    // Finalize all regions with dt (must be called after make_context)
    void finalize_with_dt(real dt) {
        for (auto& region : regions) {
            region.params.finalize_with_dt(dt);
        }
    }

    // Get the region for a cell (returns nullptr if cell has no TLS)
    const TLSRegion* get_region_for_cell(size_t i, size_t j, size_t k) const {
        if (!enabled) return nullptr;
        size_t id = idx3(i, j, k, NyT, NzT);
        if (id >= cell_region_map.size()) return nullptr;
        int region_idx = cell_region_map[id];
        if (region_idx < 0 || region_idx >= (int)regions.size()) return nullptr;
        return &regions[region_idx];
    }

    // Get region by index
    TLSRegion* get_region(size_t idx) {
        if (idx < regions.size()) return &regions[idx];
        return nullptr;
    }

    // Number of regions
    size_t num_regions() const { return regions.size(); }

    // Check if any region exists
    bool has_regions() const { return !regions.empty(); }
};

// ========== Multi-Region Update Functions ==========

// Update polarization for multi-region TLS
// Each cell uses the parameters from its assigned region
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

                // Skip if no dipoles or no region
                if (state.Ndip[id] <= 0) continue;

                // Get the region for this cell
                const TLSRegion* region = manager.get_region_for_cell(i, j, k);
                if (!region) continue;

                const auto& params = region->params;
                if (!params.coefficients_initialized) continue;

                // Three-point recursion with region-specific parameters
                real Ng_frac = (state.Ng0[id] > 0) ?
                    (state.Ng[id] - state.Nu[id]) / state.Ng0[id] : 0.0;

                real driving_factor = state.Ndip[id] * Ng_frac * Ez[id];

                real Pz_new = params.kapa_coeff * driving_factor
                            + params.pa2 * state.Pz[id]
                            + params.pa3 * state.Pz_prev[id];

                real Pz_old = state.Pz[id];
                real dPz_dt_half = (Pz_new - Pz_old) * params.inv_dt;

                state.Pz_prev[id] = Pz_old;
                state.Pz[id] = Pz_new;
                state.dPz_dt[id] = dPz_dt_half;
            }
        }
    }
}

// Update populations for multi-region TLS
inline void update_populations_multi_region(
    size_t NxT, size_t NyT, size_t NzT,
    real dt,
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
                real dP_dt = delta_P * params.inv_dt;
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