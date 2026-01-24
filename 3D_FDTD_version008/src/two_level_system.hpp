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

// ========== Two-Level System Parameters ==========
struct TwoLevelParams {
    // Atomic transition parameters (defaults from UserConfig)
    real lambda0 = UserConfig::TLS_LAMBDA0;     // Transition wavelength (m)
    real omega_a = 0.0;                         // Angular frequency ωa = 2πc/λ₀ (rad/s)

    // Decay and damping (defaults from UserConfig)
    real gamma = UserConfig::TLS_GAMMA;     // Polarization damping rate γ (s⁻¹)
    real tau = UserConfig::TLS_TAU;         // Upper level lifetime τ (s)

    // Transition dipole moment (calculated from Einstein A coefficient)
    real mu_z = 0.0;                        // |μz| = √(3πε₀ℏc³/(ωa³τ))

    // Spatial distribution of dipoles (default from UserConfig)
    real N0_total = UserConfig::TLS_N0_TOTAL;   // Total dipole density (m⁻³)

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

    // Population densities (m⁻³)
    std::vector<real> Ng;           // Ground state density Ng(r⃗,t)
    std::vector<real> Nu;           // Upper state density Nu(r⃗,t)
    std::vector<real> Ng0;          // Initial ground state density Ng⁰(r⃗) = total density

    // Polarization (C/m²) - Three-point recursion storage
    std::vector<real> Pz;           // Current polarization Pz(r⃗,t) = P^n
    std::vector<real> Pz_prev;      // Previous time step Pz(r⃗,t-dt) = P^{n-1}
    std::vector<real> dPz_dt;       // Time derivative dPz/dt (for E-field update)

    // FIX #4: Store previous Ez for stimulated term calculation
    // stim = E_avg * ΔP / (ℏωa), where E_avg = 0.5*(Ez_new + Ez_old)
    std::vector<real> Ez_old;       // Previous time step Ez for stim calculation

    // Dipole density distribution Ndip(r⃗)
    std::vector<real> Ndip;         // Only nonzero where gain medium exists

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
    // CRITICAL FIX (Issue F): N0_density is atoms/m³, must multiply by cell volume
    // Previous: Ndip[id] = N0_total (treating density as per-cell count - WRONG!)
    // Now: Ndip[id] = N0_density * cell_volume (correct atoms per cell)
    void initialize_gain_region(
        size_t i0, size_t i1,   // Half-open interval [i0, i1)
        size_t j0, size_t j1,   // Half-open interval [j0, j1)
        size_t k0, size_t k1,   // Half-open interval [k0, k1)
        real N0_density,        // Dipole DENSITY in atoms/m³ (NOT atoms per cell!)
        const GridSpacing& grid, // Grid spacing for cell volume calculation
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
        real total_atoms = 0.0;
        real min_cell_vol = 1e30, max_cell_vol = 0.0;

        for (size_t i = i0; i < i1; ++i) {       // Half-open: [i0, i1)
            for (size_t j = j0; j < j1; ++j) {   // Half-open: [j0, j1)
                for (size_t k = k0; k < k1; ++k) { // Half-open: [k0, k1)
                    size_t id = idx3(i, j, k, NyT, NzT);

                    // CRITICAL FIX: Convert density to atoms per cell
                    // N0_density is in atoms/m³, cell_volume is in m³
                    // N0_cell = N0_density × dV has units of atoms (dimensionless count)
                    real cell_volume = grid.dx[i] * grid.dy[j] * grid.dz[k];
                    real N0_cell = N0_density * cell_volume;

                    // Set dipole count per cell (not density!)
                    Ndip[id] = N0_cell;

                    // Initialize populations with partial inversion
                    Nu[id] = inversion_fraction * N0_cell;
                    Ng[id] = (1.0 - inversion_fraction) * N0_cell;
                    Ng0[id] = N0_cell;  // Store initial total per cell

                    // Track statistics
                    total_atoms += N0_cell;
                    min_cell_vol = std::min(min_cell_vol, cell_volume);
                    max_cell_vol = std::max(max_cell_vol, cell_volume);
                }
            }
        }

        // Print initialization summary
        size_t n_cells = (i1 - i0) * (j1 - j0) * (k1 - k0);
        std::cout << "\n[TLS Init] Gain region initialized:\n";
        std::cout << "  N0_density = " << N0_density << " atoms/m³\n";
        std::cout << "  Cell volume range: " << min_cell_vol * 1e27 << " - "
                  << max_cell_vol * 1e27 << " nm³\n";
        std::cout << "  Atoms per cell range: " << N0_density * min_cell_vol
                  << " - " << N0_density * max_cell_vol << "\n";
        std::cout << "  Total atoms in gain region: " << total_atoms << "\n";
        std::cout << "  Number of cells: " << n_cells << "\n";
        std::cout << "  Initial inversion fraction: " << inversion_fraction << "\n";
        std::cout << "  Total initial inversion (Nu-Ng): "
                  << total_atoms * (2.0 * inversion_fraction - 1.0) << "\n\n";
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

                // Skip if no dipoles here
                if (state.Ndip[id] < 1e-20) continue;  // Skip only truly empty cells

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
// FIX #3: Spontaneous decay is now LINEAR in Nu (was incorrectly Nu² due to sat_term)
// FIX #4: Stimulated term now uses energy-conserving form: stim = E_avg·ΔP/(ℏωa)
//
// Correct two-level equations:
//   dNu/dt = -Nu/τ - stim     (spontaneous decay + stimulated emission/absorption)
//   dNg/dt = +Nu/τ + stim     (gains from spontaneous decay + stimulated)
//
// where stim = (E_avg · ΔP) / (ℏ · ωa)
//   - E_avg = 0.5 * (Ez_new + Ez_old)  [time-centered field]
//   - ΔP = Pz_new - Pz_old             [polarization change in this step]
//
// This form ensures:
//   1. Energy conservation: field energy change = -ℏωa × population change
//   2. Nu + Ng = Ntotal is conserved (no pumping/external losses)
//   3. Correct sign: positive inversion (Nu > Ng) causes amplification

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

                // Skip if no dipoles here
                if (state.Ndip[id] < 1e-20) continue;  // Skip only truly empty cells

                real Nu_curr = state.Nu[id];
                real Ntotal = state.Ng0[id];  // Total density (conserved)

                if (Ntotal < 1e-20) continue;  // Skip only truly empty cells

                // FIX #4: Stimulated term using energy-conserving RATE form
                // stim_rate = E_avg · (dP/dt) / (ℏωa)  [units: 1/(m³·s)]
                // E_avg = 0.5 * (Ez_new + Ez_old)
                // dP/dt = ΔP/dt where ΔP = Pz_new - Pz_old
                //
                // CRITICAL FIX (Issue A): The stimulated term must be in RATE form
                // Previous code used: stim = E_avg * ΔP / (ℏω) which has units 1/m³ (not rate!)
                // Then multiplied by dt again, causing dimension error (factor of dt² instead of dt)
                // Correct: stim_rate = E_avg * (ΔP/dt) / (ℏω) has units 1/(m³·s)
                real E_avg = 0.5 * (Ez[id] + state.Ez_old[id]);
                real delta_P = state.Pz[id] - state.Pz_prev[id];
                real dP_dt = delta_P * inv_dt;  // Convert to rate [C/(m²·s)] - uses pre-computed 1/dt
                real stim_rate = E_avg * dP_dt * inv_hbar_omega;  // [1/(m³·s)]

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

// Calculate total population in upper state (Nu)
inline real compute_total_Nu(const TwoLevelState& state) {
    real total_Nu = 0.0;
    for (size_t i = 0; i < state.Nu.size(); ++i) {
        if (state.Ndip[i] > 1e-20) {
            total_Nu += state.Nu[i];
        }
    }
    return total_Nu;
}

// Calculate total population in ground state (Ng)
inline real compute_total_Ng(const TwoLevelState& state) {
    real total_Ng = 0.0;
    for (size_t i = 0; i < state.Ng.size(); ++i) {
        if (state.Ndip[i] > 1e-20) {
            total_Ng += state.Ng[i];
        }
    }
    return total_Ng;
}

// Calculate total population inversion in gain region
inline real compute_total_inversion(const TwoLevelState& state) {
    real total_inversion = 0.0;
    for (size_t i = 0; i < state.Nu.size(); ++i) {
        // Skip truly empty cells (threshold should be much smaller than atoms per cell)
        // Typical atoms per cell: 0.1 - 1.0, so use 1e-20 as threshold
        if (state.Ndip[i] > 1e-20) {
            total_inversion += (state.Nu[i] - state.Ng[i]);
        }
    }
    return total_inversion;
}

// Calculate total stored energy in polarization
// CRITICAL FIX (Issue E): Energy formula was WRONG
// Previous code: energy_density = 0.5 * Pz * Pz (comment said "0.5 * Pz * Ez")
// Correct: energy_density = 0.5 * Pz * Ez (interaction energy between field and medium)
// Note: This is the dipole interaction energy W = -P·E, with factor 0.5 for energy storage
inline real compute_polarization_energy(
    const TwoLevelState& state,
    const std::vector<real>& Ez,  // Added Ez parameter for correct calculation
    const GridSpacing& grid
) {
    real total_energy = 0.0;
    size_t NyT = state.NyT;
    size_t NzT = state.NzT;

    for (size_t i = 0; i < state.NxT; ++i) {
        for (size_t j = 0; j < state.NyT; ++j) {
            for (size_t k = 0; k < state.NzT; ++k) {
                size_t id = idx3(i, j, k, NyT, NzT);

                if (state.Ndip[id] < 1e-20) continue;  // Skip only truly empty cells

                real dV = grid.dx[i] * grid.dy[j] * grid.dz[k];

                // Correct energy density: 0.5 * Pz * Ez (field-medium interaction)
                real energy_density = 0.5 * state.Pz[id] * Ez[id];
                total_energy += energy_density * dV;
            }
        }
    }
    return total_energy;
}