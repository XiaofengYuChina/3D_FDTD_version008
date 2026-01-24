// stability_monitor.hpp — Numerical stability monitoring for FDTD simulation
//
// Monitors simulation health and detects numerical instabilities early:
// - NaN/Inf detection in fields and populations
// - Magnitude threshold checks
// - Energy conservation monitoring
// - Growth rate tracking

#pragma once

#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>
#include "global_function.hpp"
#include "two_level_system.hpp"
#include "user_config.hpp"

struct StabilityMonitor {
    // Thresholds (initialized from UserConfig)
    real max_E_threshold;
    real max_P_threshold;
    real conservation_threshold;
    real growth_rate_threshold;

    // State tracking
    real prev_max_E = 0.0;
    size_t prev_check_step = 0;
    size_t instability_count = 0;
    bool simulation_stable = true;

    // Constructor: Initialize thresholds from user config
    StabilityMonitor()
        : max_E_threshold(UserConfig::STABILITY_MAX_E)
        , max_P_threshold(UserConfig::STABILITY_MAX_P)
        , conservation_threshold(UserConfig::STABILITY_CONSERVATION_ERR)
        , growth_rate_threshold(UserConfig::STABILITY_GROWTH_RATE)
    {}

    // Constructor with explicit thresholds (for testing or override)
    StabilityMonitor(real max_E, real max_P, real conserv_err, real growth_rate)
        : max_E_threshold(max_E)
        , max_P_threshold(max_P)
        , conservation_threshold(conserv_err)
        , growth_rate_threshold(growth_rate)
    {}

    // Full stability check - returns true if simulation should continue
    bool check_stability(
        size_t step,
        const std::vector<real>& Ex, const std::vector<real>& Ey, const std::vector<real>& Ez,
        const std::vector<real>& Hx, const std::vector<real>& Hy, const std::vector<real>& Hz,
        const TwoLevelState& tls_state,
        real conservation_error
    ) {
        bool step_stable = true;

        // Check 1: NaN/Inf in E-field
        if (NumericUtils::has_nan_or_inf(Ex) ||
            NumericUtils::has_nan_or_inf(Ey) ||
            NumericUtils::has_nan_or_inf(Ez)) {
            std::cerr << "[FATAL] Step " << step << ": NaN/Inf detected in E-field!\n";
            simulation_stable = false;
            return false;
        }

        // Check 2: NaN/Inf in H-field
        if (NumericUtils::has_nan_or_inf(Hx) ||
            NumericUtils::has_nan_or_inf(Hy) ||
            NumericUtils::has_nan_or_inf(Hz)) {
            std::cerr << "[FATAL] Step " << step << ": NaN/Inf detected in H-field!\n";
            simulation_stable = false;
            return false;
        }

        // Check 3: NaN/Inf in polarization (if TLS active)
        if (!tls_state.Pz.empty()) {
            if (NumericUtils::has_nan_or_inf(tls_state.Pz) ||
                NumericUtils::has_nan_or_inf(tls_state.dPz_dt)) {
                std::cerr << "[FATAL] Step " << step << ": NaN/Inf detected in polarization!\n";
                simulation_stable = false;
                return false;
            }
        }

        // Check 4: NaN/Inf in populations (if TLS active)
        if (!tls_state.Nu.empty()) {
            if (NumericUtils::has_nan_or_inf(tls_state.Nu) ||
                NumericUtils::has_nan_or_inf(tls_state.Ng)) {
                std::cerr << "[FATAL] Step " << step << ": NaN/Inf detected in populations!\n";
                simulation_stable = false;
                return false;
            }
        }

        // Check 5: E-field magnitude threshold
        real max_E = std::max({NumericUtils::max_abs(Ex),
                              NumericUtils::max_abs(Ey),
                              NumericUtils::max_abs(Ez)});
        if (max_E > max_E_threshold) {
            std::cerr << "[WARNING] Step " << step << ": E-field exceeded threshold: "
                      << max_E << " > " << max_E_threshold << "\n";
            instability_count++;
            step_stable = false;
        }

        // Check 6: Polarization magnitude threshold (if TLS active)
        if (!tls_state.Pz.empty()) {
            real max_P = NumericUtils::max_abs(tls_state.Pz);
            if (max_P > max_P_threshold) {
                std::cerr << "[WARNING] Step " << step << ": Polarization exceeded threshold: "
                          << max_P << " > " << max_P_threshold << "\n";
                instability_count++;
                step_stable = false;
            }
        }

        // Check 7: Conservation error
        if (conservation_error > conservation_threshold) {
            std::cerr << "[WARNING] Step " << step << ": Conservation error exceeded threshold: "
                      << conservation_error << " > " << conservation_threshold << "\n";
            step_stable = false;
        }

        // Check 8: Negative populations (physical constraint, if TLS active)
        if (!tls_state.Nu.empty()) {
            for (size_t i = 0; i < tls_state.Nu.size(); ++i) {
                if (tls_state.Ndip[i] > 1e10) {
                    if (tls_state.Nu[i] < -1e10 || tls_state.Ng[i] < -1e10) {
                        std::cerr << "[WARNING] Step " << step << ": Negative population detected!\n";
                        step_stable = false;
                        break;
                    }
                }
            }
        }

        // Check 9: Growth rate (every 100 steps)
        if (step > 0 && step % 100 == 0 && prev_max_E > 1e-20) {
            real growth_rate = max_E / prev_max_E;
            if (growth_rate > growth_rate_threshold) {
                std::cerr << "[WARNING] Step " << step << ": Rapid field growth detected: "
                          << growth_rate << "x in 100 steps\n";
                step_stable = false;
            }
            prev_max_E = max_E;
        } else if (step == 0 || prev_max_E < 1e-20) {
            prev_max_E = std::max(max_E, 1e-20);
        }

        // Terminate if too many instabilities
        if (instability_count > 10) {
            std::cerr << "[FATAL] Too many instability warnings. Terminating simulation.\n";
            simulation_stable = false;
            return false;
        }

        return simulation_stable;
    }

    // Simplified check (without TLS)
    bool check_stability_em_only(
        size_t step,
        const std::vector<real>& Ex, const std::vector<real>& Ey, const std::vector<real>& Ez,
        const std::vector<real>& Hx, const std::vector<real>& Hy, const std::vector<real>& Hz
    ) {
        // Create empty TLS state for the full check
        TwoLevelState empty_tls;
        return check_stability(step, Ex, Ey, Ez, Hx, Hy, Hz, empty_tls, 0.0);
    }

    // Calculate conservation error for TLS populations
    real calculate_conservation_error(const TwoLevelState& tls_state) const {
        if (tls_state.Nu.empty()) return 0.0;

        real total_population = 0.0;
        real initial_population = 0.0;

        for (size_t idx = 0; idx < tls_state.Nu.size(); ++idx) {
            if (tls_state.Ndip[idx] > 1e10) {
                total_population += (tls_state.Nu[idx] + tls_state.Ng[idx]);
                initial_population += tls_state.Ng0[idx];
            }
        }

        if (initial_population < 1e-30) return 0.0;
        return std::abs(total_population - initial_population) / initial_population;
    }

    // Get current max E-field for progress reporting
    real get_max_E(const std::vector<real>& Ex,
                   const std::vector<real>& Ey,
                   const std::vector<real>& Ez) const {
        return std::max({NumericUtils::max_abs(Ex),
                        NumericUtils::max_abs(Ey),
                        NumericUtils::max_abs(Ez)});
    }

    // Print stability summary
    void print_summary() const {
        std::cout << "\n=== Stability Summary ===\n";
        std::cout << "Simulation " << (simulation_stable ? "STABLE" : "UNSTABLE") << "\n";
        std::cout << "Instability warnings: " << instability_count << "\n";
        std::cout << "=========================\n";
    }

    // Print configuration
    void print_config() const {
        std::cout << "Stability monitor initialized\n";
        std::cout << "  - E-field threshold: " << max_E_threshold << " V/m\n";
        std::cout << "  - Polarization threshold: " << max_P_threshold << " C/m^2\n";
        std::cout << "  - Conservation threshold: " << conservation_threshold * 100 << "%\n";
        std::cout << "  - Growth rate threshold: " << growth_rate_threshold << "x per 100 steps\n\n";
    }

    // Reset monitor state (for re-running simulations)
    void reset() {
        prev_max_E = 0.0;
        prev_check_step = 0;
        instability_count = 0;
        simulation_stable = true;
    }
};
