// isource.hpp - Source interface and common definitions
//
// This file defines the base interface for all electromagnetic sources
// in the FDTD simulation.

#pragma once

#include <cmath>
#include <numbers>
#include <cstddef>
#include <vector>
#include <memory>
#include <string>

#include "../global_function.hpp"

namespace Sources {

// ==================== Waveform types ====================
enum class Waveform {
    Ricker,                 // I(t) = I0*(1-2 a^2) e^{-a^2}, a = pi f0 (t - t0)
    GaussianModulatedSine,  // I(t) = I0*sin(2*pi*f0*(t-t_shift))*exp(-((t-t0)^2)/tau^2)
    RickerLikeGaussian2nd,  // I(t) = I0*(1 - 2((t-t0)/tau)^2)*exp(-((t-t0)^2)/tau^2)
    ContinuousWave          // I(t) = I0*sin(2*pi*f0*t) with smooth turn-on
};

// ==================== Source polarization ====================
enum class Polarization {
    Ex,  // Electric field in x-direction
    Ey,  // Electric field in y-direction
    Ez,  // Electric field in z-direction
    Hx,  // Magnetic field in x-direction
    Hy,  // Magnetic field in y-direction
    Hz   // Magnetic field in z-direction
};

// ==================== Source interface ====================
struct ISource {
    virtual ~ISource() = default;

    // Inject current at half time step (for leapfrog integration)
    // This method is called every time step to add source contributions
    virtual void inject_half_step(std::size_t n, real dt,
        std::vector<real>& Jx,
        std::vector<real>& Jy,
        std::vector<real>& Jz) = 0;

    // Get source name for logging/debugging
    virtual std::string name() const { return "ISource"; }

    // Check if source is active at given time step
    virtual bool is_active(std::size_t n, real dt) const {
        (void)n; (void)dt;
        return true;
    }
};

// ==================== Waveform computation utilities ====================

// Compute waveform value at time t
inline real compute_waveform(Waveform type, real t, real I0, real f0,
                             real tau, real t0, real t_shift = 0.0) {
    using std::exp;
    using std::sin;
    const real pi = std::numbers::pi;

    switch (type) {
    case Waveform::Ricker: {
        // Ricker wavelet (Mexican hat), f0 is the center frequency
        const real a = pi * f0 * (t - t0);
        const real a2 = a * a;
        return I0 * (1 - 2 * a2) * exp(-a2);
    }
    case Waveform::GaussianModulatedSine: {
        // Gaussian modulated sine, tau controls bandwidth
        const real x = (t - t0);
        return I0 * sin(2 * pi * f0 * (t - t_shift)) * exp(-(x * x) / (tau * tau));
    }
    case Waveform::RickerLikeGaussian2nd: {
        // Ricker-like Gaussian second derivative
        const real x = (t - t0) / tau;
        const real x2 = x * x;
        return I0 * (1 - 2 * x2) * exp(-x2);
    }
    case Waveform::ContinuousWave: {
        // Continuous wave with smooth turn-on (over tau seconds)
        const real envelope = (t < tau) ? (0.5 * (1 - std::cos(pi * t / tau))) : 1.0;
        return I0 * envelope * sin(2 * pi * f0 * t);
    }
    }
    return 0;
}

// Compute tau from bandwidth (FWHM in frequency domain)
inline real tau_from_bandwidth(real df_fwhm) {
    // For Gaussian envelope: tau = 2*sqrt(ln(2)) / (pi * df_FWHM)
    return 2.0 * std::sqrt(std::log(2.0)) / (std::numbers::pi * df_fwhm);
}

// ==================== Source configuration structure ====================
struct SourceConfig {
    // Common parameters
    real amplitude = 1e-5;      // Peak amplitude (A for current, V/m for E-field)
    real frequency = 0.0;       // Frequency (Hz), 0 means use wavelength
    real wavelength = 500e-9;   // Wavelength (m), used if frequency <= 0
    real tau = 0.0;             // Time constant (s), 0 means auto-compute
    real df_fwhm = 0.0;         // Bandwidth FWHM (Hz), used if tau <= 0
    real t0_factor = 3.0;       // t0 = t0_factor * tau_eff
    Waveform waveform = Waveform::Ricker;
    Polarization polarization = Polarization::Ez;

    // Position in physical coordinates (meters)
    real x = 0.0;
    real y = 0.0;
    real z = 0.0;

    // For plane wave: propagation direction and size
    real kx = 0.0, ky = 0.0, kz = 1.0;  // Propagation direction (normalized)
    real width = 0.0, height = 0.0;      // Plane wave extent (0 = infinite)

    // Compute effective frequency
    real get_frequency(real c0 = 299792458.0) const {
        if (frequency > 0.0) return frequency;
        if (wavelength > 0.0) return c0 / wavelength;
        return c0 / 500e-9;  // Default 500nm
    }

    // Compute effective tau
    real get_tau(real dt) const {
        if (tau > 0.0) return tau;
        if (df_fwhm > 0.0) return tau_from_bandwidth(df_fwhm);
        return 20.0 * dt;  // Default
    }

    // Compute t0
    real get_t0(real tau_eff, real f) const {
        if (t0_factor > 0.0) return t0_factor * tau_eff;
        return 3.0 / f;  // Default: 3 periods
    }
};

} // namespace Sources
