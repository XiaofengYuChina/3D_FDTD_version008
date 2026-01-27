// isource.hpp - Source interface and common definitions

#pragma once

#include <cmath>
#include <numbers>
#include <cstddef>
#include <vector>
#include <memory>
#include <string>

#include "../global_function.hpp"

namespace Sources {

enum class Waveform {
    Ricker,                 // I(t) = I0*(1-2 a^2) e^{-a^2}, a = pi f0 (t - t0)
    GaussianModulatedSine,  // I(t) = I0*sin(2*pi*f0*(t-t_shift))*exp(-((t-t0)^2)/tau^2)
    RickerLikeGaussian2nd,  // I(t) = I0*(1 - 2((t-t0)/tau)^2)*exp(-((t-t0)^2)/tau^2)
    ContinuousWave          // I(t) = I0*sin(2*pi*f0*t) with smooth turn-on
};

enum class Polarization { Ex, Ey, Ez, Hx, Hy, Hz };

struct ISource {
    virtual ~ISource() = default;

    virtual void inject_half_step(std::size_t n, real dt,
        std::vector<real>& Jx,
        std::vector<real>& Jy,
        std::vector<real>& Jz) = 0;

    virtual std::string name() const { return "ISource"; }

    virtual bool is_active(std::size_t n, real dt) const {
        (void)n; (void)dt;
        return true;
    }
};

inline real compute_waveform(Waveform type, real t, real I0, real f0,
                             real tau, real t0, real t_shift = 0.0) {
    using std::exp;
    using std::sin;
    const real pi = std::numbers::pi;

    switch (type) {
    case Waveform::Ricker: {
        const real a = pi * f0 * (t - t0);
        const real a2 = a * a;
        return I0 * (1 - 2 * a2) * exp(-a2);
    }
    case Waveform::GaussianModulatedSine: {
        const real x = (t - t0);
        return I0 * sin(2 * pi * f0 * (t - t_shift)) * exp(-(x * x) / (tau * tau));
    }
    case Waveform::RickerLikeGaussian2nd: {
        const real x = (t - t0) / tau;
        const real x2 = x * x;
        return I0 * (1 - 2 * x2) * exp(-x2);
    }
    case Waveform::ContinuousWave: {
        const real envelope = (t < tau) ? (0.5 * (1 - std::cos(pi * t / tau))) : 1.0;
        return I0 * envelope * sin(2 * pi * f0 * t);
    }
    }
    return 0;
}

// tau = 2*sqrt(ln(2)) / (pi * df_FWHM) for Gaussian envelope
inline real tau_from_bandwidth(real df_fwhm) {
    return 2.0 * std::sqrt(std::log(2.0)) / (std::numbers::pi * df_fwhm);
}

struct SourceConfig {
    real amplitude = 1e-5;
    real frequency = 0.0;
    real wavelength = 500e-9;
    real tau = 0.0;
    real df_fwhm = 0.0;
    real t0_factor = 3.0;
    Waveform waveform = Waveform::Ricker;
    Polarization polarization = Polarization::Ez;

    real x = 0.0, y = 0.0, z = 0.0;

    real kx = 0.0, ky = 0.0, kz = 1.0;
    real width = 0.0, height = 0.0;

    real get_frequency(real c0 = 299792458.0) const {
        if (frequency > 0.0) return frequency;
        if (wavelength > 0.0) return c0 / wavelength;
        return c0 / 500e-9;
    }

    real get_tau(real dt) const {
        if (tau > 0.0) return tau;
        if (df_fwhm > 0.0) return tau_from_bandwidth(df_fwhm);
        return 20.0 * dt;
    }

    real get_t0(real tau_eff, real f) const {
        if (t0_factor > 0.0) return t0_factor * tau_eff;
        return 3.0 / f;
    }
};

} // namespace Sources
