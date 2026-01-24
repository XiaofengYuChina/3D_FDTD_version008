// source.hpp —— Source

#pragma once
#include <cmath>
#include <numbers>
#include <cstddef>
#include <vector>

#include "global_function.hpp"

// -------------------- Waveform types --------------------
enum class Waveform {
    Ricker,                 // I(t) = I0*(1-2 a^2) e^{-a^2}, a = pi f0 (t - t0)
    GaussianModulatedSine,  // I(t) = I0*sin(2πf0*(t-t_shift))*exp(-((t-t0)^2)/tau^2)
    RickerLikeGaussian2nd   // I(t) = I0*(1 - 2((t-t0)/tau)^2)*exp(-((t-t0)^2)/tau^2)
};

// -------------------- Source interface --------------------
struct ISource {
    virtual ~ISource() = default;
    virtual void inject_half_step(std::size_t n, real dt,
        std::vector<real>& Jx,
        std::vector<real>& Jy,
        std::vector<real>& Jz) = 0;
};

// -------------------- Point current source (Z direction) --------------------
struct PointCurrentZ final : public ISource {
    // Position indices and grid dimensions
    std::size_t i{}, j{}, k{}, NyT{}, NzT{};

    // Grid dimensions (for calculating J = I/area)
    real cell_area{};  // CHANGED: Store area instead of dx, dy separately

    // Source parameters
    real I0{};            // Peak current (A)
    real f0{};            // Frequency (Hz)
    real tau{};           // Gaussian time constant (s) (for some waveforms, can be ignored)
    real t0{};            // Source center time (s)
    real t_shift{};       // Phase shift offset (for Gaussian modulated sine)
    Waveform type{ Waveform::Ricker };

    // Constructor: only need to pass required parameters, can be called from factory/external code, no need to repeat logic
    PointCurrentZ(std::size_t ii, std::size_t jj, std::size_t kk,
        std::size_t NyTot, std::size_t NzTot,
        real cell_area_,  // CHANGED
        real I0_, real f0_, real tau_, real t0_,
        Waveform wf = Waveform::Ricker,
        real t_shift_ = 0.0)
        : i(ii), j(jj), k(kk), NyT(NyTot), NzT(NzTot),
        cell_area(cell_area_),
        I0(I0_), f0(f0_), tau(tau_), t0(t0_),
        t_shift(t_shift_), type(wf) {}

    inline real waveform_value(real t) const {
        using std::exp;
        using std::sin;
        const real pi = std::numbers::pi;
        switch (type) {
        case Waveform::Ricker: {
            // Ricker wavelet (Mexican hat), f0 is the same as the main frequency
            const real a = pi * f0 * (t - t0);
            const real a2 = a * a;
            return I0 * (1 - 2 * a2) * exp(-a2);
        }
        case Waveform::GaussianModulatedSine: {
            // Gaussian modulated sine, can control df_FWHM through tau
            // t_shift adjusts phase to make symmetric, usually take t_shift = t0
            const real x = (t - t0);
            return I0 * sin(2 * pi * f0 * (t - t_shift)) * exp(-(x * x) / (tau * tau));
        }
        case Waveform::RickerLikeGaussian2nd: {
            // Ricker-like Gaussian second derivative: I ~ (1 - 2 (x/tau)^2) e^{-(x/tau)^2}
            const real x = (t - t0) / tau;
            const real x2 = x * x;
            return I0 * (1 - 2 * x2) * exp(-x2);
        }
        }
        return 0;
    }

    void inject_half_step(std::size_t n, real dt,
        std::vector<real>& Jx,
        std::vector<real>& Jy,
        std::vector<real>& Jz) override
    {
        (void)Jx; (void)Jy;
        const std::size_t id = idx3(i, j, k, NyT, NzT);
        const real t_half = (static_cast<real>(n) + real(0.5)) * dt;
        const real It = waveform_value(t_half); // A
        Jz[id] = It / cell_area;                // A/m^2
    }
};

// -------------------- Factory function: create point current source (Z direction) --------------------
// Notes:
// - f0 priority, if f0<=0 then use lambda0 (wavelength), f0 = c0 / lambda0, if both are invalid return nullptr
// - tau selection: tau_src>0 priority, then if df_fwhm>0 use formula calculation; otherwise default tau = 20*dt
// - t0 = (t0_factor > 0 ? t0_factor * tau : 3/f0)
// - For Gaussian modulated sine, t_shift = t0 makes phase centered; Ricker doesn't use tau but can ignore
//
inline std::unique_ptr<ISource>
make_point_current_z(std::size_t i, std::size_t j, std::size_t k,
    std::size_t NyT, std::size_t NzT,
    const GridSpacing& grid_spacing,  // CHANGED
    real I0,
    real f0,            // If f0<=0 then use lambda0
    real lambda0,       // m, if f0>0 can ignore
    real tau_src,       // s, effective if >0
    real df_fwhm,       // Hz, used when tau_src<=0
    real dt,            // s, for default tau calculation
    real t0_factor,     // >0 means t0 = t0_factor * tau
    Waveform type,
    real c0 = 299792458.0)
{
    // Frequency
    real f = f0;
    if (f <= 0.0) {
        if (lambda0 <= 0.0) return nullptr;     // Both invalid, cannot create
        f = c0 / lambda0;
    }

    // Time constant
    real tau = 0.0;
    if (tau_src > 0.0) {
        tau = tau_src;
    }
    else if (df_fwhm > 0.0) {
        tau = tau_from_df_fwhm(df_fwhm);
    }
    else {
        tau = 20.0 * dt; // Default value, sufficient to start from 0
    }

    // Source center time
    real t0 = (t0_factor > 0.0) ? (t0_factor * tau) : (3.0 / f);

    // For Gaussian modulated sine, phase adjustment: t_shift = t0, other waveforms can ignore this parameter
    const real t_shift = t0;

    // CHANGED: Compute cell area from local spacing
    real cell_area = grid_spacing.dy[j] * grid_spacing.dz[k];

    return std::make_unique<PointCurrentZ>(
        i, j, k, NyT, NzT, cell_area,
        I0, f, tau, t0, type, t_shift
    );
}
