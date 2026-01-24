// plane_wave_source.hpp - Plane wave source implementation
//
// A plane wave source injects a plane wave into the simulation domain.
// This implementation uses the soft source technique where current is
// injected along a plane perpendicular to the propagation direction.
//
// For a plane wave propagating in +z direction with E-field in x direction:
//   E_x(z, t) = E0 * f(t - z/c)
//   H_y(z, t) = E0/eta0 * f(t - z/c)
//
// Usage:
//   auto src = Sources::make_plane_wave_source(config, grid_spacing, ...);
//   src->inject_half_step(n, dt, Jx, Jy, Jz);

#pragma once

#include "isource.hpp"
#include <iostream>
#include <cmath>

namespace Sources {

// ==================== Plane Wave Direction ====================
enum class PlaneWaveDirection {
    PlusX,   // Propagating in +X direction
    MinusX,  // Propagating in -X direction
    PlusY,   // Propagating in +Y direction
    MinusY,  // Propagating in -Y direction
    PlusZ,   // Propagating in +Z direction (default)
    MinusZ   // Propagating in -Z direction
};

// ==================== Plane Wave Source ====================
// Injects a plane wave using soft source technique on a 2D plane
struct PlaneWaveSource final : public ISource {
    // Grid dimensions
    std::size_t NxT{}, NyT{}, NzT{};

    // Injection plane position (index)
    std::size_t plane_index{};
    PlaneWaveDirection direction{PlaneWaveDirection::PlusZ};
    Polarization polarization{Polarization::Ex};

    // Injection region bounds (in grid indices)
    std::size_t i_min{}, i_max{};
    std::size_t j_min{}, j_max{};
    std::size_t k_min{}, k_max{};

    // Source parameters
    real E0{};            // Peak E-field amplitude (V/m)
    real f0{};            // Frequency (Hz)
    real tau{};           // Gaussian time constant (s)
    real t0{};            // Source center time (s)
    real t_shift{};       // Phase shift offset
    Waveform waveform{Waveform::GaussianModulatedSine};

    // Grid spacing arrays (for current density calculation)
    std::vector<real> dx_array, dy_array, dz_array;

    // Physical constants
    static constexpr real c0 = 299792458.0;
    static constexpr real eps0 = 8.854187817e-12;
    static constexpr real mu0 = 1.2566370614359173e-6;
    static constexpr real eta0 = 376.730313668;  // sqrt(mu0/eps0)

    // Source name
    std::string source_name{"PlaneWaveSource"};

    // Default constructor
    PlaneWaveSource() = default;

    // Get direction name
    const char* direction_name() const {
        switch (direction) {
        case PlaneWaveDirection::PlusX:  return "+X";
        case PlaneWaveDirection::MinusX: return "-X";
        case PlaneWaveDirection::PlusY:  return "+Y";
        case PlaneWaveDirection::MinusY: return "-Y";
        case PlaneWaveDirection::PlusZ:  return "+Z";
        case PlaneWaveDirection::MinusZ: return "-Z";
        }
        return "Unknown";
    }

    // Get polarization name
    const char* polarization_name() const {
        switch (polarization) {
        case Polarization::Ex: return "Ex";
        case Polarization::Ey: return "Ey";
        case Polarization::Ez: return "Ez";
        case Polarization::Hx: return "Hx";
        case Polarization::Hy: return "Hy";
        case Polarization::Hz: return "Hz";
        }
        return "Unknown";
    }

    // Compute waveform value at time t
    inline real waveform_value(real t) const {
        return compute_waveform(waveform, t, E0, f0, tau, t0, t_shift);
    }

    // ISource interface implementation
    void inject_half_step(std::size_t n, real dt,
                          std::vector<real>& Jx,
                          std::vector<real>& Jy,
                          std::vector<real>& Jz) override
    {
        const real t_half = (static_cast<real>(n) + real(0.5)) * dt;

        // For a plane wave, we inject current density on the plane
        // J = (2/eta0) * E for soft source
        // The factor depends on the implementation (hard vs soft source)

        switch (direction) {
        case PlaneWaveDirection::PlusZ:
        case PlaneWaveDirection::MinusZ:
            inject_z_plane(t_half, Jx, Jy);
            break;
        case PlaneWaveDirection::PlusX:
        case PlaneWaveDirection::MinusX:
            inject_x_plane(t_half, Jy, Jz);
            break;
        case PlaneWaveDirection::PlusY:
        case PlaneWaveDirection::MinusY:
            inject_y_plane(t_half, Jx, Jz);
            break;
        }
    }

    std::string name() const override { return source_name; }

private:
    // Inject on XY plane (propagating in Z direction)
    void inject_z_plane(real t, std::vector<real>& Jx, std::vector<real>& Jy) {
        const real Et = waveform_value(t);
        const real sign = (direction == PlaneWaveDirection::PlusZ) ? 1.0 : -1.0;

        for (std::size_t i = i_min; i <= i_max; ++i) {
            for (std::size_t j = j_min; j <= j_max; ++j) {
                const std::size_t id = idx3(i, j, plane_index, NyT, NzT);
                const real dz_local = dz_array[plane_index];

                // Current density for soft source: J = 2*E/(eta0*dz)
                // Factor of 2 because it's a soft source (adds to both sides)
                const real J_amplitude = 2.0 * Et / (eta0 * dz_local) * sign;

                switch (polarization) {
                case Polarization::Ex:
                case Polarization::Hx:
                    Jx[id] += J_amplitude;
                    break;
                case Polarization::Ey:
                case Polarization::Hy:
                    Jy[id] += J_amplitude;
                    break;
                default:
                    // Ez polarization for Z-propagating wave needs special handling
                    break;
                }
            }
        }
    }

    // Inject on YZ plane (propagating in X direction)
    void inject_x_plane(real t, std::vector<real>& Jy, std::vector<real>& Jz) {
        const real Et = waveform_value(t);
        const real sign = (direction == PlaneWaveDirection::PlusX) ? 1.0 : -1.0;

        for (std::size_t j = j_min; j <= j_max; ++j) {
            for (std::size_t k = k_min; k <= k_max; ++k) {
                const std::size_t id = idx3(plane_index, j, k, NyT, NzT);
                const real dx_local = dx_array[plane_index];
                const real J_amplitude = 2.0 * Et / (eta0 * dx_local) * sign;

                switch (polarization) {
                case Polarization::Ey:
                case Polarization::Hy:
                    Jy[id] += J_amplitude;
                    break;
                case Polarization::Ez:
                case Polarization::Hz:
                    Jz[id] += J_amplitude;
                    break;
                default:
                    break;
                }
            }
        }
    }

    // Inject on XZ plane (propagating in Y direction)
    void inject_y_plane(real t, std::vector<real>& Jx, std::vector<real>& Jz) {
        const real Et = waveform_value(t);
        const real sign = (direction == PlaneWaveDirection::PlusY) ? 1.0 : -1.0;

        for (std::size_t i = i_min; i <= i_max; ++i) {
            for (std::size_t k = k_min; k <= k_max; ++k) {
                const std::size_t id = idx3(i, plane_index, k, NyT, NzT);
                const real dy_local = dy_array[plane_index];
                const real J_amplitude = 2.0 * Et / (eta0 * dy_local) * sign;

                switch (polarization) {
                case Polarization::Ex:
                case Polarization::Hx:
                    Jx[id] += J_amplitude;
                    break;
                case Polarization::Ez:
                case Polarization::Hz:
                    Jz[id] += J_amplitude;
                    break;
                default:
                    break;
                }
            }
        }
    }
};

// ==================== Plane Wave Configuration ====================
struct PlaneWaveConfig {
    // Wave parameters
    real amplitude = 1.0;         // Peak E-field (V/m)
    real frequency = 0.0;         // Frequency (Hz), 0 = use wavelength
    real wavelength = 500e-9;     // Wavelength (m)
    real tau = 0.0;               // Time constant (s), 0 = auto
    real df_fwhm = 0.0;           // Bandwidth FWHM (Hz)
    real t0_factor = 3.0;         // t0 = t0_factor * tau_eff
    Waveform waveform = Waveform::GaussianModulatedSine;

    // Propagation and polarization
    PlaneWaveDirection direction = PlaneWaveDirection::PlusZ;
    Polarization polarization = Polarization::Ex;

    // Injection plane position (in physical coordinates, meters)
    real injection_position = 0.0;

    // Injection region (in physical coordinates, meters)
    // If all zeros, uses full domain
    real x_min = 0.0, x_max = 0.0;
    real y_min = 0.0, y_max = 0.0;
    real z_min = 0.0, z_max = 0.0;

    // Compute effective frequency
    real get_frequency(real c0 = 299792458.0) const {
        if (frequency > 0.0) return frequency;
        if (wavelength > 0.0) return c0 / wavelength;
        return c0 / 500e-9;
    }

    // Compute effective tau
    real get_tau(real dt) const {
        if (tau > 0.0) return tau;
        if (df_fwhm > 0.0) return tau_from_bandwidth(df_fwhm);
        return 20.0 * dt;
    }
};

// ==================== Factory function ====================
inline std::unique_ptr<ISource> make_plane_wave_source(
    const PlaneWaveConfig& config,
    const GridSpacing& grid_spacing,
    std::size_t NxT, std::size_t NyT, std::size_t NzT,
    std::size_t npml,
    real dt,
    real domain_min_x = 0.0,
    real domain_min_y = 0.0,
    real domain_min_z = 0.0,
    real c0 = 299792458.0)
{
    auto src = std::make_unique<PlaneWaveSource>();

    src->NxT = NxT;
    src->NyT = NyT;
    src->NzT = NzT;
    src->direction = config.direction;
    src->polarization = config.polarization;
    src->waveform = config.waveform;

    // Copy grid spacing
    src->dx_array = grid_spacing.dx;
    src->dy_array = grid_spacing.dy;
    src->dz_array = grid_spacing.dz;

    // Compute frequency and timing
    src->f0 = config.get_frequency(c0);
    src->tau = config.get_tau(dt);
    src->t0 = (config.t0_factor > 0.0) ? (config.t0_factor * src->tau) : (3.0 / src->f0);
    src->t_shift = src->t0;
    src->E0 = config.amplitude;

    // Determine injection plane position
    real inj_pos = config.injection_position;
    switch (config.direction) {
    case PlaneWaveDirection::PlusZ:
    case PlaneWaveDirection::MinusZ:
        src->plane_index = grid_spacing.physical_to_index_z(inj_pos - domain_min_z);
        break;
    case PlaneWaveDirection::PlusX:
    case PlaneWaveDirection::MinusX:
        src->plane_index = grid_spacing.physical_to_index_x(inj_pos - domain_min_x);
        break;
    case PlaneWaveDirection::PlusY:
    case PlaneWaveDirection::MinusY:
        src->plane_index = grid_spacing.physical_to_index_y(inj_pos - domain_min_y);
        break;
    }

    // Clamp to valid range (inside PML)
    src->plane_index = std::max(npml + 1, std::min(src->plane_index,
        (config.direction == PlaneWaveDirection::PlusZ || config.direction == PlaneWaveDirection::MinusZ)
            ? NzT - npml - 2
            : (config.direction == PlaneWaveDirection::PlusX || config.direction == PlaneWaveDirection::MinusX)
                ? NxT - npml - 2
                : NyT - npml - 2));

    // Determine injection region
    auto compute_bounds = [&](real min_phys, real max_phys, std::size_t N, std::size_t npml_,
                              std::size_t& out_min, std::size_t& out_max,
                              auto to_index_func) {
        if (min_phys == 0.0 && max_phys == 0.0) {
            // Use full domain (excluding PML)
            out_min = npml_;
            out_max = N - npml_ - 1;
        } else {
            out_min = to_index_func(min_phys);
            out_max = to_index_func(max_phys);
            out_min = std::max(npml_, out_min);
            out_max = std::min(N - npml_ - 1, out_max);
        }
    };

    compute_bounds(config.x_min - domain_min_x, config.x_max - domain_min_x, NxT, npml,
                   src->i_min, src->i_max,
                   [&](real x) { return grid_spacing.physical_to_index_x(x); });

    compute_bounds(config.y_min - domain_min_y, config.y_max - domain_min_y, NyT, npml,
                   src->j_min, src->j_max,
                   [&](real y) { return grid_spacing.physical_to_index_y(y); });

    compute_bounds(config.z_min - domain_min_z, config.z_max - domain_min_z, NzT, npml,
                   src->k_min, src->k_max,
                   [&](real z) { return grid_spacing.physical_to_index_z(z); });

    // Build source name
    std::ostringstream oss;
    oss << "PlaneWave[" << src->direction_name() << "," << src->polarization_name()
        << "]@" << src->plane_index;
    src->source_name = oss.str();

    // Log source creation
    std::cout << "[PlaneWaveSource] Created:\n";
    std::cout << "  Direction: " << src->direction_name() << "\n";
    std::cout << "  Polarization: " << src->polarization_name() << "\n";
    std::cout << "  Injection plane index: " << src->plane_index << "\n";
    std::cout << "  Region: i=[" << src->i_min << "," << src->i_max << "], "
              << "j=[" << src->j_min << "," << src->j_max << "], "
              << "k=[" << src->k_min << "," << src->k_max << "]\n";
    std::cout << "  E0: " << src->E0 << " V/m\n";
    std::cout << "  Frequency: " << src->f0 * 1e-12 << " THz (lambda = "
              << c0 / src->f0 * 1e9 << " nm)\n";
    std::cout << "  t0: " << src->t0 * 1e15 << " fs, tau: " << src->tau * 1e15 << " fs\n";

    return src;
}

} // namespace Sources
