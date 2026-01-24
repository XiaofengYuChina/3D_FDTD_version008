// dipole_source.hpp - Point current (dipole) source implementation
//
// A dipole source is a point current source that injects current at a single
// grid cell. It's the simplest and most common source type for FDTD simulations.
//
// Usage:
//   auto src = Sources::make_dipole_source(config, grid_spacing, NyT, NzT, dt);
//   src->inject_half_step(n, dt, Jx, Jy, Jz);

#pragma once

#include "isource.hpp"
#include <iostream>
#include <sstream>

namespace Sources {

// ==================== Dipole Source ====================
// A point current source that injects current at a single grid cell
// Supports x, y, or z polarization
struct DipoleSource final : public ISource {
    // Grid indices and dimensions
    std::size_t i{}, j{}, k{};
    std::size_t NyT{}, NzT{};

    // Cell area for current density calculation (J = I / area)
    real cell_area{};

    // Source parameters
    real I0{};            // Peak current (A)
    real f0{};            // Frequency (Hz)
    real tau{};           // Gaussian time constant (s)
    real t0{};            // Source center time (s)
    real t_shift{};       // Phase shift offset (for Gaussian modulated sine)
    Waveform waveform{Waveform::Ricker};
    Polarization polarization{Polarization::Ez};

    // Source name for logging
    std::string source_name{"DipoleSource"};

    // Default constructor
    DipoleSource() = default;

    // Full constructor
    DipoleSource(std::size_t ii, std::size_t jj, std::size_t kk,
                 std::size_t NyTot, std::size_t NzTot,
                 real cell_area_,
                 real I0_, real f0_, real tau_, real t0_,
                 Waveform wf = Waveform::Ricker,
                 Polarization pol = Polarization::Ez,
                 real t_shift_ = 0.0)
        : i(ii), j(jj), k(kk), NyT(NyTot), NzT(NzTot),
          cell_area(cell_area_),
          I0(I0_), f0(f0_), tau(tau_), t0(t0_),
          t_shift(t_shift_), waveform(wf), polarization(pol)
    {
        std::ostringstream oss;
        oss << "DipoleSource[" << polarization_name() << "]@(" << i << "," << j << "," << k << ")";
        source_name = oss.str();
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
        return compute_waveform(waveform, t, I0, f0, tau, t0, t_shift);
    }

    // ISource interface implementation
    void inject_half_step(std::size_t n, real dt,
                          std::vector<real>& Jx,
                          std::vector<real>& Jy,
                          std::vector<real>& Jz) override
    {
        const std::size_t id = idx3(i, j, k, NyT, NzT);
        const real t_half = (static_cast<real>(n) + real(0.5)) * dt;
        const real It = waveform_value(t_half);  // Current (A)
        const real Jval = It / cell_area;        // Current density (A/m^2)

        switch (polarization) {
        case Polarization::Ex:
        case Polarization::Hx:
            Jx[id] += Jval;
            break;
        case Polarization::Ey:
        case Polarization::Hy:
            Jy[id] += Jval;
            break;
        case Polarization::Ez:
        case Polarization::Hz:
        default:
            Jz[id] += Jval;
            break;
        }
    }

    std::string name() const override { return source_name; }
};

// ==================== Factory function ====================
// Create a dipole source from configuration and grid information
//
// Parameters:
//   config       - Source configuration (position, amplitude, frequency, etc.)
//   grid_spacing - Grid spacing information for coordinate conversion
//   NyT, NzT     - Total grid dimensions (including PML)
//   dt           - Time step (for tau calculation if needed)
//   domain_min_x/y/z - Domain minimum coordinates (for physical->grid conversion)
//
inline std::unique_ptr<ISource> make_dipole_source(
    const SourceConfig& config,
    const GridSpacing& grid_spacing,
    std::size_t NyT, std::size_t NzT,
    real dt,
    real domain_min_x = 0.0,
    real domain_min_y = 0.0,
    real domain_min_z = 0.0,
    real c0 = 299792458.0)
{
    // Convert physical coordinates to grid indices
    real phys_x = config.x - domain_min_x;
    real phys_y = config.y - domain_min_y;
    real phys_z = config.z - domain_min_z;

    std::size_t i = grid_spacing.physical_to_index_x(phys_x);
    std::size_t j = grid_spacing.physical_to_index_y(phys_y);
    std::size_t k = grid_spacing.physical_to_index_z(phys_z);

    // Get frequency and timing parameters
    real f = config.get_frequency(c0);
    real tau_eff = config.get_tau(dt);
    real t0 = config.get_t0(tau_eff, f);
    real t_shift = t0;  // For Gaussian modulated sine

    // Compute cell area based on polarization
    real cell_area = 1.0;
    switch (config.polarization) {
    case Polarization::Ex:
    case Polarization::Hx:
        cell_area = grid_spacing.dy[j] * grid_spacing.dz[k];
        break;
    case Polarization::Ey:
    case Polarization::Hy:
        cell_area = grid_spacing.dx[i] * grid_spacing.dz[k];
        break;
    case Polarization::Ez:
    case Polarization::Hz:
    default:
        cell_area = grid_spacing.dx[i] * grid_spacing.dy[j];
        break;
    }

    auto src = std::make_unique<DipoleSource>(
        i, j, k, NyT, NzT, cell_area,
        config.amplitude, f, tau_eff, t0,
        config.waveform, config.polarization, t_shift
    );

    // Log source creation
    real actual_x = grid_spacing.cell_center_physical_x(i);
    real actual_y = grid_spacing.cell_center_physical_y(j);
    real actual_z = grid_spacing.cell_center_physical_z(k);

    std::cout << "[DipoleSource] Created:\n";
    std::cout << "  Polarization: " << src->polarization_name() << "\n";
    std::cout << "  Target position: (" << phys_x * 1e9 << ", "
              << phys_y * 1e9 << ", " << phys_z * 1e9 << ") nm\n";
    std::cout << "  Grid indices: (" << i << ", " << j << ", " << k << ")\n";
    std::cout << "  Actual position: (" << actual_x * 1e9 << ", "
              << actual_y * 1e9 << ", " << actual_z * 1e9 << ") nm\n";
    std::cout << "  Amplitude: " << config.amplitude << " A\n";
    std::cout << "  Frequency: " << f * 1e-12 << " THz (lambda = " << c0 / f * 1e9 << " nm)\n";
    std::cout << "  Waveform: " << static_cast<int>(config.waveform) << "\n";
    std::cout << "  t0: " << t0 * 1e15 << " fs, tau: " << tau_eff * 1e15 << " fs\n";

    return src;
}

// Convenience function for backward compatibility with old interface
inline std::unique_ptr<ISource> make_point_current(
    std::size_t i, std::size_t j, std::size_t k,
    std::size_t NyT, std::size_t NzT,
    const GridSpacing& grid_spacing,
    real I0,
    real f0,
    real lambda0,
    real tau_src,
    real df_fwhm,
    real dt,
    real t0_factor,
    Waveform type,
    Polarization pol = Polarization::Ez,
    real c0 = 299792458.0)
{
    // Compute frequency
    real f = f0;
    if (f <= 0.0) {
        if (lambda0 <= 0.0) return nullptr;
        f = c0 / lambda0;
    }

    // Compute tau
    real tau = 0.0;
    if (tau_src > 0.0) {
        tau = tau_src;
    } else if (df_fwhm > 0.0) {
        tau = tau_from_bandwidth(df_fwhm);
    } else {
        tau = 20.0 * dt;
    }

    // Compute t0
    real t0 = (t0_factor > 0.0) ? (t0_factor * tau) : (3.0 / f);
    real t_shift = t0;

    // Compute cell area based on polarization
    real cell_area = 1.0;
    switch (pol) {
    case Polarization::Ex:
    case Polarization::Hx:
        cell_area = grid_spacing.dy[j] * grid_spacing.dz[k];
        break;
    case Polarization::Ey:
    case Polarization::Hy:
        cell_area = grid_spacing.dx[i] * grid_spacing.dz[k];
        break;
    case Polarization::Ez:
    case Polarization::Hz:
    default:
        cell_area = grid_spacing.dx[i] * grid_spacing.dy[j];
        break;
    }

    return std::make_unique<DipoleSource>(
        i, j, k, NyT, NzT, cell_area,
        I0, f, tau, t0, type, pol, t_shift
    );
}

} // namespace Sources
