// dipole_source.hpp - Point current (dipole) source implementation

#pragma once

#include "isource.hpp"
#include <iostream>
#include <sstream>

namespace Sources {

struct DipoleSource final : public ISource {
    std::size_t i{}, j{}, k{};
    std::size_t NyT{}, NzT{};
    real cell_area{};

    real I0{};
    real f0{};
    real tau{};
    real t0{};
    real t_shift{};
    Waveform waveform{Waveform::Ricker};
    Polarization polarization{Polarization::Ez};
    std::string source_name{"DipoleSource"};

    DipoleSource() = default;

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

    inline real waveform_value(real t) const {
        return compute_waveform(waveform, t, I0, f0, tau, t0, t_shift);
    }

    void inject_half_step(std::size_t n, real dt,
                          std::vector<real>& Jx,
                          std::vector<real>& Jy,
                          std::vector<real>& Jz) override
    {
        const std::size_t id = idx3(i, j, k, NyT, NzT);
        const real t_half = (static_cast<real>(n) + real(0.5)) * dt;
        const real It = waveform_value(t_half);
        const real Jval = It / cell_area;

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
    real phys_x = config.x - domain_min_x;
    real phys_y = config.y - domain_min_y;
    real phys_z = config.z - domain_min_z;

    std::size_t i = grid_spacing.physical_to_index_x(phys_x);
    std::size_t j = grid_spacing.physical_to_index_y(phys_y);
    std::size_t k = grid_spacing.physical_to_index_z(phys_z);

    real f = config.get_frequency(c0);
    real tau_eff = config.get_tau(dt);
    real t0 = config.get_t0(tau_eff, f);
    real t_shift = t0;

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
    real f = f0;
    if (f <= 0.0) {
        if (lambda0 <= 0.0) return nullptr;
        f = c0 / lambda0;
    }

    real tau = 0.0;
    if (tau_src > 0.0) {
        tau = tau_src;
    } else if (df_fwhm > 0.0) {
        tau = tau_from_bandwidth(df_fwhm);
    } else {
        tau = 20.0 * dt;
    }

    real t0 = (t0_factor > 0.0) ? (t0_factor * tau) : (3.0 / f);
    real t_shift = t0;

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
