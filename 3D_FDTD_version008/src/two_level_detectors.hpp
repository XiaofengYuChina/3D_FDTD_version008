// two_level_detectors.hpp - Detectors for two-level system diagnostics

#pragma once

#include <filesystem>
#include <fstream>
#include <iostream>
#include <iomanip>
#include "global_function.hpp"
#include "two_level_system.hpp"

namespace fs = std::filesystem;

// Nu, Ng are densities [m^-3].
// Binary files store only the first and last frames (for spatial depletion analysis).
// CSV integrates with cell volume for total atom counts every saveEvery steps.
struct PopulationTimeSeriesDetector {
    fs::path det_dir;
    fs::path Nu_bin, Ng_bin, inversion_csv;
    std::ofstream f_Nu, f_Ng, f_inv;
    size_t NxT, NyT, NzT;
    size_t saveEvery, nSteps;
    real dt;
    const GridSpacing* grid_ptr;
    bool first_frame_saved{false};
    size_t last_saved_step{0};

    PopulationTimeSeriesDetector(
        const fs::path& out_root, std::string det_name,
        size_t NxT_, size_t NyT_, size_t NzT_,
        size_t saveEvery_, size_t nSteps_, real dt_,
        const GridSpacing& grid
    ) : det_dir(out_root / det_name),
        NxT(NxT_), NyT(NyT_), NzT(NzT_),
        saveEvery(saveEvery_), nSteps(nSteps_), dt(dt_),
        grid_ptr(&grid)
    {
        std::error_code ec;
        fs::create_directories(det_dir, ec);
        if (ec) std::cerr << "[ERR] create " << det_dir << ": " << ec.message() << "\n";

        Nu_bin = det_dir / "populations_Nu.bin";
        Ng_bin = det_dir / "populations_Ng.bin";
        inversion_csv = det_dir / "inversion_time.csv";

        f_Nu.open(Nu_bin, std::ios::binary);
        f_Ng.open(Ng_bin, std::ios::binary);
        f_inv.open(inversion_csv, std::ios::out);

        if (!f_Nu || !f_Ng || !f_inv)
            std::cerr << "[ERR] Failed to open population detector files\n";

        if (f_inv)
            f_inv << "t, total_Nu, total_Ng, inversion, inv_fraction\n";

        std::cout << "[PopulationTS] " << det_dir << "\n";
        std::cout << "  Binary: first + last frame only\n";
    }

    void write_binary_frame(const TwoLevelState& state) {
        const size_t N = state.Nu.size();
        for (size_t i = 0; i < N; ++i) {
            double nu_val = static_cast<double>(state.Nu[i]);
            double ng_val = static_cast<double>(state.Ng[i]);
            f_Nu.write(reinterpret_cast<const char*>(&nu_val), sizeof(double));
            f_Ng.write(reinterpret_cast<const char*>(&ng_val), sizeof(double));
        }
    }

    void record(size_t n, const TwoLevelState& state) {
        if (n % saveEvery != 0) return;
        if (!f_inv) return;

        // Binary: save only first frame
        if (!first_frame_saved && f_Nu && f_Ng) {
            write_binary_frame(state);
            first_frame_saved = true;
        }
        last_saved_step = n;

        // CSV: integrate density * dV to get total atom counts
        double total_Nu = 0.0, total_Ng = 0.0;
        for (size_t i = 0; i < NxT; ++i) {
            for (size_t j = 0; j < NyT; ++j) {
                for (size_t k = 0; k < NzT; ++k) {
                    size_t id = idx3(i, j, k, NyT, NzT);
                    if (state.Ndip[id] <= 0) continue;
                    real dV = grid_ptr->dx[i] * grid_ptr->dy[j] * grid_ptr->dz[k];
                    total_Nu += state.Nu[id] * dV;
                    total_Ng += state.Ng[id] * dV;
                }
            }
        }

        double inversion = total_Nu - total_Ng;
        double denom = total_Nu + total_Ng;
        double inv_fraction = (denom > 1e-30) ? (inversion / denom) : 0.0;

        real t = n * dt;
        f_inv << std::setprecision(18) << t << ","
              << total_Nu << "," << total_Ng << ","
              << inversion << "," << inv_fraction << "\n";
    }

    // Call at end of simulation to save the last frame binary snapshot
    void finalize(const TwoLevelState& state) {
        if (f_Nu && f_Ng) {
            write_binary_frame(state);
        }

        // Write metadata after we know how many frames were stored
        std::ofstream meta(det_dir / "metadata_populations.json");
        if (meta) {
            meta << "{\n"
                 << "  \"NxT\": " << NxT << ",\n"
                 << "  \"NyT\": " << NyT << ",\n"
                 << "  \"NzT\": " << NzT << ",\n"
                 << "  \"saveEvery\": " << saveEvery << ",\n"
                 << "  \"nSteps\": " << nSteps << ",\n"
                 << "  \"dt\": " << std::setprecision(18) << dt << ",\n"
                 << "  \"dtype\": \"float64\",\n"
                 << "  \"binary_frames\": 2,\n"
                 << "  \"binary_mode\": \"first_and_last\",\n"
                 << "  \"binary_data_units\": \"m^-3 (density)\",\n"
                 << "  \"csv_total_units\": \"atoms (integrated)\",\n"
                 << "  \"what\": \"population_densities_Nu_and_Ng_vs_time\"\n"
                 << "}\n";
        }

        if (f_Nu) f_Nu.close();
        if (f_Ng) f_Ng.close();
        if (f_inv) f_inv.close();
        std::cout << "[PopulationTS] Finalized: 2 binary frames + CSV time series\n";
    }

    ~PopulationTimeSeriesDetector() {
        if (f_Nu) f_Nu.close();
        if (f_Ng) f_Ng.close();
        if (f_inv) f_inv.close();
    }
};

// Spatial population snapshot (2D slice). Inversion data is density [m^-3].
struct PopulationSlice2D {
    fs::path det_dir;
    size_t NxT, NyT, NzT, kslice;
    size_t saveEvery;
    size_t frameId{0};

    PopulationSlice2D(
        const fs::path& out_root, std::string det_name,
        size_t NxT_, size_t NyT_, size_t NzT_, size_t kslice_,
        size_t saveEvery_
    ) : det_dir(out_root / det_name),
        NxT(NxT_), NyT(NyT_), NzT(NzT_), kslice(kslice_),
        saveEvery(saveEvery_)
    {
        std::error_code ec;
        fs::create_directories(det_dir, ec);
        if (ec) std::cerr << "[ERR] create " << det_dir << ": " << ec.message() << "\n";

        std::ofstream meta(det_dir / "metadata.json");
        if (meta) {
            meta << "{\n"
                 << "  \"NxT\": " << NxT << ",\n"
                 << "  \"NyT\": " << NyT << ",\n"
                 << "  \"NzT\": " << NzT << ",\n"
                 << "  \"kslice\": " << kslice << ",\n"
                 << "  \"saveEvery\": " << saveEvery << ",\n"
                 << "  \"dtype\": \"float64\",\n"
                 << "  \"units\": \"m^-3 (inversion density)\",\n"
                 << "  \"what\": \"population_inversion_density_xy_slice\"\n"
                 << "}\n";
        }
        std::cout << "[PopulationSlice] " << det_dir << "\n";
    }

    void save_slice(size_t n, const TwoLevelState& state) {
        if (n % saveEvery != 0) return;

        char fname[128];
        std::snprintf(fname, sizeof(fname), "inversion_%04d.raw", static_cast<int>(frameId++));

        fs::path out = det_dir / fname;
        std::ofstream ofs(out, std::ios::binary);
        if (!ofs) return;

        for (size_t i = 0; i < NxT; ++i) {
            for (size_t j = 0; j < NyT; ++j) {
                size_t id = idx3(i, j, kslice, NyT, NzT);
                double inv = static_cast<double>(state.Nu[id] - state.Ng[id]);
                ofs.write(reinterpret_cast<const char*>(&inv), sizeof(double));
            }
        }
    }
};

// Polarization diagnostics. Energy: U_int = -0.5 * integral(Pz*Ez) dV (dipole interaction energy u_int = -P*E)
struct PolarizationDiagnostics {
    fs::path det_dir;
    fs::path csv_path;
    std::ofstream f_csv;
    size_t NxT, NyT, NzT;
    size_t saveEvery;
    real dt;
    size_t i_mon, j_mon, k_mon;

    PolarizationDiagnostics(
        const fs::path& out_root, std::string det_name,
        size_t NxT_, size_t NyT_, size_t NzT_,
        size_t i_mon_, size_t j_mon_, size_t k_mon_,
        size_t saveEvery_, real dt_
    ) : det_dir(out_root / det_name),
        NxT(NxT_), NyT(NyT_), NzT(NzT_),
        i_mon(i_mon_), j_mon(j_mon_), k_mon(k_mon_),
        saveEvery(saveEvery_), dt(dt_)
    {
        std::error_code ec;
        fs::create_directories(det_dir, ec);

        csv_path = det_dir / "polarization_time.csv";
        f_csv.open(csv_path, std::ios::out);

        if (f_csv)
            f_csv << "t, Pz_center, dPz_dt_center, total_P_energy\n";

        std::cout << "[PolarizationDiag] " << det_dir << "\n";
    }

    void record(size_t n, const TwoLevelState& state,
                const std::vector<real>& Ez, const GridSpacing& grid) {
        if (n % saveEvery != 0) return;
        if (!f_csv) return;

        size_t id_mon = idx3(i_mon, j_mon, k_mon, NyT, NzT);
        real Pz_center = state.Pz[id_mon];
        real dPz_dt_center = state.dPz_dt[id_mon];

        // Dipole-field interaction energy: U_int = -0.5 * integral(Pz*Ez) dV
        // Units: [C/m^2] * [V/m] * [m^3] = [J]
        double P_energy = 0.0;
        for (size_t i = 0; i < NxT; ++i) {
            for (size_t j = 0; j < NyT; ++j) {
                for (size_t k = 0; k < NzT; ++k) {
                    size_t id = idx3(i, j, k, NyT, NzT);
                    if (state.Ndip[id] <= 0) continue;
                    real dV = grid.dx[i] * grid.dy[j] * grid.dz[k];
                    P_energy += -0.5 * state.Pz[id] * Ez[id] * dV;
                }
            }
        }

        real t = n * dt;
        f_csv << std::setprecision(18) << t << ","
              << Pz_center << "," << dPz_dt_center << ","
              << P_energy << "\n";
    }

    ~PolarizationDiagnostics() {
        if (f_csv) f_csv.close();
    }
};

// Records global population statistics over time for the gain region.
// Integrates densities with cell volume. Generates CSV and Python plotting script.
struct TLSGlobalHistory {
    fs::path out_dir;
    fs::path csv_path;
    std::ofstream f_csv;
    size_t saveEvery;
    real dt;
    bool auto_plot;
    const GridSpacing* grid_ptr;

    TLSGlobalHistory(
        const fs::path& out_root,
        const TwoLevelState& state,
        const GridSpacing& grid,
        size_t saveEvery_,
        real dt_,
        bool enable_auto_plot = false
    ) : out_dir(out_root),
        saveEvery(saveEvery_),
        dt(dt_),
        auto_plot(enable_auto_plot),
        grid_ptr(&grid)
    {
        std::error_code ec;
        fs::create_directories(out_dir, ec);
        if (ec) std::cerr << "[ERR] create " << out_dir << ": " << ec.message() << "\n";

        csv_path = out_dir / "tls_global_history.csv";
        f_csv.open(csv_path, std::ios::out);

        if (f_csv) {
            f_csv << "step,t_fs,Nu_total,Ng_total,inversion_total,Nu_percent,Ng_percent,inv_norm\n";
            f_csv.flush();
        } else {
            std::cerr << "[ERR] Failed to open " << csv_path << "\n";
        }

        std::cout << "[TLSGlobalHistory] " << csv_path << "\n";
        std::cout << "  Gain region: [" << state.gain_i0 << "," << state.gain_i1 << ") x ["
                  << state.gain_j0 << "," << state.gain_j1 << ") x ["
                  << state.gain_k0 << "," << state.gain_k1 << ")\n";
        std::cout << "  Auto-plot: " << (auto_plot ? "enabled" : "disabled") << "\n";
    }

    void record(size_t n, const TwoLevelState& state) {
        if (n % saveEvery != 0) return;
        if (!f_csv) return;

        const size_t i0 = state.gain_i0, i1 = state.gain_i1;
        const size_t j0 = state.gain_j0, j1 = state.gain_j1;
        const size_t k0 = state.gain_k0, k1 = state.gain_k1;

        double Nu_total = 0.0, Ng_total = 0.0;
        for (size_t i = i0; i < i1; ++i) {
            for (size_t j = j0; j < j1; ++j) {
                for (size_t k = k0; k < k1; ++k) {
                    size_t id = idx3(i, j, k, state.NyT, state.NzT);
                    real dV = grid_ptr->dx[i] * grid_ptr->dy[j] * grid_ptr->dz[k];
                    Nu_total += state.Nu[id] * dV;
                    Ng_total += state.Ng[id] * dV;
                }
            }
        }

        double Ntotal = Nu_total + Ng_total;
        double inversion_total = Nu_total - Ng_total;

        double Nu_percent = 0.0, Ng_percent = 0.0;
        double inv_norm = 0.5;  // Default: equilibrium
        if (Ntotal > 1e-30) {
            Nu_percent = Nu_total / Ntotal;
            Ng_percent = Ng_total / Ntotal;
            inv_norm = (inversion_total / Ntotal + 1.0) / 2.0;  // 1=all Nu, 0=all Ng, 0.5=equilibrium
        }

        double t_fs = n * dt * 1e15;
        f_csv << n << ","
              << std::setprecision(6) << std::fixed << t_fs << ","
              << std::setprecision(12) << std::scientific
              << Nu_total << "," << Ng_total << "," << inversion_total << ","
              << std::setprecision(8) << std::fixed
              << Nu_percent << "," << Ng_percent << "," << inv_norm << "\n";
        f_csv.flush();
    }

    void write_plot_script() {
        fs::path script_path = out_dir / "plot_tls_global_history.py";
        std::ofstream f_py(script_path);

        if (!f_py) {
            std::cerr << "[ERR] Failed to create " << script_path << "\n";
            return;
        }

        f_py << R"(#!/usr/bin/env python3
"""
TLS Global Population History Plotter

Population data (Nu_total, Ng_total) are integrated values: integral N(r) dV [atoms]
This ensures grid-independent results.

Generates:
  1. tls_Nu_Ng.png           - Nu_total(t) and Ng_total(t)
  2. tls_inversion.png       - inversion_total(t) = Nu - Ng
  3. tls_percent_and_norm.png - population fractions vs time
"""

import csv
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

script_dir = os.path.dirname(os.path.abspath(__file__))
csv_path = os.path.join(script_dir, 'tls_global_history.csv')

print(f"Reading {csv_path}...")
data = {
    'step': [], 't_fs': [], 'Nu_total': [], 'Ng_total': [],
    'inversion_total': [], 'Nu_percent': [], 'Ng_percent': [], 'inv_norm': []
}

with open(csv_path, 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        data['step'].append(int(row['step']))
        data['t_fs'].append(float(row['t_fs']))
        data['Nu_total'].append(float(row['Nu_total']))
        data['Ng_total'].append(float(row['Ng_total']))
        data['inversion_total'].append(float(row['inversion_total']))
        data['Nu_percent'].append(float(row['Nu_percent']))
        data['Ng_percent'].append(float(row['Ng_percent']))
        data['inv_norm'].append(float(row['inv_norm']))

for key in data:
    data[key] = np.array(data[key])

n_points = len(data['step'])
print(f"Loaded {n_points} data points")
print(f"Time range: {data['t_fs'].min():.2f} - {data['t_fs'].max():.2f} fs")
sum_pct = data['Nu_percent'] + data['Ng_percent']
print(f"Nu_percent + Ng_percent check: min={sum_pct.min():.6f}, max={sum_pct.max():.6f}")
print(f"inv_norm range: [{data['inv_norm'].min():.6f}, {data['inv_norm'].max():.6f}]")

t = data['t_fs']

# Plot 1: Nu_total and Ng_total
fig1, ax1 = plt.subplots(figsize=(10, 6))
ax1.plot(t, data['Nu_total'], 'b-', linewidth=1.5, label='$N_u$ (upper level)')
ax1.plot(t, data['Ng_total'], 'r-', linewidth=1.5, label='$N_g$ (ground level)')
ax1.set_xlabel('Time (fs)', fontsize=12)
ax1.set_ylabel('Integrated Population (atoms)', fontsize=12)
ax1.set_title('Two-Level System: Global Population Dynamics', fontsize=14)
ax1.legend(loc='best', fontsize=11)
ax1.grid(True, alpha=0.3)
ax1.ticklabel_format(style='scientific', axis='y', scilimits=(0,0))
fig1.tight_layout()
out1 = os.path.join(script_dir, 'tls_Nu_Ng.png')
fig1.savefig(out1, dpi=150)
print(f"Saved: {out1}")
plt.close(fig1)

# Plot 2: Inversion
fig2, ax2 = plt.subplots(figsize=(10, 6))
inv = data['inversion_total']
ax2.plot(t, inv, 'g-', linewidth=1.5)
ax2.axhline(y=0, color='k', linestyle='--', linewidth=0.8, alpha=0.5)
ax2.fill_between(t, inv, 0, where=(inv > 0), color='green', alpha=0.2, label='Gain (Nu > Ng)')
ax2.fill_between(t, inv, 0, where=(inv < 0), color='red', alpha=0.2, label='Absorption (Ng > Nu)')
ax2.set_xlabel('Time (fs)', fontsize=12)
ax2.set_ylabel('Integrated Inversion (atoms)', fontsize=12)
ax2.set_title('Two-Level System: Population Inversion', fontsize=14)
ax2.legend(loc='best', fontsize=11)
ax2.grid(True, alpha=0.3)
ax2.ticklabel_format(style='scientific', axis='y', scilimits=(0,0))
fig2.tight_layout()
out2 = os.path.join(script_dir, 'tls_inversion.png')
fig2.savefig(out2, dpi=150)
print(f"Saved: {out2}")
plt.close(fig2)

# Plot 3: Percentages and Normalized Inversion
fig3, ax3 = plt.subplots(figsize=(10, 6))
ax3.plot(t, data['Nu_percent'], 'b-', linewidth=1.5, label='$N_u / N_{total}$')
ax3.plot(t, data['Ng_percent'], 'r-', linewidth=1.5, label='$N_g / N_{total}$')
ax3.plot(t, data['inv_norm'], 'g--', linewidth=1.5, label='inv_norm (1=all Nu, 0=all Ng)')
ax3.axhline(y=0.5, color='k', linestyle=':', linewidth=0.8, alpha=0.5, label='Equilibrium (0.5)')
ax3.set_xlabel('Time (fs)', fontsize=12)
ax3.set_ylabel('Fractional Population', fontsize=12)
ax3.set_title('Two-Level System: Normalized Population Fractions', fontsize=14)
ax3.set_ylim(-0.05, 1.05)
ax3.legend(loc='best', fontsize=10)
ax3.grid(True, alpha=0.3)
fig3.tight_layout()
out3 = os.path.join(script_dir, 'tls_percent_and_norm.png')
fig3.savefig(out3, dpi=150)
print(f"Saved: {out3}")
plt.close(fig3)

print("Done!")
)";

        f_py.close();
        std::cout << "[TLSGlobalHistory] Python script: " << script_path << "\n";
    }

    void generate_plots() {
        write_plot_script();

        if (!auto_plot) {
            fs::path script_path = out_dir / "plot_tls_global_history.py";
            std::cout << "[TLSGlobalHistory] To generate plots, run:\n";
            std::cout << "  python3 " << script_path << "\n";
            return;
        }

        fs::path script_path = out_dir / "plot_tls_global_history.py";
        std::string cmd = "python3 " + script_path.string() + " 2>&1";
        std::cout << "[TLSGlobalHistory] Running: " << cmd << "\n";
        int ret = std::system(cmd.c_str());

        if (ret == 0) {
            std::cout << "[TLSGlobalHistory] Plots generated successfully\n";
        } else {
            std::cerr << "[TLSGlobalHistory] Plot generation failed (ret=" << ret << ")\n";
            std::cerr << "  You can run manually: python3 " << script_path << "\n";
        }
    }

    void finalize() {
        if (f_csv) f_csv.close();
        generate_plots();
    }

    ~TLSGlobalHistory() {
        if (f_csv) f_csv.close();
    }
};
