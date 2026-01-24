// detectors.hpp —— Detectors

#pragma once

#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cstdio>       // std::snprintf
#include <iomanip>      // std::setprecision
#include <cstddef>
#include <cassert>

#include "global_function.hpp"

namespace fs = std::filesystem;

// ======================= 2D slice (XY @ fixed k) frame sequence =======================
struct Detector2D {
    fs::path det_dir;                       // Output directory
    std::size_t NxT{}, NyT{}, NzT{}, kslice{};
    std::size_t saveEvery{}, frameId{ 0 };

    // Only for writing metadata (to help Python restore coordinates from physical region/PML)
    std::size_t Nx_phys{}, Ny_phys{}, Nz_phys{}, npml{};

    // Physical coordinate information
    real pml_offset_x{}, pml_offset_y{}, pml_offset_z{};  // PML offsets (meters)
    real z_slice_physical{};  // Z-coordinate of slice in physical coords (meters)

    // Physical domain size (user-configured, for visualization)
    real Lx_phys{}, Ly_phys{}, Lz_phys{};  // Physical domain dimensions (meters)

    // Output control
    std::string framePattern{ "ez_%04d.raw" };
    bool writeFloat64{ true };  // true=write float64(double), false=write float32(float)

    Detector2D(const fs::path& out_root,
        std::string det_name,
        std::size_t NxT_, std::size_t NyT_, std::size_t NzT_, std::size_t kslice_,
        std::size_t saveEvery_, std::size_t nSteps_,
        std::size_t Nx_phys_, std::size_t Ny_phys_, std::size_t Nz_phys_,
        std::size_t npml_,
        const GridSpacing& grid,  // Grid spacing for physical coordinate calculation
        real Lx_phys_ = 0, real Ly_phys_ = 0, real Lz_phys_ = 0,  // User-configured physical domain size
        std::string framePattern_ = "ez_%04d.raw",
        bool writeFloat64_ = true)
        : det_dir(out_root / det_name),
        NxT(NxT_), NyT(NyT_), NzT(NzT_), kslice(kslice_),
        saveEvery(saveEvery_),
        Nx_phys(Nx_phys_), Ny_phys(Ny_phys_), Nz_phys(Nz_phys_), npml(npml_),
        Lx_phys(Lx_phys_), Ly_phys(Ly_phys_), Lz_phys(Lz_phys_),
        framePattern(std::move(framePattern_)),
        writeFloat64(writeFloat64_)
    {
        // Calculate physical coordinate offsets
        pml_offset_x = grid.pml_offset_x();
        pml_offset_y = grid.pml_offset_y();
        pml_offset_z = grid.pml_offset_z();
        // Z-coordinate of slice in physical coordinates
        z_slice_physical = grid.cell_center_physical_z(kslice);

        std::error_code ec;
        fs::create_directories(det_dir, ec);
        if (ec) {
            std::cerr << "[ERR] create " << det_dir << ": " << ec.message() << "\n";
        }
        write_metadata_(nSteps_);
        std::cout << "[Detector2D] " << det_dir << "\n";
    }

    // Save one frame every saveEvery steps
    template<typename Real>
    void try_save(std::size_t n, const std::vector<Real>& fieldVecXYZ) {
        if (n % saveEvery) return;
        char fname[128];
        std::snprintf(fname, sizeof(fname), framePattern.c_str(), static_cast<int>(frameId++));
        fs::path out = det_dir / fname;
        save_slice_xy_at_k_(fieldVecXYZ, out.string());
    }

private:
    void write_metadata_(std::size_t nSteps_) const {
        const fs::path json_path = det_dir / "metadata.json";
        std::ofstream ofs(json_path, std::ios::binary);
        if (!ofs) {
            std::cerr << "[ERR] open " << json_path << " failed.\n";
            return;
        }
        ofs << std::setprecision(15);
        ofs << "{\n"
            << "  \"NxT\": " << NxT << ",\n"
            << "  \"NyT\": " << NyT << ",\n"
            << "  \"NzT\": " << NzT << ",\n"
            << "  \"Nx\": " << Nx_phys << ",\n"
            << "  \"Ny\": " << Ny_phys << ",\n"
            << "  \"Nz\": " << Nz_phys << ",\n"
            << "  \"npml\": " << npml << ",\n"
            << "  \"kslice\": " << kslice << ",\n"
            << "  \"saveEvery\": " << saveEvery << ",\n"
            << "  \"nSteps\": " << nSteps_ << ",\n"
            << "  \"framePattern\": \"" << framePattern << "\",\n"
            << "  \"dtype\": \"" << (writeFloat64 ? "float64" : "float32") << "\",\n"
            << "  \"pml_offset_x_m\": " << pml_offset_x << ",\n"
            << "  \"pml_offset_y_m\": " << pml_offset_y << ",\n"
            << "  \"pml_offset_z_m\": " << pml_offset_z << ",\n"
            << "  \"z_slice_physical_m\": " << z_slice_physical << ",\n"
            << "  \"z_slice_physical_nm\": " << z_slice_physical * 1e9 << ",\n"
            << "  \"Lx_phys_m\": " << Lx_phys << ",\n"
            << "  \"Ly_phys_m\": " << Ly_phys << ",\n"
            << "  \"Lz_phys_m\": " << Lz_phys << ",\n"
            << "  \"Lx_phys_nm\": " << Lx_phys * 1e9 << ",\n"
            << "  \"Ly_phys_nm\": " << Ly_phys * 1e9 << ",\n"
            << "  \"Lz_phys_nm\": " << Lz_phys * 1e9 << ",\n"
            << "  \"note\": \"row-major (i then j) for XY plane at fixed k=kslice\",\n"
            << "  \"coord_note\": \"Physical coords: origin at physical domain boundary (after PML)\"\n"
            << "}\n";
    }

    // Write XY@kslice slice
    template<typename Real>
    void save_slice_xy_at_k_(const std::vector<Real>& F,
        const std::string& filename) const
    {
        std::ofstream ofs(filename, std::ios::binary);
        if (!ofs) return;

        if (writeFloat64) {
            // Write double (float64)
            for (std::size_t i = 0; i < NxT; ++i)
                for (std::size_t j = 0; j < NyT; ++j) {
                    const double v = static_cast<double>(F[idx3(i, j, kslice, NyT, NzT)]);
                    ofs.write(reinterpret_cast<const char*>(&v), sizeof(double));
                }
        }
        else {
            // Write float (float32)
            for (std::size_t i = 0; i < NxT; ++i)
                for (std::size_t j = 0; j < NyT; ++j) {
                    const float v = static_cast<float>(F[idx3(i, j, kslice, NyT, NzT)]);
                    ofs.write(reinterpret_cast<const char*>(&v), sizeof(float));
                }
        }
    }
};

// ======================= Single point Ez time series detector (center or specified point) =======================
struct Probe1D {
    fs::path det_dir;      // e.g. frames/<run_tag>/Ez_probe
    fs::path bin_path;     // ez_center_ts.bin
    fs::path meta_path;    // metadata_ts.json

    std::size_t NxT{}, NyT{}, NzT{};
    std::size_t i0{}, j0{}, k0{};
    std::size_t saveEvery{};
    std::size_t nSteps{};
    real dt{};

    // Physical coordinate information
    real probe_x_physical{}, probe_y_physical{}, probe_z_physical{};
    real pml_offset_x{}, pml_offset_y{}, pml_offset_z{};

    std::ofstream ofs;     // Binary write handle (float64)

    Probe1D(const fs::path& out_root,
        std::string det_name,
        std::size_t NxT_, std::size_t NyT_, std::size_t NzT_,
        std::size_t i0_, std::size_t j0_, std::size_t k0_,
        std::size_t saveEvery_, std::size_t nSteps_, real dt_,
        const GridSpacing& grid)  // Grid spacing for physical coordinate calculation
        : det_dir(out_root / det_name),
        NxT(NxT_), NyT(NyT_), NzT(NzT_),
        i0(i0_), j0(j0_), k0(k0_),
        saveEvery(saveEvery_), nSteps(nSteps_), dt(dt_)
    {
        // Calculate physical coordinates
        pml_offset_x = grid.pml_offset_x();
        pml_offset_y = grid.pml_offset_y();
        pml_offset_z = grid.pml_offset_z();
        probe_x_physical = grid.cell_center_physical_x(i0);
        probe_y_physical = grid.cell_center_physical_y(j0);
        probe_z_physical = grid.cell_center_physical_z(k0);

        std::error_code ec;
        fs::create_directories(det_dir, ec);
        if (ec) std::cerr << "[ERR] create " << det_dir << ": " << ec.message() << "\n";

        bin_path = det_dir / "ez_center_ts.bin";
        meta_path = det_dir / "metadata_ts.json";

        ofs.open(bin_path, std::ios::binary);
        if (!ofs) std::cerr << "[ERR] open " << bin_path << " failed.\n";

        // Simple metadata for Python to read and plot
        std::ofstream m(meta_path, std::ios::binary);
        if (m) {
            m << std::setprecision(15);
            m << "{\n"
                << "  \"NxT\": " << NxT << ",\n"
                << "  \"NyT\": " << NyT << ",\n"
                << "  \"NzT\": " << NzT << ",\n"
                << "  \"i0\": " << i0 << ",\n"
                << "  \"j0\": " << j0 << ",\n"
                << "  \"k0\": " << k0 << ",\n"
                << "  \"saveEvery\": " << saveEvery << ",\n"
                << "  \"nSteps\": " << nSteps << ",\n"
                << "  \"dt\": " << std::setprecision(18) << dt << ",\n"
                << "  \"dtype\": \"float64\",\n"
                << "  \"pml_offset_x_m\": " << pml_offset_x << ",\n"
                << "  \"pml_offset_y_m\": " << pml_offset_y << ",\n"
                << "  \"pml_offset_z_m\": " << pml_offset_z << ",\n"
                << "  \"probe_x_physical_m\": " << probe_x_physical << ",\n"
                << "  \"probe_y_physical_m\": " << probe_y_physical << ",\n"
                << "  \"probe_z_physical_m\": " << probe_z_physical << ",\n"
                << "  \"probe_x_physical_nm\": " << probe_x_physical * 1e9 << ",\n"
                << "  \"probe_y_physical_nm\": " << probe_y_physical * 1e9 << ",\n"
                << "  \"probe_z_physical_nm\": " << probe_z_physical * 1e9 << ",\n"
                << "  \"note\": \"time series of Ez(i0,j0,k0); stored every saveEvery steps\",\n"
                << "  \"coord_note\": \"Physical coords: origin at physical domain boundary (after PML)\"\n"
                << "}\n";
        }
        else {
            std::cerr << "[ERR] open " << meta_path << " failed.\n";
        }

        std::cout << "[Probe1D] " << det_dir << "\n";
    }

    template<typename Real>
    void try_record(std::size_t n, const std::vector<Real>& Ez) {
        if (!ofs) return;
        if (n % saveEvery) return;
        const std::size_t id = idx3(i0, j0, k0, NyT, NzT);
        const real v = static_cast<real>(Ez[id]);        // Always store as float64
        ofs.write(reinterpret_cast<const char*>(&v), sizeof(real));
    }

    ~Probe1D() {
        if (ofs) ofs.close();
    }
};

// ======================= 2D slice refractive index n(x,y) @ fixed k =======================
struct RefrIndex2D {
    fs::path det_dir;                                   // Output directory: frames/<run_tag>/n_center
    std::size_t NxT{}, NyT{}, NzT{}, kslice{};
    std::size_t Nx_phys{}, Ny_phys{}, Nz_phys{}, npml{};
    std::string framePattern{ "n_%04d.raw" };
    bool writeFloat64{ true };
    std::size_t frameId{ 0 };

    // Physical coordinate information
    real pml_offset_x{}, pml_offset_y{}, pml_offset_z{};
    real z_slice_physical{};

    // Physical domain size (user-configured, for visualization)
    real Lx_phys{}, Ly_phys{}, Lz_phys{};

    RefrIndex2D(const fs::path& out_root,
        std::string det_name,       // Suggested "n_center"
        std::size_t NxT_, std::size_t NyT_, std::size_t NzT_, std::size_t kslice_,
        std::size_t Nx_phys_, std::size_t Ny_phys_, std::size_t Nz_phys_,
        std::size_t npml_,
        const GridSpacing& grid,  // Grid spacing for physical coordinate calculation
        real Lx_phys_ = 0, real Ly_phys_ = 0, real Lz_phys_ = 0,  // User-configured physical domain size
        std::string framePattern_ = "n_%04d.raw",
        bool writeFloat64_ = true)
        : det_dir(out_root / det_name),
        NxT(NxT_), NyT(NyT_), NzT(NzT_), kslice(kslice_),
        Nx_phys(Nx_phys_), Ny_phys(Ny_phys_), Nz_phys(Nz_phys_), npml(npml_),
        Lx_phys(Lx_phys_), Ly_phys(Ly_phys_), Lz_phys(Lz_phys_),
        framePattern(std::move(framePattern_)), writeFloat64(writeFloat64_)
    {
        // Calculate physical coordinate offsets
        pml_offset_x = grid.pml_offset_x();
        pml_offset_y = grid.pml_offset_y();
        pml_offset_z = grid.pml_offset_z();
        z_slice_physical = grid.cell_center_physical_z(kslice);

        std::error_code ec;
        fs::create_directories(det_dir, ec);
        if (ec) std::cerr << "[ERR] create " << det_dir << ": " << ec.message() << "\n";
        write_metadata_();
        std::cout << "[RefrIndex2D] " << det_dir << "\n";
    }

    // Extract z=kslice slice from a 3D scalar field (size NxT*NyT*NzT) and save as one frame
    template<typename Real>
    void save_from_scalar3d(const std::vector<Real>& scalar3d) {
        if (scalar3d.size() != NxT * NyT * NzT) {
            std::cerr << "[ERR] n_grid size mismatch.\n";
            return;
        }
        char fname[128];
        std::snprintf(fname, sizeof(fname), framePattern.c_str(), static_cast<int>(frameId++));
        fs::path out = det_dir / fname;

        std::ofstream ofs(out, std::ios::binary);
        if (!ofs) {
            std::cerr << "[ERR] open " << out << " failed.\n";
            return;
        }

        if (writeFloat64) {
            for (std::size_t i = 0; i < NxT; ++i)
                for (std::size_t j = 0; j < NyT; ++j) {
                    const double v = static_cast<double>(scalar3d[idx3(i, j, kslice, NyT, NzT)]);
                    ofs.write(reinterpret_cast<const char*>(&v), sizeof(double));
                }
        }
        else {
            for (std::size_t i = 0; i < NxT; ++i)
                for (std::size_t j = 0; j < NyT; ++j) {
                    const float v = static_cast<float>(scalar3d[idx3(i, j, kslice, NyT, NzT)]);
                    ofs.write(reinterpret_cast<const char*>(&v), sizeof(float));
                }
        }
        // std::cout << "[SAVE n] " << fs::absolute(out) << "\n";
    }

private:
    void write_metadata_() const {
        const fs::path json_path = det_dir / "metadata.json";
        std::ofstream ofs(json_path, std::ios::binary);
        if (!ofs) {
            std::cerr << "[ERR] open " << json_path << " failed.\n";
            return;
        }
        ofs << std::setprecision(15);
        ofs << "{\n"
            << "  \"NxT\": " << NxT << ",\n"
            << "  \"NyT\": " << NyT << ",\n"
            << "  \"NzT\": " << NzT << ",\n"
            << "  \"Nx\": " << Nx_phys << ",\n"
            << "  \"Ny\": " << Ny_phys << ",\n"
            << "  \"Nz\": " << Nz_phys << ",\n"
            << "  \"npml\": " << npml << ",\n"
            << "  \"kslice\": " << kslice << ",\n"
            << "  \"framePattern\": \"" << framePattern << "\",\n"
            << "  \"dtype\": \"" << (writeFloat64 ? "float64" : "float32") << "\",\n"
            << "  \"pml_offset_x_m\": " << pml_offset_x << ",\n"
            << "  \"pml_offset_y_m\": " << pml_offset_y << ",\n"
            << "  \"pml_offset_z_m\": " << pml_offset_z << ",\n"
            << "  \"z_slice_physical_m\": " << z_slice_physical << ",\n"
            << "  \"z_slice_physical_nm\": " << z_slice_physical * 1e9 << ",\n"
            << "  \"Lx_phys_m\": " << Lx_phys << ",\n"
            << "  \"Ly_phys_m\": " << Ly_phys << ",\n"
            << "  \"Lz_phys_m\": " << Lz_phys << ",\n"
            << "  \"Lx_phys_nm\": " << Lx_phys * 1e9 << ",\n"
            << "  \"Ly_phys_nm\": " << Ly_phys * 1e9 << ",\n"
            << "  \"Lz_phys_nm\": " << Lz_phys * 1e9 << ",\n"
            << "  \"what\": \"refractive_index_n_xy_at_fixed_k\",\n"
            << "  \"note\": \"row-major (i then j) for XY plane at fixed k=kslice\",\n"
            << "  \"coord_note\": \"Physical coords: origin at physical domain boundary (after PML)\"\n"
            << "}\n";
    }
};

// ======================= Box energy flux detector (time-centered version) =======================

#ifndef DETECTORS_BOX_POYNTING_HPP
#define DETECTORS_BOX_POYNTING_HPP

struct BoxPoyntingDetector {
    fs::path det_dir;                // Output directory: frames/<run_tag>/<det_name>
    std::size_t NxT{}, NyT{}, NzT{}; // Total grid
    // Box interior voxel closed interval [i0,i1]×[j0,j1]×[k0,k1]
    std::size_t i0{}, i1{}, j0{}, j1{}, k0{}, k1{};
    std::size_t saveEvery{}, nSteps{};
    real dt{};
    std::vector<real> dx_array_, dy_array_, dz_array_;

    // Cumulative energy (J)
    real Eout{ 0 };

    // Output files
    fs::path power_csv;
    fs::path energy_csv;
    std::ofstream f_power, f_energy;

    bool has_prev{ false };
    std::vector<real> Eprev_x, Eprev_y, Eprev_z;
    std::vector<real> Hprev_x, Hprev_y, Hprev_z;

    BoxPoyntingDetector(const fs::path& out_root,
        std::string det_name,
        std::size_t NxT_, std::size_t NyT_, std::size_t NzT_,
        std::size_t i0_, std::size_t i1_,
        std::size_t j0_, std::size_t j1_,
        std::size_t k0_, std::size_t k1_,
        std::size_t saveEvery_, std::size_t nSteps_,
        const GridSpacing& grid_spacing, real dt_)
        : det_dir(out_root / det_name),
        NxT(NxT_), NyT(NyT_), NzT(NzT_),
        i0(i0_), i1(i1_), j0(j0_), j1(j1_), k0(k0_), k1(k1_),
        saveEvery(saveEvery_), nSteps(nSteps_), dt(dt_)
    {
        dx_array_ = grid_spacing.dx;
        dy_array_ = grid_spacing.dy;
        dz_array_ = grid_spacing.dz; 

        std::error_code ec;
        fs::create_directories(det_dir, ec);
        if (ec) std::cerr << "[ERR] create " << det_dir << ": " << ec.message() << "\n";

        // Boundary/out-of-bounds protection
        assert(i0 >= 1 && j0 >= 1 && k0 >= 1);
        assert(i1 + 1 < NxT && j1 + 1 < NyT && k1 + 1 < NzT);
        assert(i0 <= i1 && j0 <= j1 && k0 <= k1);

        power_csv = det_dir / "power_time.csv";
        energy_csv = det_dir / "energy_time.csv";
        f_power.open(power_csv, std::ios::out | std::ios::trunc);
        f_energy.open(energy_csv, std::ios::out | std::ios::trunc);
        if (!f_power)  std::cerr << "[ERR] open " << power_csv << " failed.\n";
        if (!f_energy) std::cerr << "[ERR] open " << energy_csv << " failed.\n";

        if (f_power)  f_power << "t, Px_pos, Px_neg, Py_pos, Py_neg, Pz_pos, Pz_neg, P_total\n";
        if (f_energy) f_energy << "t, E_out_cumulative\n";

        std::ofstream meta(det_dir / "metadata.json", std::ios::binary);
        if (meta) {
            meta << "{\n"
                << "  \"NxT\": " << NxT << ",\n"
                << "  \"NyT\": " << NyT << ",\n"
                << "  \"NzT\": " << NzT << ",\n"
                << "  \"i0\": " << i0 << ", \"i1\": " << i1 << ",\n"
                << "  \"j0\": " << j0 << ", \"j1\": " << j1 << ",\n"
                << "  \"k0\": " << k0 << ", \"k1\": " << k1 << ",\n"
                << "  \"saveEvery\": " << saveEvery << ",\n"
                << "  \"nSteps\": " << nSteps << ",\n"
                // << "  \"dx\": " << std::setprecision(18) << dx << ",\n"
                // << "  \"dy\": " << std::setprecision(18) << dy << ",\n"
                // << "  \"dz\": " << std::setprecision(18) << dz << ",\n"
                << "  \"note_spacing\": \"non-uniform grid, see spacing arrays\",\n"
                << "  \"dt\": " << std::setprecision(18) << dt << ",\n"
                << "  \"what\": \"time-centered power through 6 faces and cumulative outward energy\",\n"
                << "  \"note\": \"+X/-X/+Y/-Y/+Z/-Z faces; outward normal positive; time-centered\"\n"
                << "}\n";
        }

        std::cout << "[BoxPoynting] " << det_dir << " (time-centered)\n";
    }

    template<typename Real>
    void try_record(std::size_t n,
        const std::vector<Real>& Ex,
        const std::vector<Real>& Ey,
        const std::vector<Real>& Ez,
        const std::vector<Real>& Hx,
        const std::vector<Real>& Hy,
        const std::vector<Real>& Hz)
    {
        // --- Calculate "current" six-face power ---
        real Px_pos_now, Px_neg_now, Py_pos_now, Py_neg_now, Pz_pos_now, Pz_neg_now, Ptot_now;
        compute_flux_all(Ex, Ey, Ez, Hx, Hy, Hz,
            Px_pos_now, Px_neg_now, Py_pos_now, Py_neg_now, Pz_pos_now, Pz_neg_now, Ptot_now);

        // --- If "previous time step" exists, do time-centered average ---
        real Px_pos_out = Px_pos_now, Px_neg_out = Px_neg_now;
        real Py_pos_out = Py_pos_now, Py_neg_out = Py_neg_now;
        real Pz_pos_out = Pz_pos_now, Pz_neg_out = Pz_neg_now;
        real Ptot_out = Ptot_now;

        if (has_prev) {
            real Px_pos_prev, Px_neg_prev, Py_pos_prev, Py_neg_prev, Pz_pos_prev, Pz_neg_prev, Ptot_prev;
            compute_flux_all(Eprev_x, Eprev_y, Eprev_z, Hprev_x, Hprev_y, Hprev_z,
                Px_pos_prev, Px_neg_prev, Py_pos_prev, Py_neg_prev, Pz_pos_prev, Pz_neg_prev, Ptot_prev);
            Px_pos_out = static_cast<real>(0.5) * (Px_pos_now + Px_pos_prev);
            Px_neg_out = static_cast<real>(0.5) * (Px_neg_now + Px_neg_prev);
            Py_pos_out = static_cast<real>(0.5) * (Py_pos_now + Py_pos_prev);
            Py_neg_out = static_cast<real>(0.5) * (Py_neg_now + Py_neg_prev);
            Pz_pos_out = static_cast<real>(0.5) * (Pz_pos_now + Pz_pos_prev);
            Pz_neg_out = static_cast<real>(0.5) * (Pz_neg_now + Pz_neg_prev);
            Ptot_out = static_cast<real>(0.5) * (Ptot_now + Ptot_prev);
        }

        // --- Write only when needed, but "energy integral" still proceeds at sampling interval ---
        if (f_power && (n % saveEvery == 0)) {
            real t = n * dt;
            f_power << t << "," << Px_pos_out << "," << Px_neg_out << ","
                << Py_pos_out << "," << Py_neg_out << ","
                << Pz_pos_out << "," << Pz_neg_out << "," << Ptot_out << "\n";
        }
        if (f_energy && (n % saveEvery == 0)) {
            real dt_eff = static_cast<real>(saveEvery) * dt; // Approximate integral with recording interval
            Eout += Ptot_out * dt_eff;
            real t = n * dt;
            f_energy << t << "," << Eout << "\n";
        }

        // --- Update cache (update on every call to ensure next time can be centered) ---
        if (!has_prev) {
            Eprev_x.resize(Ex.size()); Eprev_y.resize(Ey.size()); Eprev_z.resize(Ez.size());
            Hprev_x.resize(Hx.size()); Hprev_y.resize(Hy.size()); Hprev_z.resize(Hz.size());
        }
        // Copy current fields to prev (note: cost is proportional to domain size, but only one memory copy)
        std::copy(Ex.begin(), Ex.end(), Eprev_x.begin());
        std::copy(Ey.begin(), Ey.end(), Eprev_y.begin());
        std::copy(Ez.begin(), Ez.end(), Eprev_z.begin());
        std::copy(Hx.begin(), Hx.end(), Hprev_x.begin());
        std::copy(Hy.begin(), Hy.end(), Hprev_y.begin());
        std::copy(Hz.begin(), Hz.end(), Hprev_z.begin());
        has_prev = true;
    }

private:
    // Calculate "six-face power and total power" once
    template<typename RX, typename RY, typename RZ>
    void compute_flux_all(const std::vector<RX>& Ex,
        const std::vector<RY>& Ey,
        const std::vector<RZ>& Ez,
        const std::vector<RX>& Hx,
        const std::vector<RY>& Hy,
        const std::vector<RZ>& Hz,
        real& Px_pos, real& Px_neg,
        real& Py_pos, real& Py_neg,
        real& Pz_pos, real& Pz_neg,
        real& Ptot) const
    {
        Px_pos = flux_x_face(i1 + 1, +1, Ex, Ey, Ez, Hx, Hy, Hz);
        Px_neg = flux_x_face(i0, -1, Ex, Ey, Ez, Hx, Hy, Hz);
        Py_pos = flux_y_face(j1 + 1, +1, Ex, Ey, Ez, Hx, Hy, Hz);
        Py_neg = flux_y_face(j0, -1, Ex, Ey, Ez, Hx, Hy, Hz);
        Pz_pos = flux_z_face(k1 + 1, +1, Ex, Ey, Ez, Hx, Hy, Hz);
        Pz_neg = flux_z_face(k0, -1, Ex, Ey, Ez, Hx, Hy, Hz);
        Ptot = Px_pos + Px_neg + Py_pos + Py_neg + Pz_pos + Pz_neg;
    }

    // The following three flux_* functions are consistent with original (only type templates made more generic)

    template<typename RX, typename RY, typename RZ>
    real flux_x_face(std::size_t iface, int sign,
        const std::vector<RX>& Ex,
        const std::vector<RY>& Ey,
        const std::vector<RZ>& Ez,
        const std::vector<RX>& Hx,
        const std::vector<RY>& Hy,
        const std::vector<RZ>& Hz) const
    {
        const std::size_t iL = iface - 1, iR = iface;
        real P = 0;
        for (std::size_t j = j0; j <= j1; ++j)
            for (std::size_t k = k0; k <= k1; ++k) {
                real dy_local = dy_array_[j];
                real dz_local = dz_array_[k];

                std::size_t idL = idx3(iL, j, k, NyT, NzT);
                std::size_t idR = idx3(iR, j, k, NyT, NzT);

                real Ey_c = 0.5 * (Ey[idL] + Ey[idR]);
                real Ez_c = 0.5 * (Ez[idL] + Ez[idR]);
                real Hy_c = 0.5 * (Hy[idL] + Hy[idR]);
                real Hz_c = 0.5 * (Hz[idL] + Hz[idR]);

                real Sx = Ey_c * Hz_c - Ez_c * Hy_c;   // W/m^2
                P += Sx * dy_local * dz_local * sign;
            }
        return P;
    }

    template<typename RX, typename RY, typename RZ>
    real flux_y_face(std::size_t jface, int sign,
        const std::vector<RX>& Ex,
        const std::vector<RY>& Ey,
        const std::vector<RZ>& Ez,
        const std::vector<RX>& Hx,
        const std::vector<RY>& Hy,
        const std::vector<RZ>& Hz) const
    {
        const std::size_t jB = jface - 1, jF = jface;
        real P = 0;
        for (std::size_t i = i0; i <= i1; ++i)
            for (std::size_t k = k0; k <= k1; ++k) {
                real dx_local = dx_array_[i];
                real dz_local = dz_array_[k];

                std::size_t idB = idx3(i, jB, k, NyT, NzT);
                std::size_t idF = idx3(i, jF, k, NyT, NzT);

                real Ex_c = 0.5 * (Ex[idB] + Ex[idF]);
                real Ez_c = 0.5 * (Ez[idB] + Ez[idF]);
                real Hx_c = 0.5 * (Hx[idB] + Hx[idF]);
                real Hz_c = 0.5 * (Hz[idB] + Hz[idF]);

                real Sy = Ez_c * Hx_c - Ex_c * Hz_c;
                P += Sy * dx_local * dz_local * sign; 
            }
        return P;
    }

    template<typename RX, typename RY, typename RZ>
    real flux_z_face(std::size_t kface, int sign,
        const std::vector<RX>& Ex,
        const std::vector<RY>& Ey,
        const std::vector<RZ>& Ez,
        const std::vector<RX>& Hx,
        const std::vector<RY>& Hy,
        const std::vector<RZ>& Hz) const
    {
        const std::size_t kD = kface - 1, kU = kface;
        real P = 0;
        for (std::size_t i = i0; i <= i1; ++i)
            for (std::size_t j = j0; j <= j1; ++j) {
                real dy_local = dy_array_[j];
                real dx_local = dx_array_[i];

                std::size_t idD = idx3(i, j, kD, NyT, NzT);
                std::size_t idU = idx3(i, j, kU, NyT, NzT);

                real Ey_c = 0.5 * (Ey[idD] + Ey[idU]);
                real Ex_c = 0.5 * (Ex[idD] + Ex[idU]);
                real Hy_c = 0.5 * (Hy[idD] + Hy[idU]);
                real Hx_c = 0.5 * (Hx[idD] + Hx[idU]);

                real Sz = Ex_c * Hy_c - Ey_c * Hx_c;
                P += Sz * dy_local * dx_local * sign;
            }
        return P;
    }
};

#endif // DETECTORS_BOX_POYNTING_HPP

// ======================= Total electromagnetic energy detector inside box =======================

#ifndef DETECTORS_BOX_ENERGY_HPP
#define DETECTORS_BOX_ENERGY_HPP

struct BoxEnergyDetector {
    fs::path det_dir;                // frames/<run_tag>/EM_total_energy
    std::size_t NxT{}, NyT{}, NzT{};
    std::size_t i0{}, i1{}, j0{}, j1{}, k0{}, k1{};
    std::size_t saveEvery{}, nSteps{};
    real dt{};
    std::vector<real> dx_array_, dy_array_, dz_array_;

    // Reference material coefficients (from MaterialGrids::bEx... etc.)
    const std::vector<real>& bEx, & bEy, & bEz;
    const std::vector<real>& bHx, & bHy, & bHz;

    // Cache recovered eps/mu (same dimensions as Yee points)
    bool coeff_ready{ false };
    std::vector<real> epsEx, epsEy, epsEz;
    std::vector<real> muHx, muHy, muHz;

    // Output
    fs::path energy_csv;
    std::ofstream f_energy;

    BoxEnergyDetector(const fs::path& out_root,
        std::size_t NxT_, std::size_t NyT_, std::size_t NzT_,
        std::size_t i0_, std::size_t i1_,
        std::size_t j0_, std::size_t j1_,
        std::size_t k0_, std::size_t k1_,
        std::size_t saveEvery_, std::size_t nSteps_,
        const GridSpacing& grid_spacing,
        real dt_,
        const std::vector<real>& bEx_,
        const std::vector<real>& bEy_,
        const std::vector<real>& bEz_,
        const std::vector<real>& bHx_,
        const std::vector<real>& bHy_,
        const std::vector<real>& bHz_)
        : det_dir(out_root / "EM_total_energy"),
        NxT(NxT_), NyT(NyT_), NzT(NzT_),
        i0(i0_), i1(i1_), j0(j0_), j1(j1_), k0(k0_), k1(k1_),
        saveEvery(saveEvery_), nSteps(nSteps_),
        dt(dt_),
        bEx(bEx_), bEy(bEy_), bEz(bEz_),
        bHx(bHx_), bHy(bHy_), bHz(bHz_)
        {
            dx_array_ = grid_spacing.dx;
            dy_array_ = grid_spacing.dy;
            dz_array_ = grid_spacing.dz;
            std::error_code ec;
            fs::create_directories(det_dir, ec);
            if (ec) std::cerr << "[ERR] create " << det_dir << ": " << ec.message() << "\n";

            energy_csv = det_dir / "energy_time.csv";
            f_energy.open(energy_csv, std::ios::out | std::ios::trunc);
            if (!f_energy) std::cerr << "[ERR] open " << energy_csv << " failed.\n";
            if (f_energy) f_energy << "t, U_box\n";

            // metadata
            std::ofstream meta(det_dir / "metadata.json", std::ios::binary);
            if (meta) {
                meta << "{\n"
                    << "  \"NxT\": " << NxT << ",\n"
                    << "  \"NyT\": " << NyT << ",\n"
                    << "  \"NzT\": " << NzT << ",\n"
                    << "  \"i0\": " << i0 << ", \"i1\": " << i1 << ",\n"
                    << "  \"j0\": " << j0 << ", \"j1\": " << j1 << ",\n"
                    << "  \"k0\": " << k0 << ", \"k1\": " << k1 << ",\n"
                    << "  \"saveEvery\": " << saveEvery << ",\n"
                    << "  \"nSteps\": " << nSteps << ",\n"
                    // << "  \"dx\": " << std::setprecision(18) << dx << ",\n"
                    // << "  \"dy\": " << std::setprecision(18) << dy << ",\n"
                    // << "  \"dz\": " << std::setprecision(18) << dz << ",\n"
                    << "  \"note_spacing\": \"non-uniform grid, see spacing arrays\",\n"
                    << "  \"dt\": " << std::setprecision(18) << dt << ",\n"
                    << "  \"what\": \"total EM energy inside a box (sum over Yee cells)\",\n"
                    << "  \"eps_mu_from\": \"lossless approximation: eps≈dt/bE, mu≈dt/bH\"\n"
                    << "}\n";
            }

            std::cout << "[BoxEnergy] " << det_dir << "\n";
    }

    // One-time recovery of eps/mu from bE/bH (valid when sigma≈0)
    void prepare_coefficients_() {
        if (coeff_ready) return;
        const size_t N = NxT * NyT * NzT;
        epsEx.resize(N); epsEy.resize(N); epsEz.resize(N);
        muHx.resize(N);  muHy.resize(N);  muHz.resize(N);

        auto safe_inv = [](real x)->real { return (std::abs(x) > 1e-30) ? (1.0 / x) : 0.0; };

        for (size_t n = 0; n < N; ++n) {
            epsEx[n] = dt * safe_inv(bEx[n]);   // eps ≈ dt / bE
            epsEy[n] = dt * safe_inv(bEy[n]);
            epsEz[n] = dt * safe_inv(bEz[n]);

            muHx[n] = dt * safe_inv(bHx[n]);    // mu  ≈ dt / bH
            muHy[n] = dt * safe_inv(bHy[n]);
            muHz[n] = dt * safe_inv(bHz[n]);
        }
        coeff_ready = true;
    }

    template<typename Real>
    void try_record(std::size_t n,
        const std::vector<Real>& Ex,
        const std::vector<Real>& Ey,
        const std::vector<Real>& Ez,
        const std::vector<Real>& Hx,
        const std::vector<Real>& Hy,
        const std::vector<Real>& Hz) {
        if (!f_energy) return;
        if (n % saveEvery) return;

        if (!coeff_ready) prepare_coefficients_();

        //const real dV = dx * dy * dz;
        double U = 0.0L;  // Accumulate energy (use double to reduce error)

        for (std::size_t i = i0; i <= i1; ++i)
            for (std::size_t j = j0; j <= j1; ++j)
                for (std::size_t k = k0; k <= k1; ++k) {
                    const std::size_t id = idx3(i, j, k, NyT, NzT);

                    const real dV = dx_array_[i] * dy_array_[j] * dz_array_[k];

                    const double Ex2 = (double)Ex[id] * (double)Ex[id];
                    const double Ey2 = (double)Ey[id] * (double)Ey[id];
                    const double Ez2 = (double)Ez[id] * (double)Ez[id];
                    const double Hx2 = (double)Hx[id] * (double)Hx[id];
                    const double Hy2 = (double)Hy[id] * (double)Hy[id];
                    const double Hz2 = (double)Hz[id] * (double)Hz[id];

                    const double e_part =
                        0.5L * ((double)epsEx[id] * Ex2 +
                            (double)epsEy[id] * Ey2 +
                            (double)epsEz[id] * Ez2);

                    const double h_part =
                        0.5L * ((double)muHx[id] * Hx2 +
                            (double)muHy[id] * Hy2 +
                            (double)muHz[id] * Hz2);

                    U += (e_part + h_part) * (double)dV;
                }

        const real t = n * dt;
        f_energy << std::setprecision(18) << t << "," << (double)U << "\n";
    }
};

#endif // DETECTORS_BOX_ENERGY_HPP
