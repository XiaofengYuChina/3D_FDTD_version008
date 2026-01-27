// mesh_1d_profile_export.hpp - Export 1D mesh grid profiles to CSV files
// Exports dx(x) and dy(y) distributions for the physical domain (excluding PML).

#pragma once

#include <fstream>
#include <string>
#include <iostream>
#include <iomanip>
#include "global_function.hpp"

namespace mesh_export {

inline bool dump_mesh_1d_profiles(const GridSpacing& grid, size_t npml,
                                   const std::string& out_dir = "") {
    std::string prefix = out_dir.empty() ? "" : (out_dir + "/");
    real phys_start_x = grid.x_bounds[npml];
    real phys_start_y = grid.y_bounds[npml];

    // Export dx profile
    {
        std::string filename = prefix + "mesh_dx_profile.csv";
        std::ofstream f(filename);
        if (!f.is_open()) {
            std::cerr << "ERROR: Cannot open " << filename << " for writing\n";
            return false;
        }

        f << std::setprecision(10);
        f << "x_center_nm,dx_nm\n";

        for (size_t i = npml; i < grid.dx.size() - npml; ++i) {
            real x_center_total = 0.5 * (grid.x_bounds[i] + grid.x_bounds[i + 1]);
            real x_center_phys = x_center_total - phys_start_x;
            f << (x_center_phys * 1e9) << "," << (grid.dx[i] * 1e9) << "\n";
        }

        f.close();
        std::cout << "  Exported: " << filename << " ("
                  << (grid.dx.size() - 2 * npml) << " cells)\n";
    }

    // Export dy profile
    {
        std::string filename = prefix + "mesh_dy_profile.csv";
        std::ofstream f(filename);
        if (!f.is_open()) {
            std::cerr << "ERROR: Cannot open " << filename << " for writing\n";
            return false;
        }

        f << std::setprecision(10);
        f << "y_center_nm,dy_nm\n";

        for (size_t j = npml; j < grid.dy.size() - npml; ++j) {
            real y_center_total = 0.5 * (grid.y_bounds[j] + grid.y_bounds[j + 1]);
            real y_center_phys = y_center_total - phys_start_y;
            f << (y_center_phys * 1e9) << "," << (grid.dy[j] * 1e9) << "\n";
        }

        f.close();
        std::cout << "  Exported: " << filename << " ("
                  << (grid.dy.size() - 2 * npml) << " cells)\n";
    }

    return true;
}

} // namespace mesh_export
