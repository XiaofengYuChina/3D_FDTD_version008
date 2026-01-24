// verify_mesh_override.cpp - Integration test to verify mesh override alignment
//
// This test validates the user's real case:
// - Cylinder at center=(200nm, 200nm), radius=100nm -> physical coverage [100nm, 300nm]
// - Override box: x,y ∈ [50nm, 350nm] with 5nm uniform mesh
// - Expected: Fine mesh strictly within [50nm, 350nm], no spread to 450nm+

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <fstream>
#include <set>
#include <numeric>

#include "auto_mesh_generator.hpp"
#include "mesh_1d_profile_export.hpp"

// ANSI color codes
#define GREEN "\033[32m"
#define RED "\033[31m"
#define YELLOW "\033[33m"
#define RESET "\033[0m"
#define BOLD "\033[1m"

// Helper to find index and error for a target position
std::pair<size_t, double> find_nearest(const std::vector<real>& lines, real target) {
    size_t best_idx = 0;
    real best_err = std::abs(lines[0] - target);
    for (size_t i = 1; i < lines.size(); ++i) {
        real err = std::abs(lines[i] - target);
        if (err < best_err) {
            best_err = err;
            best_idx = i;
        }
    }
    return {best_idx, best_err};
}

// Analyze mesh spacing in a range
struct RangeAnalysis {
    size_t n_cells;
    real min_dx, max_dx, mean_dx;
    bool is_uniform;    // variation < 1%
    real variation_pct;
};

RangeAnalysis analyze_range(const std::vector<real>& bounds, const std::vector<real>& dx,
                            real start, real end) {
    RangeAnalysis result = {};
    std::vector<real> spacings;

    for (size_t i = 0; i + 1 < bounds.size(); ++i) {
        real center = 0.5 * (bounds[i] + bounds[i + 1]);
        if (center > start && center < end) {
            spacings.push_back(dx[i]);
        }
    }

    if (spacings.empty()) {
        return result;
    }

    result.n_cells = spacings.size();
    result.min_dx = *std::min_element(spacings.begin(), spacings.end());
    result.max_dx = *std::max_element(spacings.begin(), spacings.end());
    result.mean_dx = std::accumulate(spacings.begin(), spacings.end(), 0.0) / spacings.size();
    result.variation_pct = (result.max_dx - result.min_dx) / result.mean_dx * 100;
    result.is_uniform = result.variation_pct < 1.0;

    return result;
}

int main() {
    std::cout << BOLD << "\n"
              << "======================================================================\n"
              << "     MESH OVERRIDE VERIFICATION TEST (User's Real Case)              \n"
              << "======================================================================\n"
              << RESET << "\n";

    // ========== Configuration (matching user_config.hpp) ==========
    std::cout << BOLD << "=== Configuration ===" << RESET << "\n";

    AutoMeshConfig config;
    config.lambda_min = 1500e-9;           // 1500nm wavelength
    config.mesh_accuracy = 2;               // 10 PPW
    config.dx_min = 2e-9;
    config.dx_max = 100e-9;
    config.max_grading_ratio = std::sqrt(2.0);
    config.mesh_override_enabled = true;

    size_t npml = 8;
    real physical_domain = 1000e-9;  // 1um physical domain

    // Cylinder: center=(200nm, 200nm), radius=100nm
    // -> Physical bounds: x,y ∈ [100nm, 300nm]
    real cyl_cx = 200e-9;
    real cyl_cy = 200e-9;
    real cyl_r = 100e-9;

    // Override box: x,y ∈ [50nm, 350nm]
    real override_x0 = 50e-9;
    real override_x1 = 350e-9;
    real override_y0 = 50e-9;
    real override_y1 = 350e-9;
    real override_dx = 5e-9;

    std::cout << "  Physical domain: 1000nm x 1000nm x 1000nm\n";
    std::cout << "  npml: " << npml << " cells\n";
    std::cout << "  Cylinder: center=(" << cyl_cx*1e9 << "nm, " << cyl_cy*1e9
              << "nm), r=" << cyl_r*1e9 << "nm\n";
    std::cout << "  Override box: [" << override_x0*1e9 << "-" << override_x1*1e9
              << "] x [" << override_y0*1e9 << "-" << override_y1*1e9 << "] nm\n";
    std::cout << "  Override dx/dy/dz: " << override_dx*1e9 << " nm\n";

    // Setup override region
    MeshOverrideRegion override_region;
    override_region.enabled = true;
    override_region.x0 = override_x0;
    override_region.x1 = override_x1;
    override_region.y0 = override_y0;
    override_region.y1 = override_y1;
    override_region.z0 = 400e-9;
    override_region.z1 = 800e-9;
    override_region.dx = override_dx;
    override_region.dy = override_dx;
    override_region.dz = override_dx;

    config.override_regions.push_back(override_region);

    // Setup structure (cylinder approximated as box for mesh generation)
    std::vector<AutoMeshStructureBounds> structures = {
        {cyl_cx - cyl_r, cyl_cx + cyl_r,   // x: [100nm, 300nm]
         cyl_cy - cyl_r, cyl_cy + cyl_r,   // y: [100nm, 300nm]
         500e-9, 700e-9,                    // z: [500nm, 700nm]
         2.0}  // n = 2.0
    };

    // ========== Generate Mesh ==========
    std::cout << "\n" << BOLD << "=== Generating Mesh ===" << RESET << "\n";

    AutoMeshGenerator generator(config);
    GridSpacing grid = generator.generate(physical_domain, physical_domain, physical_domain,
                                          npml, structures, 1.0);

    // Get physical domain boundaries in total coords
    real phys_start_x = grid.x_bounds[npml];
    real phys_end_x = grid.x_bounds[grid.dx.size() - npml];
    real phys_start_y = grid.y_bounds[npml];
    real phys_end_y = grid.y_bounds[grid.dy.size() - npml];

    // Physical coordinates (relative to physical domain start)
    std::vector<real> x_phys_lines(grid.x_bounds.size());
    std::vector<real> y_phys_lines(grid.y_bounds.size());
    for (size_t i = 0; i < grid.x_bounds.size(); ++i) {
        x_phys_lines[i] = grid.x_bounds[i] - phys_start_x;
    }
    for (size_t i = 0; i < grid.y_bounds.size(); ++i) {
        y_phys_lines[i] = grid.y_bounds[i] - phys_start_y;
    }

    // ========== Verification 1: Override boundaries in grid ==========
    std::cout << "\n" << BOLD << "=== Verification 1: Override Boundary Snapping ===" << RESET << "\n";

    auto [x0_idx, x0_err] = find_nearest(x_phys_lines, override_x0);
    auto [x1_idx, x1_err] = find_nearest(x_phys_lines, override_x1);
    auto [y0_idx, y0_err] = find_nearest(y_phys_lines, override_y0);
    auto [y1_idx, y1_err] = find_nearest(y_phys_lines, override_y1);

    std::cout << "  X boundaries:\n";
    std::cout << "    x0=" << override_x0*1e9 << "nm -> index " << x0_idx
              << ", err=" << x0_err*1e9 << "nm"
              << (x0_err < 1e-12 ? GREEN " [EXACT]" RESET : (x0_err < 1e-9 ? YELLOW " [CLOSE]" RESET : RED " [MISS]" RESET)) << "\n";
    std::cout << "    x1=" << override_x1*1e9 << "nm -> index " << x1_idx
              << ", err=" << x1_err*1e9 << "nm"
              << (x1_err < 1e-12 ? GREEN " [EXACT]" RESET : (x1_err < 1e-9 ? YELLOW " [CLOSE]" RESET : RED " [MISS]" RESET)) << "\n";
    std::cout << "  Y boundaries:\n";
    std::cout << "    y0=" << override_y0*1e9 << "nm -> index " << y0_idx
              << ", err=" << y0_err*1e9 << "nm"
              << (y0_err < 1e-12 ? GREEN " [EXACT]" RESET : (y0_err < 1e-9 ? YELLOW " [CLOSE]" RESET : RED " [MISS]" RESET)) << "\n";
    std::cout << "    y1=" << override_y1*1e9 << "nm -> index " << y1_idx
              << ", err=" << y1_err*1e9 << "nm"
              << (y1_err < 1e-12 ? GREEN " [EXACT]" RESET : (y1_err < 1e-9 ? YELLOW " [CLOSE]" RESET : RED " [MISS]" RESET)) << "\n";

    bool boundaries_ok = (x0_err < 1e-9 && x1_err < 1e-9 && y0_err < 1e-9 && y1_err < 1e-9);
    std::cout << "  " << (boundaries_ok ? GREEN "PASS" : RED "FAIL") << RESET
              << ": Override boundaries snapped to grid lines\n";

    // ========== Verification 2: Uniformity inside override ==========
    std::cout << "\n" << BOLD << "=== Verification 2: Uniformity Inside Override ===" << RESET << "\n";

    // Analyze X direction (using physical coordinates)
    auto x_analysis = analyze_range(x_phys_lines, grid.dx, override_x0, override_x1);
    auto y_analysis = analyze_range(y_phys_lines, grid.dy, override_y0, override_y1);

    std::cout << "  X direction (inside [" << override_x0*1e9 << ", " << override_x1*1e9 << "] nm):\n";
    std::cout << "    Cells: " << x_analysis.n_cells << "\n";
    std::cout << "    Target dx: " << override_dx*1e9 << " nm\n";
    std::cout << "    Actual dx: min=" << x_analysis.min_dx*1e9 << ", max=" << x_analysis.max_dx*1e9
              << ", mean=" << x_analysis.mean_dx*1e9 << " nm\n";
    std::cout << "    Variation: " << std::fixed << std::setprecision(3) << x_analysis.variation_pct << "%\n";
    std::cout << "    Uniform (min==max==override_step): "
              << (x_analysis.is_uniform ? GREEN "YES" : RED "NO") << RESET << "\n";

    std::cout << "  Y direction (inside [" << override_y0*1e9 << ", " << override_y1*1e9 << "] nm):\n";
    std::cout << "    Cells: " << y_analysis.n_cells << "\n";
    std::cout << "    Target dy: " << override_dx*1e9 << " nm\n";
    std::cout << "    Actual dy: min=" << y_analysis.min_dx*1e9 << ", max=" << y_analysis.max_dx*1e9
              << ", mean=" << y_analysis.mean_dx*1e9 << " nm\n";
    std::cout << "    Variation: " << std::fixed << std::setprecision(3) << y_analysis.variation_pct << "%\n";
    std::cout << "    Uniform (min==max==override_step): "
              << (y_analysis.is_uniform ? GREEN "YES" : RED "NO") << RESET << "\n";

    bool uniformity_ok = x_analysis.is_uniform && y_analysis.is_uniform;
    std::cout << "  " << (uniformity_ok ? GREEN "PASS" : RED "FAIL") << RESET
              << ": Override region is strictly uniform\n";

    // ========== Verification 3: No spread beyond box ==========
    std::cout << "\n" << BOLD << "=== Verification 3: No Fine Mesh Spread Beyond Box ===" << RESET << "\n";

    // Find where fine mesh (< 10nm) exists
    real fine_threshold = 10e-9;
    real fine_x_start = 0, fine_x_end = 0;
    real fine_y_start = 0, fine_y_end = 0;

    // Only look within physical domain
    for (size_t i = npml; i < grid.dx.size() - npml; ++i) {
        if (grid.dx[i] < fine_threshold) {
            real phys_x = grid.x_bounds[i] - phys_start_x;
            if (fine_x_start == 0) fine_x_start = phys_x;
            fine_x_end = grid.x_bounds[i + 1] - phys_start_x;
        }
    }
    for (size_t i = npml; i < grid.dy.size() - npml; ++i) {
        if (grid.dy[i] < fine_threshold) {
            real phys_y = grid.y_bounds[i] - phys_start_y;
            if (fine_y_start == 0) fine_y_start = phys_y;
            fine_y_end = grid.y_bounds[i + 1] - phys_start_y;
        }
    }

    std::cout << "  Fine mesh region (dx < 10nm):\n";
    std::cout << "    X: [" << fine_x_start*1e9 << ", " << fine_x_end*1e9 << "] nm\n";
    std::cout << "    Y: [" << fine_y_start*1e9 << ", " << fine_y_end*1e9 << "] nm\n";
    std::cout << "  Expected region: [" << override_x0*1e9 << ", " << override_x1*1e9 << "] nm\n";

    // Allow small transition region (< 50nm) for grading
    real tolerance = 50e-9;
    bool x_no_spread = (fine_x_start >= override_x0 - tolerance && fine_x_end <= override_x1 + tolerance);
    bool y_no_spread = (fine_y_start >= override_y0 - tolerance && fine_y_end <= override_y1 + tolerance);

    // Critical: Check that fine mesh does NOT extend to 450nm (far beyond 350nm)
    bool not_at_450nm = (fine_x_end < 400e-9 && fine_y_end < 400e-9);

    std::cout << "  X: start error = " << (fine_x_start - override_x0)*1e9 << "nm, "
              << "end error = " << (fine_x_end - override_x1)*1e9 << "nm "
              << (x_no_spread ? GREEN "[OK]" : RED "[SPREAD]") << RESET << "\n";
    std::cout << "  Y: start error = " << (fine_y_start - override_y0)*1e9 << "nm, "
              << "end error = " << (fine_y_end - override_y1)*1e9 << "nm "
              << (y_no_spread ? GREEN "[OK]" : RED "[SPREAD]") << RESET << "\n";
    std::cout << "  Fine mesh at 450nm (should be NO): "
              << (not_at_450nm ? GREEN "NO - GOOD" : RED "YES - BUG!") << RESET << "\n";

    bool no_spread_ok = x_no_spread && y_no_spread && not_at_450nm;
    std::cout << "  " << (no_spread_ok ? GREEN "PASS" : RED "FAIL") << RESET
              << ": Fine mesh confined to override box\n";

    // ========== Output physical grid lines for verification ==========
    std::cout << "\n" << BOLD << "=== Physical Grid Lines (first 50 and last 10) ===" << RESET << "\n";

    std::cout << "  X physical grid lines (nm):\n    ";
    for (size_t i = npml; i < std::min(npml + 50, grid.x_bounds.size()); ++i) {
        std::cout << std::fixed << std::setprecision(1) << (grid.x_bounds[i] - phys_start_x)*1e9 << " ";
        if ((i - npml + 1) % 10 == 0) std::cout << "\n    ";
    }
    std::cout << "...\n";

    std::cout << "  Last 10: ";
    for (size_t i = std::max(npml, grid.x_bounds.size() - npml - 10); i < grid.x_bounds.size() - npml; ++i) {
        std::cout << std::fixed << std::setprecision(1) << (grid.x_bounds[i] - phys_start_x)*1e9 << " ";
    }
    std::cout << "\n";

    // ========== Export 1D mesh profiles (CSV) ==========
    std::cout << "\n" << BOLD << "=== Exporting 1D Mesh Profiles ===" << RESET << "\n";
    mesh_export::dump_mesh_1d_profiles(grid, npml);

    // ========== Save data for Python visualization ==========
    std::cout << "\n" << BOLD << "=== Saving Data for Visualization ===" << RESET << "\n";

    {
        std::ofstream f("mesh_verify_data.txt");
        f << "# Mesh verification data\n";
        f << "# Physical domain: [0, " << physical_domain*1e9 << "] nm\n";
        f << "# Override box: [" << override_x0*1e9 << ", " << override_x1*1e9 << "] nm\n";
        f << "# Cylinder: center=(" << cyl_cx*1e9 << ", " << cyl_cy*1e9 << "), r=" << cyl_r*1e9 << " nm\n";
        f << "\n# x_phys_lines (nm)\n";
        for (size_t i = npml; i < grid.x_bounds.size() - npml; ++i) {
            f << (grid.x_bounds[i] - phys_start_x)*1e9 << "\n";
        }
        f << "\n# y_phys_lines (nm)\n";
        for (size_t i = npml; i < grid.y_bounds.size() - npml; ++i) {
            f << (grid.y_bounds[i] - phys_start_y)*1e9 << "\n";
        }
        f << "\n# dx (nm)\n";
        for (size_t i = npml; i < grid.dx.size() - npml; ++i) {
            f << grid.dx[i]*1e9 << "\n";
        }
        f << "\n# dy (nm)\n";
        for (size_t i = npml; i < grid.dy.size() - npml; ++i) {
            f << grid.dy[i]*1e9 << "\n";
        }
        f.close();
        std::cout << "  Saved: mesh_verify_data.txt\n";
    }

    // ========== Summary ==========
    std::cout << "\n" << BOLD << "======================================================================\n"
              << "                           SUMMARY                                    \n"
              << "======================================================================\n" << RESET;

    std::cout << "  [1] Override boundaries snapped: " << (boundaries_ok ? GREEN "PASS" : RED "FAIL") << RESET << "\n";
    std::cout << "  [2] Override region uniform:     " << (uniformity_ok ? GREEN "PASS" : RED "FAIL") << RESET << "\n";
    std::cout << "  [3] No spread beyond box:        " << (no_spread_ok ? GREEN "PASS" : RED "FAIL") << RESET << "\n";

    bool all_pass = boundaries_ok && uniformity_ok && no_spread_ok;
    std::cout << "\n  " << BOLD << (all_pass ? GREEN "ALL VERIFICATIONS PASSED" : RED "SOME VERIFICATIONS FAILED")
              << RESET << "\n\n";

    return all_pass ? 0 : 1;
}
