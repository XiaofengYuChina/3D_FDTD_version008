// test_auto_mesh.cpp - Unit tests for the Auto Non-uniform Mesh Generator
//
// This test program validates:
// 1. Two-material interface snapping and mesh refinement
// 2. Single override box uniformity and boundary snapping
// 3. Multiple overlapping overrides (finer mesh wins)
// 4. Global grading ratio constraint (sqrt(2) by default)
// 5. PPW-based cell size calculation
// 6. Interface positions appearing as grid lines

#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include <iomanip>
#include <string>
#include <algorithm>
#include <set>
#include <numeric>

#include "auto_mesh_generator.hpp"

// ANSI color codes for output
#define GREEN "\033[32m"
#define RED "\033[31m"
#define YELLOW "\033[33m"
#define RESET "\033[0m"
#define BOLD "\033[1m"

struct TestResult {
    std::string name;
    bool passed;
    std::string message;
};

std::vector<TestResult> all_results;

void report_test(const std::string& name, bool passed, const std::string& msg = "") {
    all_results.push_back({name, passed, msg});
    if (passed) {
        std::cout << GREEN << "[PASS] " << RESET << name;
    } else {
        std::cout << RED << "[FAIL] " << RESET << name;
    }
    if (!msg.empty()) {
        std::cout << " - " << msg;
    }
    std::cout << "\n";
}

// Helper: Check if a position is present in grid lines (within tolerance)
bool position_in_grid(const std::vector<real>& bounds, real pos, real tol = 1e-12) {
    for (real b : bounds) {
        if (std::abs(b - pos) < tol) return true;
    }
    return false;
}

// Helper: Find spacing at a position
real spacing_at_position(const std::vector<real>& bounds, const std::vector<real>& dx, real pos) {
    for (size_t i = 0; i + 1 < bounds.size(); ++i) {
        if (pos >= bounds[i] && pos < bounds[i + 1]) {
            return dx[i];
        }
    }
    return dx.back();
}

// Helper: Get maximum grading ratio in a spacing array
real max_grading_ratio(const std::vector<real>& dx) {
    real max_r = 1.0;
    for (size_t i = 1; i < dx.size(); ++i) {
        real r = std::max(dx[i] / dx[i - 1], dx[i - 1] / dx[i]);
        max_r = std::max(max_r, r);
    }
    return max_r;
}

// ==================== Test 1: Two-Material Interface Case ====================
bool test_two_material_interface() {
    std::cout << "\n" << BOLD << "=== Test 1: Two-Material Interface Case ===" << RESET << "\n";

    AutoMeshConfig config;
    config.lambda_min = 500e-9;          // 500 nm wavelength
    config.mesh_accuracy = 2;             // 10 ppw
    config.dx_min = 5e-9;
    config.dx_max = 100e-9;
    config.max_grading_ratio = std::sqrt(2.0);  // sqrt(2)

    AutoMeshGenerator generator(config);

    // Two-material case: air (n=1) and glass (n=1.5)
    // Structure from 400nm to 600nm (in physical domain)
    real struct_start = 400e-9;
    real struct_end = 600e-9;
    real struct_n = 1.5;

    std::vector<AutoMeshStructureBounds> structures = {
        {struct_start, struct_end, struct_start, struct_end, struct_start, struct_end, struct_n}
    };

    GridSpacing grid = generator.generate(1000e-9, 1000e-9, 1000e-9, 0, structures, 1.0);

    // Test 1a: Interface snapping - check that interface positions are grid lines
    bool start_snapped = position_in_grid(grid.x_bounds, struct_start);
    bool end_snapped = position_in_grid(grid.x_bounds, struct_end);

    report_test("Interface start snapped to grid", start_snapped,
                "Position " + std::to_string(struct_start * 1e9) + " nm");
    report_test("Interface end snapped to grid", end_snapped,
                "Position " + std::to_string(struct_end * 1e9) + " nm");

    // Test 1b: High-n region gets finer mesh
    // Expected: in glass (n=1.5), dx ~ lambda_min / (n * ppw) = 500nm / (1.5 * 10) = 33.3nm
    // In air: dx ~ lambda_min / ppw = 500nm / 10 = 50nm
    real dx_in_glass = spacing_at_position(grid.x_bounds, grid.dx, 500e-9);
    real dx_in_air = spacing_at_position(grid.x_bounds, grid.dx, 100e-9);

    std::cout << "  dx in glass (n=1.5): " << dx_in_glass * 1e9 << " nm\n";
    std::cout << "  dx in air (n=1): " << dx_in_air * 1e9 << " nm\n";
    std::cout << "  Ratio (air/glass): " << dx_in_air / dx_in_glass << "\n";

    // Glass should have finer or equal mesh compared to air
    bool glass_finer = dx_in_glass <= dx_in_air * 1.2;
    report_test("High-n region has finer mesh", glass_finer,
                "Glass dx=" + std::to_string(dx_in_glass * 1e9) + "nm, Air dx=" + std::to_string(dx_in_air * 1e9) + "nm");

    // Test 1c: Grading ratio constraint
    real max_ratio = max_grading_ratio(grid.dx);
    bool ratio_ok = max_ratio <= config.max_grading_ratio * 1.05;  // 5% tolerance
    report_test("Grading ratio <= sqrt(2)", ratio_ok,
                "Max ratio: " + std::to_string(max_ratio) + " (target: " + std::to_string(config.max_grading_ratio) + ")");

    return start_snapped && end_snapped && glass_finer && ratio_ok;
}

// ==================== Test 2: Single Override Box ====================
bool test_single_override_box() {
    std::cout << "\n" << BOLD << "=== Test 2: Single Override Box ===" << RESET << "\n";

    AutoMeshConfig config;
    config.lambda_min = 500e-9;
    config.mesh_accuracy = 2;
    config.dx_min = 2e-9;
    config.dx_max = 100e-9;
    config.max_grading_ratio = std::sqrt(2.0);
    config.mesh_override_enabled = true;

    // Override region: [300nm-700nm] with 10nm uniform mesh
    MeshOverrideRegion override_region;
    override_region.enabled = true;
    override_region.x0 = 300e-9;
    override_region.x1 = 700e-9;
    override_region.y0 = 300e-9;
    override_region.y1 = 700e-9;
    override_region.z0 = 300e-9;
    override_region.z1 = 700e-9;
    override_region.dx = 10e-9;
    override_region.dy = 10e-9;
    override_region.dz = 10e-9;

    config.override_regions.push_back(override_region);

    AutoMeshGenerator generator(config);

    std::vector<AutoMeshStructureBounds> structures;  // No structures, just override
    GridSpacing grid = generator.generate(1000e-9, 1000e-9, 1000e-9, 0, structures, 1.0);

    // Test 2a: Override boundary snapping
    bool start_snapped = position_in_grid(grid.x_bounds, override_region.x0);
    bool end_snapped = position_in_grid(grid.x_bounds, override_region.x1);

    report_test("Override start boundary snapped", start_snapped,
                "Position " + std::to_string(override_region.x0 * 1e9) + " nm");
    report_test("Override end boundary snapped", end_snapped,
                "Position " + std::to_string(override_region.x1 * 1e9) + " nm");

    // Test 2b: Inside override is uniform
    std::vector<real> override_spacings;
    for (size_t i = 0; i + 1 < grid.x_bounds.size(); ++i) {
        real center = 0.5 * (grid.x_bounds[i] + grid.x_bounds[i + 1]);
        if (center > override_region.x0 && center < override_region.x1) {
            override_spacings.push_back(grid.dx[i]);
        }
    }

    if (!override_spacings.empty()) {
        real min_sp = *std::min_element(override_spacings.begin(), override_spacings.end());
        real max_sp = *std::max_element(override_spacings.begin(), override_spacings.end());
        real mean_sp = std::accumulate(override_spacings.begin(), override_spacings.end(), 0.0) / override_spacings.size();
        real variation = (max_sp - min_sp) / mean_sp;

        std::cout << "  Inside override: " << override_spacings.size() << " cells\n";
        std::cout << "  Target dx: " << override_region.dx * 1e9 << " nm\n";
        std::cout << "  Actual dx: " << min_sp * 1e9 << " - " << max_sp * 1e9 << " nm (mean: " << mean_sp * 1e9 << " nm)\n";
        std::cout << "  Variation: " << variation * 100 << "%\n";

        // Allow 5% variation for uniformity check
        bool uniform = variation < 0.05;
        report_test("Override region is uniform", uniform,
                    "Variation: " + std::to_string(variation * 100) + "% (target: <5%)");

        // Check that actual spacing is close to target
        bool close_to_target = std::abs(mean_sp - override_region.dx) / override_region.dx < 0.15;
        report_test("Override spacing close to target", close_to_target,
                    "Mean: " + std::to_string(mean_sp * 1e9) + "nm, Target: " + std::to_string(override_region.dx * 1e9) + "nm");
    }

    // Test 2c: Override takes priority over auto
    // Outside override (at 100nm), spacing should be auto-calculated
    real dx_outside = spacing_at_position(grid.x_bounds, grid.dx, 100e-9);
    real dx_inside = spacing_at_position(grid.x_bounds, grid.dx, 500e-9);

    std::cout << "  dx outside override (100nm): " << dx_outside * 1e9 << " nm\n";
    std::cout << "  dx inside override (500nm): " << dx_inside * 1e9 << " nm\n";

    // Inside should be close to override target (10nm)
    bool override_applied = std::abs(dx_inside - override_region.dx) / override_region.dx < 0.2;
    report_test("Override takes priority", override_applied);

    // Test 2d: Smooth transition at boundary
    real max_ratio = max_grading_ratio(grid.dx);
    bool ratio_ok = max_ratio <= config.max_grading_ratio * 1.1;
    report_test("Smooth transition (ratio <= sqrt(2))", ratio_ok,
                "Max ratio: " + std::to_string(max_ratio));

    return start_snapped && end_snapped;
}

// ==================== Test 3: Multiple Overlapping Overrides ====================
bool test_multiple_overrides() {
    std::cout << "\n" << BOLD << "=== Test 3: Multiple Overlapping Overrides ===" << RESET << "\n";

    AutoMeshConfig config;
    config.lambda_min = 500e-9;
    config.mesh_accuracy = 2;
    config.dx_min = 2e-9;
    config.dx_max = 100e-9;
    config.max_grading_ratio = std::sqrt(2.0);
    config.mesh_override_enabled = true;

    // Override 1: Large region [200nm-800nm] with 20nm mesh
    MeshOverrideRegion override1;
    override1.enabled = true;
    override1.x0 = 200e-9;
    override1.x1 = 800e-9;
    override1.y0 = 200e-9;
    override1.y1 = 800e-9;
    override1.z0 = 200e-9;
    override1.z1 = 800e-9;
    override1.dx = 20e-9;
    override1.dy = 20e-9;
    override1.dz = 20e-9;

    // Override 2: Smaller region [400nm-600nm] with 5nm mesh (finer)
    MeshOverrideRegion override2;
    override2.enabled = true;
    override2.x0 = 400e-9;
    override2.x1 = 600e-9;
    override2.y0 = 400e-9;
    override2.y1 = 600e-9;
    override2.z0 = 400e-9;
    override2.z1 = 600e-9;
    override2.dx = 5e-9;
    override2.dy = 5e-9;
    override2.dz = 5e-9;

    config.override_regions.push_back(override1);
    config.override_regions.push_back(override2);

    AutoMeshGenerator generator(config);

    std::vector<AutoMeshStructureBounds> structures;
    GridSpacing grid = generator.generate(1000e-9, 1000e-9, 1000e-9, 0, structures, 1.0);

    // Test 3a: Check spacing in overlap region (should use finer 5nm)
    real dx_in_overlap = spacing_at_position(grid.x_bounds, grid.dx, 500e-9);
    std::cout << "  dx in overlap region (500nm): " << dx_in_overlap * 1e9 << " nm\n";
    std::cout << "  Expected (finer wins): ~5 nm\n";

    bool finer_wins = dx_in_overlap < 15e-9;  // Should be closer to 5nm than 20nm
    report_test("Finer mesh wins in overlap", finer_wins,
                "Actual: " + std::to_string(dx_in_overlap * 1e9) + "nm, Expected: ~5nm");

    // Test 3b: Check spacing in override1-only region (should use 20nm)
    real dx_in_override1_only = spacing_at_position(grid.x_bounds, grid.dx, 300e-9);
    std::cout << "  dx in override1-only region (300nm): " << dx_in_override1_only * 1e9 << " nm\n";
    std::cout << "  Expected: ~20 nm\n";

    bool override1_applied = std::abs(dx_in_override1_only - override1.dx) / override1.dx < 0.3;
    report_test("Override1 applied outside overlap", override1_applied,
                "Actual: " + std::to_string(dx_in_override1_only * 1e9) + "nm, Expected: ~20nm");

    // Test 3c: All override boundaries should be snapped
    bool all_boundaries_snapped = true;
    if (!position_in_grid(grid.x_bounds, override1.x0)) all_boundaries_snapped = false;
    if (!position_in_grid(grid.x_bounds, override1.x1)) all_boundaries_snapped = false;
    if (!position_in_grid(grid.x_bounds, override2.x0)) all_boundaries_snapped = false;
    if (!position_in_grid(grid.x_bounds, override2.x1)) all_boundaries_snapped = false;

    report_test("All override boundaries snapped", all_boundaries_snapped);

    return finer_wins && all_boundaries_snapped;
}

// ==================== Test 4: Global Grading Ratio Check ====================
bool test_global_grading_ratio() {
    std::cout << "\n" << BOLD << "=== Test 4: Global Grading Ratio Check ===" << RESET << "\n";

    AutoMeshConfig config;
    config.lambda_min = 500e-9;
    config.mesh_accuracy = 3;  // Higher accuracy for more cells
    config.dx_min = 3e-9;
    config.dx_max = 80e-9;
    config.max_grading_ratio = std::sqrt(2.0);

    AutoMeshGenerator generator(config);

    // Complex structure with high contrast
    std::vector<AutoMeshStructureBounds> structures = {
        {200e-9, 300e-9, 200e-9, 300e-9, 200e-9, 300e-9, 3.5},  // High-n region
        {600e-9, 800e-9, 600e-9, 800e-9, 600e-9, 800e-9, 2.0}   // Medium-n region
    };

    GridSpacing grid = generator.generate(1000e-9, 1000e-9, 1000e-9, 0, structures, 1.0);

    // Check grading ratio for all axes
    real max_ratio_x = max_grading_ratio(grid.dx);
    real max_ratio_y = max_grading_ratio(grid.dy);
    real max_ratio_z = max_grading_ratio(grid.dz);

    std::cout << "  Max grading ratio X: " << max_ratio_x << "\n";
    std::cout << "  Max grading ratio Y: " << max_ratio_y << "\n";
    std::cout << "  Max grading ratio Z: " << max_ratio_z << "\n";
    std::cout << "  Target: <= " << config.max_grading_ratio << " (sqrt(2))\n";

    // Allow 10% tolerance due to cell count rounding
    real tolerance = 1.1;
    bool x_ok = max_ratio_x <= config.max_grading_ratio * tolerance;
    bool y_ok = max_ratio_y <= config.max_grading_ratio * tolerance;
    bool z_ok = max_ratio_z <= config.max_grading_ratio * tolerance;

    report_test("X grading ratio constraint", x_ok,
                std::to_string(max_ratio_x) + " <= " + std::to_string(config.max_grading_ratio * tolerance));
    report_test("Y grading ratio constraint", y_ok,
                std::to_string(max_ratio_y) + " <= " + std::to_string(config.max_grading_ratio * tolerance));
    report_test("Z grading ratio constraint", z_ok,
                std::to_string(max_ratio_z) + " <= " + std::to_string(config.max_grading_ratio * tolerance));

    // Verify grid is strictly increasing
    bool x_monotonic = true, y_monotonic = true, z_monotonic = true;
    for (size_t i = 1; i < grid.x_bounds.size(); ++i) {
        if (grid.x_bounds[i] <= grid.x_bounds[i - 1]) x_monotonic = false;
    }
    for (size_t i = 1; i < grid.y_bounds.size(); ++i) {
        if (grid.y_bounds[i] <= grid.y_bounds[i - 1]) y_monotonic = false;
    }
    for (size_t i = 1; i < grid.z_bounds.size(); ++i) {
        if (grid.z_bounds[i] <= grid.z_bounds[i - 1]) z_monotonic = false;
    }

    report_test("X grid strictly increasing", x_monotonic);
    report_test("Y grid strictly increasing", y_monotonic);
    report_test("Z grid strictly increasing", z_monotonic);

    return x_ok && y_ok && z_ok && x_monotonic && y_monotonic && z_monotonic;
}

// ==================== Test 5: PPW Mapping Verification ====================
bool test_ppw_mapping() {
    std::cout << "\n" << BOLD << "=== Test 5: PPW Mapping Verification ===" << RESET << "\n";

    // Verify PPW mapping: acc=1->6, acc=2->10, acc=3->14, acc>=3: 14+(acc-3)*4
    struct PPWTestCase {
        int accuracy;
        real expected_ppw;
    };

    std::vector<PPWTestCase> test_cases = {
        {1, 6.0},
        {2, 10.0},
        {3, 14.0},
        {4, 18.0},
        {5, 22.0},
        {8, 34.0}
    };

    bool all_correct = true;
    for (const auto& tc : test_cases) {
        AutoMeshConfig config;
        config.mesh_accuracy = tc.accuracy;
        real actual_ppw = config.get_ppw_target();

        bool correct = std::abs(actual_ppw - tc.expected_ppw) < 0.01;
        all_correct &= correct;

        std::cout << "  Accuracy " << tc.accuracy << ": expected " << tc.expected_ppw
                  << " ppw, got " << actual_ppw << " ppw"
                  << (correct ? " [OK]" : " [FAIL]") << "\n";
    }

    report_test("PPW mapping correct", all_correct);

    // Test that PPW affects mesh density
    AutoMeshConfig config_low;
    config_low.lambda_min = 500e-9;
    config_low.mesh_accuracy = 1;  // 6 ppw
    config_low.dx_min = 2e-9;
    config_low.dx_max = 200e-9;

    AutoMeshConfig config_high;
    config_high.lambda_min = 500e-9;
    config_high.mesh_accuracy = 3;  // 14 ppw
    config_high.dx_min = 2e-9;
    config_high.dx_max = 200e-9;

    AutoMeshGenerator gen_low(config_low);
    AutoMeshGenerator gen_high(config_high);

    std::vector<AutoMeshStructureBounds> structures;
    GridSpacing grid_low = gen_low.generate(1000e-9, 1000e-9, 1000e-9, 0, structures, 1.0);
    GridSpacing grid_high = gen_high.generate(1000e-9, 1000e-9, 1000e-9, 0, structures, 1.0);

    std::cout << "\n  Low accuracy (acc=1): " << grid_low.dx.size() << " cells\n";
    std::cout << "  High accuracy (acc=3): " << grid_high.dx.size() << " cells\n";

    bool higher_accuracy_more_cells = grid_high.dx.size() > grid_low.dx.size();
    report_test("Higher accuracy produces finer mesh", higher_accuracy_more_cells,
                "Low: " + std::to_string(grid_low.dx.size()) + ", High: " + std::to_string(grid_high.dx.size()));

    return all_correct && higher_accuracy_more_cells;
}

// ==================== Test 6: Demo Output ====================
bool test_demo_output() {
    std::cout << "\n" << BOLD << "=== Test 6: Demo Output ===" << RESET << "\n";

    AutoMeshConfig config;
    config.lambda_min = 500e-9;
    config.mesh_accuracy = 2;
    config.dx_min = 5e-9;
    config.dx_max = 80e-9;
    config.max_grading_ratio = std::sqrt(2.0);
    config.mesh_override_enabled = true;

    // Add an override region
    MeshOverrideRegion override_region;
    override_region.enabled = true;
    override_region.x0 = 350e-9;
    override_region.x1 = 650e-9;
    override_region.y0 = 350e-9;
    override_region.y1 = 650e-9;
    override_region.z0 = 350e-9;
    override_region.z1 = 650e-9;
    override_region.dx = 8e-9;
    override_region.dy = 8e-9;
    override_region.dz = 8e-9;
    config.override_regions.push_back(override_region);

    AutoMeshGenerator generator(config);

    // Structure with different refractive index
    std::vector<AutoMeshStructureBounds> structures = {
        {400e-9, 600e-9, 400e-9, 600e-9, 400e-9, 600e-9, 2.5}
    };

    GridSpacing grid = generator.generate(1000e-9, 1000e-9, 1000e-9, 0, structures, 1.0);

    // Print demo statistics
    std::cout << "\n" << BOLD << "  Demo Case Statistics:" << RESET << "\n";
    std::cout << "  =====================\n";
    std::cout << "  Configuration:\n";
    std::cout << "    Lambda_min: " << config.lambda_min * 1e9 << " nm\n";
    std::cout << "    Mesh accuracy: " << config.mesh_accuracy << " (PPW: " << config.get_ppw_target() << ")\n";
    std::cout << "    Max grading ratio: " << config.max_grading_ratio << " (sqrt(2))\n";
    std::cout << "    Override region: [" << override_region.x0 * 1e9 << "-" << override_region.x1 * 1e9 << "] nm, dx=" << override_region.dx * 1e9 << " nm\n";
    std::cout << "    Structure (n=2.5): [400-600] nm\n";

    std::cout << "\n  Results:\n";
    std::cout << "    Total cell counts: Nx=" << grid.dx.size()
              << ", Ny=" << grid.dy.size()
              << ", Nz=" << grid.dz.size() << "\n";

    real dx_min = *std::min_element(grid.dx.begin(), grid.dx.end());
    real dx_max = *std::max_element(grid.dx.begin(), grid.dx.end());
    real dy_min = *std::min_element(grid.dy.begin(), grid.dy.end());
    real dy_max = *std::max_element(grid.dy.begin(), grid.dy.end());
    real dz_min = *std::min_element(grid.dz.begin(), grid.dz.end());
    real dz_max = *std::max_element(grid.dz.begin(), grid.dz.end());

    std::cout << "    dx range: " << dx_min * 1e9 << " - " << dx_max * 1e9 << " nm\n";
    std::cout << "    dy range: " << dy_min * 1e9 << " - " << dy_max * 1e9 << " nm\n";
    std::cout << "    dz range: " << dz_min * 1e9 << " - " << dz_max * 1e9 << " nm\n";

    real max_ratio_x = max_grading_ratio(grid.dx);
    real max_ratio_y = max_grading_ratio(grid.dy);
    real max_ratio_z = max_grading_ratio(grid.dz);

    std::cout << "    Max grading ratio (X/Y/Z): " << max_ratio_x << " / " << max_ratio_y << " / " << max_ratio_z << "\n";

    // Check uniformity inside override
    std::vector<real> override_dx;
    for (size_t i = 0; i + 1 < grid.x_bounds.size(); ++i) {
        real center = 0.5 * (grid.x_bounds[i] + grid.x_bounds[i + 1]);
        if (center > override_region.x0 && center < override_region.x1) {
            override_dx.push_back(grid.dx[i]);
        }
    }

    if (!override_dx.empty()) {
        real mean_dx = std::accumulate(override_dx.begin(), override_dx.end(), 0.0) / override_dx.size();
        real min_dx_ov = *std::min_element(override_dx.begin(), override_dx.end());
        real max_dx_ov = *std::max_element(override_dx.begin(), override_dx.end());
        real variation = (max_dx_ov - min_dx_ov) / mean_dx * 100;

        std::cout << "\n  Override Region Check:\n";
        std::cout << "    Cells inside override: " << override_dx.size() << "\n";
        std::cout << "    Target spacing: " << override_region.dx * 1e9 << " nm\n";
        std::cout << "    Actual mean spacing: " << mean_dx * 1e9 << " nm\n";
        std::cout << "    Spacing variation: " << variation << "%\n";
    }

    report_test("Demo output generated", true);
    return true;
}

// ==================== Test 7: Two-Step Mesh Generation and Override Alignment ====================
// This test verifies the two-step mesh generation design:
// 1. Physical domain mesh generated first with override regions in physical coordinates
// 2. PML extended based on boundary cell sizes
bool test_pml_override_alignment() {
    std::cout << "\n" << BOLD << "=== Test 7: Two-Step Mesh Generation and Override Alignment ===" << RESET << "\n";
    std::cout << "  (Regression test for mesh override coordinate mapping bug)\n\n";

    // Configuration similar to user_config.hpp
    AutoMeshConfig config;
    config.lambda_min = 1500e-9;          // 1500 nm wavelength
    config.mesh_accuracy = 2;             // 10 ppw
    config.dx_min = 2e-9;
    config.dx_max = 100e-9;
    config.max_grading_ratio = std::sqrt(2.0);
    config.mesh_override_enabled = true;

    // PML configuration
    size_t npml = 8;

    std::cout << "  Config: npml=" << npml << ", dx_max=" << config.dx_max * 1e9 << "nm\n";

    // User-defined override region in PHYSICAL domain coordinates
    // In the new design, overrides are directly in physical coordinates (no conversion)
    real phys_override_x0 = 50e-9;   // 50nm physical
    real phys_override_x1 = 350e-9;  // 350nm physical

    // Override region directly in physical coordinates (no pml_offset conversion)
    MeshOverrideRegion override_region;
    override_region.enabled = true;
    override_region.x0 = phys_override_x0;  // 50nm physical
    override_region.x1 = phys_override_x1;  // 350nm physical
    override_region.y0 = phys_override_x0;
    override_region.y1 = phys_override_x1;
    override_region.z0 = 400e-9;
    override_region.z1 = 800e-9;
    override_region.dx = 5e-9;
    override_region.dy = 5e-9;
    override_region.dz = 5e-9;

    config.override_regions.push_back(override_region);

    std::cout << "  Override in physical domain: [" << phys_override_x0 * 1e9 << ", "
              << phys_override_x1 * 1e9 << "] nm (direct, no conversion)\n";

    AutoMeshGenerator generator(config);

    // Physical domain size: 1000nm
    // In the new design, we pass PHYSICAL domain size (not total domain)
    real physical_domain = 1000e-9;

    // Structure at (200nm, 200nm) physical with r=100nm
    // In the new design, structure bounds are also in PHYSICAL coordinates
    std::vector<AutoMeshStructureBounds> structures = {
        {100e-9, 300e-9,   // x: 100-300nm physical
         100e-9, 300e-9,   // y: 100-300nm physical
         500e-9, 700e-9,   // z: 500-700nm physical
         2.0}
    };

    // Generate mesh: Step 1 = physical mesh, Step 2 = PML extension
    GridSpacing grid = generator.generate(physical_domain, physical_domain, physical_domain, npml, structures, 1.0);

    // Test 7a: Verify physical domain starts at index npml
    // After PML extension, physical domain [0, 1000nm] is mapped to indices [npml, N-npml]
    real physical_start_x = grid.x_bounds[npml];
    real physical_start_y = grid.y_bounds[npml];
    real physical_end_x = grid.x_bounds[grid.dx.size() - npml];
    real physical_end_y = grid.y_bounds[grid.dy.size() - npml];

    std::cout << "\n  Physical Domain Boundary Verification:\n";
    std::cout << "    Expected physical domain: [0, " << physical_domain * 1e9 << "] nm\n";
    std::cout << "    Actual X: [" << physical_start_x * 1e9 << ", " << physical_end_x * 1e9 << "] nm\n";
    std::cout << "    Actual Y: [" << physical_start_y * 1e9 << ", " << physical_end_y * 1e9 << "] nm\n";

    // Physical domain size should be 1000nm
    real domain_size_x = physical_end_x - physical_start_x;
    real domain_size_y = physical_end_y - physical_start_y;
    bool domain_size_correct_x = std::abs(domain_size_x - physical_domain) < 1e-9;
    bool domain_size_correct_y = std::abs(domain_size_y - physical_domain) < 1e-9;

    report_test("Physical domain size correct (X)", domain_size_correct_x,
                "Size: " + std::to_string(domain_size_x * 1e9) + " nm");
    report_test("Physical domain size correct (Y)", domain_size_correct_y,
                "Size: " + std::to_string(domain_size_y * 1e9) + " nm");

    // Test 7b: Verify override region boundaries are in grid (in total domain coordinates)
    // Override [50nm, 350nm] physical -> [physical_start + 50nm, physical_start + 350nm] total
    real override_start_total = physical_start_x + phys_override_x0;
    real override_end_total = physical_start_x + phys_override_x1;

    bool override_start_in_grid = position_in_grid(grid.x_bounds, override_start_total, 1e-9);
    bool override_end_in_grid = position_in_grid(grid.x_bounds, override_end_total, 1e-9);

    report_test("Override start boundary in grid", override_start_in_grid,
                "Position: " + std::to_string(override_start_total * 1e9) + " nm (total domain)");
    report_test("Override end boundary in grid", override_end_in_grid,
                "Position: " + std::to_string(override_end_total * 1e9) + " nm (total domain)");

    // Test 7c: Verify override region is at correct physical position
    // Find the fine mesh region (dx < 10nm) within the physical domain

    real fine_region_start = 0;
    real fine_region_end = 0;
    real threshold_dx = 10e-9;  // Cells < 10nm are considered "fine"

    for (size_t i = npml; i < grid.dx.size() - npml; ++i) {
        if (grid.dx[i] < threshold_dx) {
            if (fine_region_start == 0) {
                fine_region_start = grid.x_bounds[i];
            }
            fine_region_end = grid.x_bounds[i + 1];
        }
    }

    // Convert to physical coordinates (relative to physical domain start)
    real fine_region_start_phys = fine_region_start - physical_start_x;
    real fine_region_end_phys = fine_region_end - physical_start_x;

    std::cout << "\n  Fine Mesh Region (dx < 10nm) in Physical Coords:\n";
    std::cout << "    Expected: [" << phys_override_x0 * 1e9 << ", " << phys_override_x1 * 1e9 << "] nm\n";
    std::cout << "    Actual:   [" << fine_region_start_phys * 1e9 << ", " << fine_region_end_phys * 1e9 << "] nm\n";

    // The fine region should start near 50nm and end near 350nm (with some tolerance for grading transition)
    real start_error = std::abs(fine_region_start_phys - phys_override_x0);
    real end_error = std::abs(fine_region_end_phys - phys_override_x1);

    // Allow up to 100nm error for grading transition cells
    // Before the fix, error was ~400nm (PML offset miscalculation)
    bool start_close = start_error < 100e-9;
    bool end_close = end_error < 100e-9;

    report_test("Fine mesh starts near expected position", start_close,
                "Error: " + std::to_string(start_error * 1e9) + " nm (tolerance: 100nm)");
    report_test("Fine mesh ends near expected position", end_close,
                "Error: " + std::to_string(end_error * 1e9) + " nm (tolerance: 100nm)");

    // Test 7d: Verify PML region uses boundary cell size (UNIFORM strategy)
    std::cout << "\n  PML Region Mesh (should use boundary cell size):\n";

    // Get the cell size at physical domain boundary (left side)
    real boundary_dx_left = grid.dx[npml];  // First cell of physical domain
    real boundary_dx_right = grid.dx[grid.dx.size() - npml - 1];  // Last cell of physical domain

    std::cout << "    Left boundary dx: " << boundary_dx_left * 1e9 << " nm\n";
    std::cout << "    Right boundary dx: " << boundary_dx_right * 1e9 << " nm\n";

    std::vector<real> pml_dx_left;
    for (size_t i = 0; i < npml && i < grid.dx.size(); ++i) {
        pml_dx_left.push_back(grid.dx[i]);
    }

    if (!pml_dx_left.empty()) {
        real pml_dx_mean = std::accumulate(pml_dx_left.begin(), pml_dx_left.end(), 0.0) / pml_dx_left.size();
        real pml_dx_var = 0;
        for (real d : pml_dx_left) {
            pml_dx_var = std::max(pml_dx_var, std::abs(d - pml_dx_mean) / pml_dx_mean);
        }

        std::cout << "    Left PML cells: " << pml_dx_left.size() << "\n";
        std::cout << "    Mean dx: " << pml_dx_mean * 1e9 << " nm\n";
        std::cout << "    Max variation: " << pml_dx_var * 100 << "%\n";

        // For UNIFORM strategy, PML should be uniform with boundary cell size
        bool pml_uniform = pml_dx_var < 0.05;  // < 5% variation
        bool pml_matches_boundary = std::abs(pml_dx_mean - boundary_dx_left) / boundary_dx_left < 0.1;

        report_test("PML region is uniform", pml_uniform,
                    "Variation: " + std::to_string(pml_dx_var * 100) + "%");
        report_test("PML matches boundary cell size", pml_matches_boundary,
                    "PML mean: " + std::to_string(pml_dx_mean * 1e9) + "nm, Boundary: " + std::to_string(boundary_dx_left * 1e9) + "nm");
    }

    return domain_size_correct_x && domain_size_correct_y && start_close && end_close;
}

// ==================== Main ====================

int main() {
    std::cout << BOLD << "\n"
              << "======================================================================\n"
              << "     AUTO NON-UNIFORM MESH GENERATOR - UNIT TEST SUITE               \n"
              << "======================================================================\n"
              << RESET;

    bool t1 = test_two_material_interface();
    bool t2 = test_single_override_box();
    bool t3 = test_multiple_overrides();
    bool t4 = test_global_grading_ratio();
    bool t5 = test_ppw_mapping();
    bool t6 = test_demo_output();
    bool t7 = test_pml_override_alignment();

    std::cout << BOLD << "\n"
              << "======================================================================\n"
              << "                       TEST SUMMARY                                   \n"
              << "======================================================================\n"
              << RESET;

    int passed = 0, failed = 0;
    for (const auto& r : all_results) {
        if (r.passed) passed++;
        else failed++;
    }

    std::cout << "\nTotal tests: " << all_results.size() << "\n";
    std::cout << GREEN << "Passed: " << passed << RESET << "\n";
    if (failed > 0) {
        std::cout << RED << "Failed: " << failed << RESET << "\n";
        std::cout << "\nFailed tests:\n";
        for (const auto& r : all_results) {
            if (!r.passed) {
                std::cout << RED << "  - " << r.name << RESET;
                if (!r.message.empty()) std::cout << ": " << r.message;
                std::cout << "\n";
            }
        }
    }

    bool all_passed = (failed == 0);
    std::cout << "\n" << (all_passed ? GREEN : RED) << BOLD
              << "Overall: " << (all_passed ? "ALL TESTS PASSED" : "SOME TESTS FAILED")
              << RESET << "\n\n";

    return all_passed ? 0 : 1;
}
