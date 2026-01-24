// auto_mesh_generator.hpp - Automatic Non-uniform Mesh Generator
//
// Two-step mesh generation approach (Physical first, PML later):
//
// Step 1: Generate Physical Domain mesh
//   - PPW-based local cell size: D <= lambda_min / (n_local * ppw_target)
//   - Interface snapping: material boundaries are grid lines
//   - Mesh override regions: user-defined uniform mesh (in physical coordinates)
//   - Grading constraint: adjacent cell ratio <= sqrt(2)
//
// Step 2: Extend with PML
//   - Read boundary cell size from physical mesh edges
//   - Extend npml cells on each side
//   - Strategy A (default): uniform dx = boundary dx (most stable for CPML)
//   - Strategy B (optional): graded PML respecting ratio <= g
//
// Key design principle:
//   Override regions are ALWAYS defined in physical domain coordinates.
//   No coordinate conversion needed - override applies directly to physical mesh.
//   PML is added as extension, so physical domain starts at index npml.

#pragma once

#include <vector>
#include <cmath>
#include <algorithm>
#include <set>
#include <map>
#include <iostream>
#include <limits>
#include <numeric>
#include <cassert>

#include "global_function.hpp"

// Forward declarations
struct Material;
struct Shape;
struct StructureItem;

// StructureBounds - defined here if auto_mesh_generator.hpp is not included
#ifndef AUTO_MESH_GENERATOR_HPP_STRUCTURE_BOUNDS_DEFINED
#define AUTO_MESH_GENERATOR_HPP_STRUCTURE_BOUNDS_DEFINED
struct StructureBounds {
    real x_min, x_max;
    real y_min, y_max;
    real z_min, z_max;
    real n_material;
};
#endif

// ========== Mesh Override Region ==========
// User-defined uniform mesh region - coordinates are in PHYSICAL domain
struct MeshOverrideRegion {
    bool enabled = false;
    real x0 = 0, x1 = 0;               // X bounds (physical domain coordinates)
    real y0 = 0, y1 = 0;               // Y bounds (physical domain coordinates)
    real z0 = 0, z1 = 0;               // Z bounds (physical domain coordinates)
    real dx = 10e-9;                   // Target mesh size in X direction
    real dy = 10e-9;                   // Target mesh size in Y direction
    real dz = 10e-9;                   // Target mesh size in Z direction
    int priority = 0;                  // Higher priority overrides lower (default: finer wins)
};

// ========== PML Extension Strategy ==========
enum class PMLStrategy {
    UNIFORM,    // Strategy A: PML uses uniform dx = boundary dx (most stable)
    GRADED      // Strategy B: PML uses graded dx respecting ratio <= g
};

// ========== Auto Mesh Configuration ==========
struct AutoMeshConfig {
    // Wavelength parameters
    real lambda_min = 300e-9;           // Minimum wavelength (m) = c0 / f_max
    real f_max = 0;                     // Maximum frequency (Hz), alternative to lambda_min

    // Mesh accuracy level (determines PPW)
    // acc=1 -> 6 ppw, acc=2 -> 10 ppw, acc=3 -> 14 ppw, acc>=3: 14 + (acc-3)*4
    int mesh_accuracy = 2;              // Default: accuracy level 2 -> 10 ppw

    // Direct PPW override (if > 0, overrides mesh_accuracy)
    real ppw_override = 0;              // Points per wavelength override

    // Mesh size limits
    real dx_min = 1e-9;                 // Minimum allowed spacing (m)
    real dx_max = 100e-9;               // Maximum allowed spacing (m)

    // Grading control
    real max_grading_ratio = 1.41421356237;  // sqrt(2), max ratio between adjacent cells

    // PML configuration
    size_t npml = 0;                    // Number of PML cells
    PMLStrategy pml_strategy = PMLStrategy::UNIFORM;  // How to generate PML mesh

    // Mesh override regions (coordinates are in PHYSICAL domain)
    bool mesh_override_enabled = false;
    std::vector<MeshOverrideRegion> override_regions;

    // Get effective PPW based on mesh_accuracy
    real get_ppw_target() const {
        if (ppw_override > 0) return ppw_override;

        // Mapping: acc=1 -> 6, acc=2 -> 10, acc=3 -> 14, acc>=3: 14 + (acc-3)*4
        if (mesh_accuracy <= 1) return 6.0;
        if (mesh_accuracy == 2) return 10.0;
        return 14.0 + (mesh_accuracy - 3) * 4.0;
    }

    // Get lambda_min from f_max if specified
    real get_lambda_min() const {
        if (f_max > 0) {
            return PhysConst::C0 / f_max;
        }
        return lambda_min;
    }
};

// ========== Auto Mesh Structure Bounding Box ==========
// Coordinates are in PHYSICAL domain
struct AutoMeshStructureBounds {
    real x_min, x_max;
    real y_min, y_max;
    real z_min, z_max;
    real n_material;     // Refractive index of structure
};

// ========== 1D Constraint Point ==========
struct GridConstraint1D {
    real position;
    bool is_interface;      // True if this is a material interface
    bool is_override_bound; // True if this is an override region boundary
    real n_left;            // Refractive index on left side
    real n_right;           // Refractive index on right side
};

// ========== 1D Override Segment ==========
struct OverrideSegment1D {
    real start;
    real end;
    real dx_override;
    int priority;
};

// ========== Auto Mesh Generator Class ==========
class AutoMeshGenerator {
public:
    AutoMeshConfig config;

    AutoMeshGenerator() = default;
    explicit AutoMeshGenerator(const AutoMeshConfig& cfg) : config(cfg) {}

    // ===== Main entry point: Generate complete 3D grid =====
    // physical_domain_x/y/z: size of the physical domain (excluding PML)
    // npml: number of PML cells to add on each side
    // structures: geometry bounds (in physical coordinates)
    GridSpacing generate(
        real physical_domain_x, real physical_domain_y, real physical_domain_z,
        size_t npml,
        const std::vector<AutoMeshStructureBounds>& structures,
        real background_n = 1.0
    ) {
        config.npml = npml;

        std::cout << "\n[AutoMesh] ========================================\n";
        std::cout << "[AutoMesh] Two-Step Mesh Generation (Physical first, PML later)\n";
        std::cout << "[AutoMesh] Physical domain: " << physical_domain_x*1e9 << " x "
                  << physical_domain_y*1e9 << " x " << physical_domain_z*1e9 << " nm\n";
        std::cout << "[AutoMesh] Lambda_min: " << config.get_lambda_min()*1e9 << " nm\n";
        std::cout << "[AutoMesh] Mesh accuracy: " << config.mesh_accuracy
                  << " (PPW target: " << config.get_ppw_target() << ")\n";
        std::cout << "[AutoMesh] Max grading ratio: " << config.max_grading_ratio << "\n";
        std::cout << "[AutoMesh] Structure count: " << structures.size() << "\n";
        std::cout << "[AutoMesh] PML: npml=" << npml << ", strategy="
                  << (config.pml_strategy == PMLStrategy::UNIFORM ? "UNIFORM" : "GRADED") << "\n";

        // Print mesh override regions (in physical coordinates)
        if (config.mesh_override_enabled && !config.override_regions.empty()) {
            size_t active_count = 0;
            for (const auto& r : config.override_regions) {
                if (r.enabled) active_count++;
            }
            std::cout << "[AutoMesh] Mesh Override Regions (physical coords): " << active_count << " active\n";
            for (size_t i = 0; i < config.override_regions.size(); ++i) {
                const auto& r = config.override_regions[i];
                if (r.enabled) {
                    std::cout << "  Region " << i << ": ["
                              << r.x0*1e9 << "-" << r.x1*1e9 << "] x ["
                              << r.y0*1e9 << "-" << r.y1*1e9 << "] x ["
                              << r.z0*1e9 << "-" << r.z1*1e9 << "] nm, "
                              << "dx=" << r.dx*1e9 << ", dy=" << r.dy*1e9
                              << ", dz=" << r.dz*1e9 << " nm\n";
                }
            }
        } else {
            std::cout << "[AutoMesh] Mesh Override Regions: disabled\n";
        }

        // ===== STEP 1: Generate Physical Domain mesh =====
        std::cout << "\n[AutoMesh] === Step 1: Generating Physical Domain mesh ===\n";

        std::cout << "[AutoMesh] Generating X physical grid...\n";
        auto x_phys_lines = generate_physical_mesh_1d(structures, 0, physical_domain_x, background_n);

        std::cout << "[AutoMesh] Generating Y physical grid...\n";
        auto y_phys_lines = generate_physical_mesh_1d(structures, 1, physical_domain_y, background_n);

        std::cout << "[AutoMesh] Generating Z physical grid...\n";
        auto z_phys_lines = generate_physical_mesh_1d(structures, 2, physical_domain_z, background_n);

        // Print physical mesh stats
        print_physical_mesh_stats(x_phys_lines, y_phys_lines, z_phys_lines);

        // Verify override boundaries in physical mesh
        verify_override_boundaries(x_phys_lines, y_phys_lines, z_phys_lines);

        // ===== STEP 2: Extend with PML =====
        std::cout << "\n[AutoMesh] === Step 2: Extending with PML ===\n";

        auto x_total_lines = extend_with_pml_1d(x_phys_lines, npml, "X");
        auto y_total_lines = extend_with_pml_1d(y_phys_lines, npml, "Y");
        auto z_total_lines = extend_with_pml_1d(z_phys_lines, npml, "Z");

        // Build GridSpacing from total grid lines
        GridSpacing grid;
        grid.x_bounds = x_total_lines;
        grid.y_bounds = y_total_lines;
        grid.z_bounds = z_total_lines;

        // Compute dx/dy/dz arrays from bounds
        grid.dx.resize(x_total_lines.size() - 1);
        grid.dy.resize(y_total_lines.size() - 1);
        grid.dz.resize(z_total_lines.size() - 1);

        for (size_t i = 0; i < grid.dx.size(); ++i) {
            grid.dx[i] = x_total_lines[i + 1] - x_total_lines[i];
        }
        for (size_t i = 0; i < grid.dy.size(); ++i) {
            grid.dy[i] = y_total_lines[i + 1] - y_total_lines[i];
        }
        for (size_t i = 0; i < grid.dz.size(); ++i) {
            grid.dz[i] = z_total_lines[i + 1] - z_total_lines[i];
        }

        grid.npml = npml;

        // Print final statistics
        print_final_mesh_stats(grid, npml, x_phys_lines.size() - 1,
                               y_phys_lines.size() - 1, z_phys_lines.size() - 1);

        std::cout << "[AutoMesh] ========================================\n\n";

        return grid;
    }

private:
    // ===== Step 1: Generate Physical Domain mesh for one axis =====
    std::vector<real> generate_physical_mesh_1d(
        const std::vector<AutoMeshStructureBounds>& structures,
        int axis,  // 0=x, 1=y, 2=z
        real domain_size,
        real background_n
    ) {
        // Collect hard constraints (must-include grid lines)
        auto constraints = collect_constraints_1d(structures, axis, domain_size, background_n);

        // Collect override segments
        auto overrides = collect_override_segments_1d(axis, domain_size);

        // Build grid using segment-based approach
        auto grid_lines = build_segmented_grid_1d(constraints, overrides, structures,
                                                   axis, domain_size, background_n);

        // Validate
        validate_grid_1d(grid_lines, constraints, overrides, axis, "physical");

        return grid_lines;
    }

    // ===== Step 2: Extend physical mesh with PML =====
    std::vector<real> extend_with_pml_1d(
        const std::vector<real>& phys_lines,
        size_t npml,
        const char* axis_name
    ) {
        if (npml == 0 || phys_lines.size() < 2) {
            return phys_lines;
        }

        // Read boundary cell sizes from physical mesh
        real dx_left = phys_lines[1] - phys_lines[0];
        real dx_right = phys_lines.back() - phys_lines[phys_lines.size() - 2];

        std::vector<real> left_pml_lines;
        std::vector<real> right_pml_lines;

        if (config.pml_strategy == PMLStrategy::UNIFORM) {
            // Strategy A: Uniform PML (most stable for CPML)
            // Each PML cell has the same size as the boundary cell

            // Left PML: extends from -npml*dx_left to 0
            real left_pml_thickness = npml * dx_left;
            for (size_t i = 0; i < npml; ++i) {
                real pos = -left_pml_thickness + i * dx_left;
                left_pml_lines.push_back(pos);
            }
            // Note: phys_lines[0] = 0 will be the connecting point

            // Right PML: extends from phys_end to phys_end + npml*dx_right
            real phys_end = phys_lines.back();
            for (size_t i = 1; i <= npml; ++i) {
                real pos = phys_end + i * dx_right;
                right_pml_lines.push_back(pos);
            }

        } else {
            // Strategy B: Graded PML (allows coarsening, more efficient)
            const real g = config.max_grading_ratio;

            // Left PML: start from boundary, grow outward
            real dx = dx_left;
            real pos = 0;
            std::vector<real> left_tmp;
            for (size_t i = 0; i < npml; ++i) {
                pos -= dx;
                left_tmp.push_back(pos);
                dx = std::min(dx * g, config.dx_max);
            }
            // Reverse to get increasing order
            for (auto it = left_tmp.rbegin(); it != left_tmp.rend(); ++it) {
                left_pml_lines.push_back(*it);
            }

            // Right PML: start from boundary, grow outward
            dx = dx_right;
            pos = phys_lines.back();
            for (size_t i = 0; i < npml; ++i) {
                pos += dx;
                right_pml_lines.push_back(pos);
                dx = std::min(dx * g, config.dx_max);
            }
        }

        // Combine: left_pml + phys + right_pml
        std::vector<real> total_lines;
        total_lines.reserve(left_pml_lines.size() + phys_lines.size() + right_pml_lines.size());

        for (real x : left_pml_lines) {
            total_lines.push_back(x);
        }
        for (real x : phys_lines) {
            total_lines.push_back(x);
        }
        for (real x : right_pml_lines) {
            total_lines.push_back(x);
        }

        // Shift so that total domain starts at 0
        real offset = -total_lines[0];
        for (real& x : total_lines) {
            x += offset;
        }

        // Report
        real left_pml_thickness = total_lines[npml] - total_lines[0];
        real right_pml_thickness = total_lines.back() - total_lines[total_lines.size() - 1 - npml];
        real phys_start = total_lines[npml];
        real phys_end = total_lines[total_lines.size() - 1 - npml];

        std::cout << "  " << axis_name << ": left_PML=" << left_pml_thickness*1e9
                  << "nm (" << npml << " cells, dx=" << dx_left*1e9 << "nm), "
                  << "right_PML=" << right_pml_thickness*1e9
                  << "nm (" << npml << " cells, dx=" << dx_right*1e9 << "nm)\n";
        std::cout << "      Physical domain in total: [" << phys_start*1e9 << ", "
                  << phys_end*1e9 << "] nm\n";

        return total_lines;
    }

    // ===== Build segmented grid for physical domain =====
    // ALGORITHM v6: Simple segment-based approach
    //
    // 1. Divide domain into segments based on override regions
    // 2. For each segment, generate cells with proper grading
    // 3. Post-process to fix any remaining violations
    //
    std::vector<real> build_segmented_grid_1d(
        const std::vector<GridConstraint1D>& constraints,
        const std::vector<OverrideSegment1D>& overrides,
        const std::vector<AutoMeshStructureBounds>& structures,
        int axis,
        real domain_size,
        real background_n
    ) {
        const real g = config.max_grading_ratio;

        // Helper: get auto dx at position (based on material)
        auto get_auto_dx = [&](real pos) -> real {
            real n_local = sample_n_at_position(structures, axis, pos, background_n);
            real lambda_eff = config.get_lambda_min() / n_local;
            real auto_dx = lambda_eff / config.get_ppw_target();
            return std::clamp(auto_dx, config.dx_min, config.dx_max);
        };

        // Build list of regions: [start, end, dx, is_override]
        struct Region {
            real start, end, target_dx;
            bool is_override;
        };
        std::vector<Region> regions;

        // Collect all boundary points (override regions + structure interfaces)
        std::set<real> boundaries;
        boundaries.insert(0.0);
        boundaries.insert(domain_size);

        // Add override region boundaries
        for (const auto& ov : overrides) {
            boundaries.insert(std::max(0.0, std::min(ov.start, domain_size)));
            boundaries.insert(std::max(0.0, std::min(ov.end, domain_size)));
        }

        // Add structure interface boundaries (for proper boundary snapping)
        for (const auto& s : structures) {
            real s_min, s_max;
            switch (axis) {
                case 0: s_min = s.x_min; s_max = s.x_max; break;
                case 1: s_min = s.y_min; s_max = s.y_max; break;
                case 2: s_min = s.z_min; s_max = s.z_max; break;
                default: continue;
            }
            if (s_min > 0 && s_min < domain_size) {
                boundaries.insert(s_min);
            }
            if (s_max > 0 && s_max < domain_size) {
                boundaries.insert(s_max);
            }
        }

        // Create regions from boundaries
        std::vector<real> bounds(boundaries.begin(), boundaries.end());
        for (size_t i = 0; i + 1 < bounds.size(); ++i) {
            real start = bounds[i];
            real end = bounds[i + 1];
            real mid = (start + end) / 2;

            // Check if this region is an override
            bool is_override = false;
            real override_dx = 0;
            for (const auto& ov : overrides) {
                if (mid >= ov.start && mid <= ov.end) {
                    is_override = true;
                    override_dx = ov.dx_override;
                    break;
                }
            }

            real target_dx = is_override ? override_dx : get_auto_dx(mid);
            regions.push_back({start, end, target_dx, is_override});
        }

        // Generate grid
        std::vector<real> grid;
        grid.push_back(0.0);
        real prev_dx = regions[0].target_dx;

        for (size_t ri = 0; ri < regions.size(); ++ri) {
            const auto& reg = regions[ri];
            real pos = grid.back();
            real target = reg.end;

            if (reg.is_override) {
                // OVERRIDE REGION: strictly uniform cells
                // Transitions happen in adjacent non-override regions

                real ov_length = reg.end - reg.start;
                size_t n_cells = std::max(size_t(1), (size_t)std::round(ov_length / reg.target_dx));
                real actual_dx = ov_length / n_cells;

                // Generate uniform cells from reg.start to reg.end
                for (size_t i = 1; i <= n_cells; ++i) {
                    real cell_pos = reg.start + i * actual_dx;
                    if (i == n_cells) cell_pos = reg.end;  // Snap to exact boundary
                    grid.push_back(cell_pos);
                }
                prev_dx = actual_dx;
                pos = reg.end;

            } else {
                // NON-OVERRIDE REGION: handle transitions properly

                bool prev_is_override = (ri > 0 && regions[ri - 1].is_override);
                bool next_is_override = (ri + 1 < regions.size() && regions[ri + 1].is_override);

                real region_length = reg.end - reg.start;
                real this_region_target_dx = reg.target_dx;

                // CASE 1: Sandwiched between two adjacent override regions - use graded mesh
                // Grade up from prev override dx, reach peak, grade down to next override dx
                if (prev_is_override && next_is_override) {
                    real prev_ov_dx = regions[ri - 1].target_dx;
                    real next_ov_dx = regions[ri + 1].target_dx;
                    real start_dx = prev_dx;
                    real end_dx = next_ov_dx;

                    // Calculate peak dx based on gap length and material
                    real half_gap = region_length / 2;
                    real peak_dx = start_dx;
                    real temp_length = 0;
                    while (temp_length < half_gap && peak_dx < this_region_target_dx) {
                        peak_dx *= g;
                        temp_length += peak_dx;
                    }
                    peak_dx = std::min(peak_dx, this_region_target_dx);

                    // Build graded sequence: up phase + down phase
                    std::vector<real> up_phase, down_phase;

                    // Grade up from start_dx to peak
                    real temp_dx = start_dx;
                    while (temp_dx < peak_dx * 0.95) {
                        temp_dx = std::min(temp_dx * g, peak_dx);
                        up_phase.push_back(temp_dx);
                    }

                    // Grade down from peak to end_dx
                    temp_dx = peak_dx;
                    while (temp_dx > end_dx * 1.05) {
                        temp_dx = std::max(temp_dx / g, end_dx);
                        down_phase.push_back(temp_dx);
                    }

                    // Calculate total graded length
                    real graded_length = 0;
                    for (real dx : up_phase) graded_length += dx;
                    for (real dx : down_phase) graded_length += dx;

                    // If graded fits, use it; otherwise use simpler approach
                    if (graded_length > 0 && graded_length <= region_length * 1.3) {
                        real scale = region_length / graded_length;

                        std::vector<real> dx_sequence;
                        dx_sequence.insert(dx_sequence.end(), up_phase.begin(), up_phase.end());
                        dx_sequence.insert(dx_sequence.end(), down_phase.begin(), down_phase.end());

                        real current_pos = reg.start;
                        for (size_t i = 0; i < dx_sequence.size(); ++i) {
                            real scaled_dx = dx_sequence[i] * scale;
                            current_pos += scaled_dx;
                            if (i == dx_sequence.size() - 1) current_pos = reg.end;
                            grid.push_back(current_pos);
                        }
                        prev_dx = dx_sequence.back() * scale;
                        continue;
                    }

                    // Fallback: if grading doesn't work, use uniform at finer dx
                    real fine_dx = std::min(start_dx, end_dx);
                    size_t n_cells = std::max(size_t(1), (size_t)std::round(region_length / fine_dx));
                    real actual_dx = region_length / n_cells;
                    for (size_t i = 1; i <= n_cells; ++i) {
                        real cell_pos = reg.start + i * actual_dx;
                        if (i == n_cells) cell_pos = reg.end;
                        grid.push_back(cell_pos);
                    }
                    prev_dx = actual_dx;
                    continue;
                }

                // CASE 2: Structure region (not air) - use uniform mesh at target dx
                // This prevents over-refinement in structure regions
                // A structure region has target dx significantly below the air target (dx_max)
                // and is not close to the override dx
                real air_target_dx = config.get_lambda_min() / background_n / config.get_ppw_target();
                bool is_structure_region = (this_region_target_dx < air_target_dx * 0.6) &&
                                           (this_region_target_dx > config.dx_min * 2);
                bool coming_from_finer = (prev_dx < this_region_target_dx * 0.7);

                if (is_structure_region && !next_is_override) {
                    // For structure regions, use uniform mesh at the material's target dx
                    // But first, handle transition from finer mesh if needed

                    if (coming_from_finer) {
                        // Need to grade UP from prev_dx toward this_region_target_dx
                        // Calculate how many grading steps needed
                        std::vector<real> transition_dx;
                        real temp_dx = prev_dx;
                        while (temp_dx < this_region_target_dx * 0.95) {
                            temp_dx *= g;
                            if (temp_dx > this_region_target_dx) temp_dx = this_region_target_dx;
                            transition_dx.push_back(temp_dx);
                        }

                        // Calculate transition length needed
                        real trans_length = 0;
                        for (real dx : transition_dx) trans_length += dx;

                        if (trans_length < region_length) {
                            // Add transition cells
                            for (real dx : transition_dx) {
                                pos += dx;
                                if (pos > reg.end - config.dx_min) break;
                                grid.push_back(pos);
                                prev_dx = dx;
                            }
                        }
                    }

                    // Fill remaining region with uniform cells at target dx
                    real remaining = reg.end - pos;
                    if (remaining > config.dx_min) {
                        // Use target dx, but respect grading from last cell
                        real use_dx = this_region_target_dx;
                        if (use_dx / prev_dx > g) use_dx = prev_dx * g;
                        if (prev_dx / use_dx > g) use_dx = prev_dx / g;

                        size_t n_cells = std::max(size_t(1), (size_t)std::round(remaining / use_dx));
                        real actual_dx = remaining / n_cells;

                        // Verify grading constraint
                        if (actual_dx / prev_dx > g || prev_dx / actual_dx > g) {
                            // Adjust to satisfy grading
                            if (prev_dx < actual_dx) {
                                actual_dx = std::min(actual_dx, prev_dx * g);
                            } else {
                                actual_dx = std::max(actual_dx, prev_dx / g);
                            }
                            n_cells = std::max(size_t(1), (size_t)std::round(remaining / actual_dx));
                            actual_dx = remaining / n_cells;
                        }

                        for (size_t i = 1; i <= n_cells; ++i) {
                            real cell_pos = pos + i * actual_dx;
                            if (i == n_cells) cell_pos = reg.end;
                            grid.push_back(cell_pos);
                        }
                        prev_dx = actual_dx;
                    }
                    continue;
                }

                // CASE 3: Air region or transitioning to override - use graded mesh
                //
                // Key insight: Build an explicit dx sequence that:
                // 1. Starts from entry_dx (inherited from previous region)
                // 2. Grows toward target_dx respecting grading constraint
                // 3. Stays at stable_dx (close to target_dx) as long as possible
                // 4. Shrinks toward exit_target_dx if needed
                // 5. Ends close to exit_target_dx to avoid grading violation at boundary

                real entry_dx = prev_dx;
                real exit_target_dx = reg.target_dx;
                bool need_exit_transition = false;

                if (next_is_override) {
                    exit_target_dx = regions[ri + 1].target_dx;
                    need_exit_transition = true;
                } else if (ri + 1 < regions.size()) {
                    real next_target = regions[ri + 1].target_dx;
                    if (next_target < this_region_target_dx) {
                        exit_target_dx = next_target;
                        need_exit_transition = true;
                    }
                }

                // Build explicit dx sequences for entry and exit transitions
                std::vector<real> entry_seq;  // Growing sequence
                std::vector<real> exit_seq;   // Shrinking sequence

                // Entry sequence: grow from entry_dx to target_dx
                if (entry_dx < this_region_target_dx * 0.99) {
                    real temp_dx = entry_dx;
                    while (temp_dx < this_region_target_dx * 0.99) {
                        temp_dx = std::min(temp_dx * g, this_region_target_dx);
                        temp_dx = std::clamp(temp_dx, config.dx_min, config.dx_max);
                        entry_seq.push_back(temp_dx);
                    }
                }

                // Exit sequence: shrink from target_dx to exit_target_dx
                if (need_exit_transition && exit_target_dx < this_region_target_dx * 0.99) {
                    real temp_dx = this_region_target_dx;
                    while (temp_dx > exit_target_dx * 1.01) {
                        temp_dx = std::max(temp_dx / g, exit_target_dx);
                        temp_dx = std::clamp(temp_dx, config.dx_min, config.dx_max);
                        exit_seq.push_back(temp_dx);
                    }
                    // Ensure last element is close to exit_target_dx
                    if (exit_seq.empty() || exit_seq.back() > exit_target_dx * 1.2) {
                        exit_seq.push_back(exit_target_dx);
                    }
                }

                // Calculate lengths
                real entry_length = 0;
                for (real dx : entry_seq) entry_length += dx;
                real exit_length = 0;
                for (real dx : exit_seq) exit_length += dx;

                real stable_length = region_length - entry_length - exit_length;

                // If not enough room, scale down to fit
                if (stable_length < 0) {
                    // Need to squeeze - reduce stable region and use finer mesh
                    real total_trans = entry_length + exit_length;
                    real scale = region_length / total_trans * 0.95;  // Leave some margin

                    // Scale both sequences
                    for (real& dx : entry_seq) dx *= scale;
                    for (real& dx : exit_seq) dx *= scale;

                    entry_length *= scale;
                    exit_length *= scale;
                    stable_length = region_length - entry_length - exit_length;
                }

                // Determine stable dx - should always be the target dx for this region
                // Transitions happen in entry/exit phases, not in stable region
                real stable_dx = std::clamp(this_region_target_dx, config.dx_min, config.dx_max);

                // Generate the grid
                real current_dx = entry_dx;

                // Phase 1: Entry transition
                for (real dx : entry_seq) {
                    // Ensure grading constraint
                    real next_dx = std::clamp(dx, current_dx / g, current_dx * g);
                    next_dx = std::clamp(next_dx, config.dx_min, config.dx_max);

                    if (pos + next_dx > target - config.dx_min) break;

                    pos += next_dx;
                    grid.push_back(pos);
                    current_dx = next_dx;
                    prev_dx = next_dx;
                }

                // Phase 2: Stable region
                if (stable_length > stable_dx * 0.5 && pos < target - exit_length - config.dx_min) {
                    real stable_end_pos = target - exit_length;
                    real remaining_stable = stable_end_pos - pos;

                    if (remaining_stable > 0) {
                        // Ensure grading from last entry cell
                        real use_dx = stable_dx;
                        if (use_dx / current_dx > g) use_dx = current_dx * g;
                        if (current_dx / use_dx > g) use_dx = current_dx / g;
                        use_dx = std::clamp(use_dx, config.dx_min, config.dx_max);

                        size_t n_stable = std::max(size_t(1), (size_t)std::round(remaining_stable / use_dx));
                        real actual_stable_dx = remaining_stable / n_stable;

                        for (size_t i = 1; i <= n_stable; ++i) {
                            real cell_pos = pos + i * actual_stable_dx;
                            if (i == n_stable) cell_pos = stable_end_pos;
                            grid.push_back(cell_pos);
                        }
                        pos = stable_end_pos;
                        current_dx = actual_stable_dx;
                        prev_dx = actual_stable_dx;
                    }
                }

                // Phase 3: Exit transition
                if (!exit_seq.empty() && pos < target - config.dx_min * 0.5) {
                    for (size_t i = 0; i < exit_seq.size(); ++i) {
                        real dx = exit_seq[i];

                        // Ensure grading constraint
                        real next_dx = std::clamp(dx, current_dx / g, current_dx * g);
                        next_dx = std::clamp(next_dx, config.dx_min, config.dx_max);

                        real remaining = target - pos;
                        if (remaining < next_dx * 0.5) break;

                        // Last cell should end at target
                        if (i == exit_seq.size() - 1 || remaining < next_dx * 1.5) {
                            // Fill remaining with appropriate dx
                            int n = std::max(1, (int)std::round(remaining / next_dx));
                            real final_dx = remaining / n;

                            // Check grading
                            if (final_dx / current_dx > g || current_dx / final_dx > g) {
                                n = std::max(2, (int)std::ceil(remaining / (current_dx * g)));
                                final_dx = remaining / n;
                            }

                            for (int j = 0; j < n; ++j) {
                                pos += final_dx;
                                grid.push_back(pos);
                            }
                            prev_dx = final_dx;
                            break;
                        }

                        pos += next_dx;
                        grid.push_back(pos);
                        current_dx = next_dx;
                        prev_dx = next_dx;
                    }
                } else if (exit_seq.empty() && pos < target - config.dx_min * 0.5) {
                    // No exit transition needed - fill with stable dx
                    while (pos < target - config.dx_min * 0.5) {
                        real remaining = target - pos;
                        real next_dx = stable_dx;

                        // Ensure grading
                        next_dx = std::clamp(next_dx, current_dx / g, current_dx * g);
                        next_dx = std::clamp(next_dx, config.dx_min, config.dx_max);

                        if (remaining < next_dx * 1.5) {
                            int n = std::max(1, (int)std::round(remaining / next_dx));
                            real final_dx = remaining / n;
                            if (final_dx / current_dx > g || current_dx / final_dx > g) {
                                n = std::max(2, (int)std::ceil(remaining / (current_dx * g)));
                                final_dx = remaining / n;
                            }
                            for (int j = 0; j < n; ++j) {
                                pos += final_dx;
                                grid.push_back(pos);
                            }
                            prev_dx = final_dx;
                            break;
                        }

                        pos += next_dx;
                        grid.push_back(pos);
                        current_dx = next_dx;
                        prev_dx = next_dx;
                    }
                }

                // Snap to target if needed
                if (std::abs(grid.back() - target) > config.dx_min * 0.1 && grid.back() < target) {
                    grid.push_back(target);
                    prev_dx = target - grid[grid.size() - 2];
                }
            }
        }

        // Snap to domain end
        grid.back() = domain_size;

        // Post-process to fix any remaining grading violations
        grid = fix_grading_by_subdivision(grid, overrides, g);

        // Remove duplicates
        std::vector<real> clean;
        clean.push_back(grid[0]);
        for (size_t i = 1; i < grid.size(); ++i) {
            if (grid[i] > clean.back() + config.dx_min * 0.001) {
                clean.push_back(grid[i]);
            }
        }

        return clean;
    }

    // Fix grading violations by subdividing cells
    std::vector<real> fix_grading_by_subdivision(
        std::vector<real> grid,
        const std::vector<OverrideSegment1D>& overrides,
        real g
    ) {
        const int max_iter = 300;
        const real tol = config.dx_min * 0.01;

        for (int iter = 0; iter < max_iter; ++iter) {
            // Find worst violation
            real worst_ratio = 1.0;
            size_t worst_i = 0;

            for (size_t i = 0; i + 2 < grid.size(); ++i) {
                real dx1 = grid[i + 1] - grid[i];
                real dx2 = grid[i + 2] - grid[i + 1];
                if (dx1 < tol || dx2 < tol) continue;
                real ratio = std::max(dx1 / dx2, dx2 / dx1);
                if (ratio > worst_ratio) {
                    worst_ratio = ratio;
                    worst_i = i;
                }
            }

            if (worst_ratio <= g * 1.005) break;

            // Subdivide the larger cell
            real dx1 = grid[worst_i + 1] - grid[worst_i];
            real dx2 = grid[worst_i + 2] - grid[worst_i + 1];
            size_t split_i = (dx1 > dx2) ? worst_i : worst_i + 1;

            real split_dx = grid[split_i + 1] - grid[split_i];
            if (split_dx < config.dx_min * 1.5) break;

            real mid = 0.5 * (grid[split_i] + grid[split_i + 1]);
            grid.insert(grid.begin() + split_i + 1, mid);
        }

        return grid;
    }


    // ===== Collect hard constraints for physical mesh =====
    std::vector<GridConstraint1D> collect_constraints_1d(
        const std::vector<AutoMeshStructureBounds>& structures,
        int axis,
        real domain_size,
        real background_n
    ) {
        std::vector<GridConstraint1D> constraints;

        // Domain boundaries are always constraints
        constraints.push_back({0.0, false, false, background_n, background_n});
        constraints.push_back({domain_size, false, false, background_n, background_n});

        // Material interface positions
        std::set<real> interface_positions;
        for (const auto& s : structures) {
            real s_min, s_max;
            switch (axis) {
                case 0: s_min = s.x_min; s_max = s.x_max; break;
                case 1: s_min = s.y_min; s_max = s.y_max; break;
                case 2: s_min = s.z_min; s_max = s.z_max; break;
                default: continue;
            }

            if (s_min > 0 && s_min < domain_size) {
                interface_positions.insert(s_min);
            }
            if (s_max > 0 && s_max < domain_size) {
                interface_positions.insert(s_max);
            }
        }

        for (real pos : interface_positions) {
            real n_left = sample_n_at_position(structures, axis, pos - 1e-12, background_n);
            real n_right = sample_n_at_position(structures, axis, pos + 1e-12, background_n);
            constraints.push_back({pos, true, false, n_left, n_right});
        }

        // Override region boundaries (in physical coordinates)
        if (config.mesh_override_enabled) {
            for (const auto& r : config.override_regions) {
                if (!r.enabled) continue;

                real r_min, r_max;
                switch (axis) {
                    case 0: r_min = r.x0; r_max = r.x1; break;
                    case 1: r_min = r.y0; r_max = r.y1; break;
                    case 2: r_min = r.z0; r_max = r.z1; break;
                    default: continue;
                }

                r_min = std::max(0.0, std::min(r_min, domain_size));
                r_max = std::max(0.0, std::min(r_max, domain_size));

                if (r_min < r_max) {
                    real n_at_min = sample_n_at_position(structures, axis, r_min, background_n);
                    real n_at_max = sample_n_at_position(structures, axis, r_max, background_n);

                    if (r_min > 0) {
                        constraints.push_back({r_min, false, true, n_at_min, n_at_min});
                    }
                    if (r_max < domain_size) {
                        constraints.push_back({r_max, false, true, n_at_max, n_at_max});
                    }
                }
            }
        }

        // Sort and merge duplicates
        std::sort(constraints.begin(), constraints.end(),
                  [](const GridConstraint1D& a, const GridConstraint1D& b) {
                      return a.position < b.position;
                  });

        std::vector<GridConstraint1D> merged;
        for (const auto& c : constraints) {
            if (merged.empty() || std::abs(c.position - merged.back().position) > 1e-15) {
                merged.push_back(c);
            } else {
                merged.back().is_interface |= c.is_interface;
                merged.back().is_override_bound |= c.is_override_bound;
                if (c.is_interface) {
                    merged.back().n_left = c.n_left;
                    merged.back().n_right = c.n_right;
                }
            }
        }

        return merged;
    }

    // ===== Collect override segments (in physical coordinates) =====
    std::vector<OverrideSegment1D> collect_override_segments_1d(int axis, real domain_size) {
        std::vector<OverrideSegment1D> segments;

        if (!config.mesh_override_enabled) return segments;

        for (const auto& r : config.override_regions) {
            if (!r.enabled) continue;

            real r_min, r_max, r_dx;
            switch (axis) {
                case 0: r_min = r.x0; r_max = r.x1; r_dx = r.dx; break;
                case 1: r_min = r.y0; r_max = r.y1; r_dx = r.dy; break;
                case 2: r_min = r.z0; r_max = r.z1; r_dx = r.dz; break;
                default: continue;
            }

            r_min = std::max(0.0, std::min(r_min, domain_size));
            r_max = std::max(0.0, std::min(r_max, domain_size));

            if (r_min < r_max) {
                r_dx = std::clamp(r_dx, config.dx_min, config.dx_max);
                segments.push_back({r_min, r_max, r_dx, r.priority});
            }
        }

        return segments;
    }

    // ===== Validate grid =====
    void validate_grid_1d(
        const std::vector<real>& grid_lines,
        const std::vector<GridConstraint1D>& constraints,
        const std::vector<OverrideSegment1D>& overrides,
        int axis,
        const char* mesh_type
    ) {
        const char* axis_name[] = {"X", "Y", "Z"};

        // Check monotonicity
        for (size_t i = 1; i < grid_lines.size(); ++i) {
            if (grid_lines[i] <= grid_lines[i - 1]) {
                std::cerr << "[AutoMesh] WARNING: Non-monotonic " << mesh_type
                          << " grid on " << axis_name[axis] << " at index " << i << "\n";
            }
        }

        // Check grading ratio
        real max_ratio = 1.0;
        for (size_t i = 1; i + 1 < grid_lines.size(); ++i) {
            real dx_prev = grid_lines[i] - grid_lines[i - 1];
            real dx_curr = grid_lines[i + 1] - grid_lines[i];
            real ratio = std::max(dx_curr / dx_prev, dx_prev / dx_curr);
            max_ratio = std::max(max_ratio, ratio);
        }

        std::cout << "  " << axis_name[axis] << " (" << mesh_type << "): "
                  << grid_lines.size() - 1 << " cells, max grading: " << max_ratio;
        if (max_ratio > config.max_grading_ratio * 1.01) {
            std::cout << " (EXCEEDS " << config.max_grading_ratio << ")";
        }
        std::cout << "\n";

        // Check uniform spacing inside override regions
        for (const auto& ov : overrides) {
            std::vector<real> override_spacings;
            for (size_t i = 0; i + 1 < grid_lines.size(); ++i) {
                real center = 0.5 * (grid_lines[i] + grid_lines[i + 1]);
                if (center > ov.start && center < ov.end) {
                    override_spacings.push_back(grid_lines[i + 1] - grid_lines[i]);
                }
            }

            if (!override_spacings.empty()) {
                real min_sp = *std::min_element(override_spacings.begin(), override_spacings.end());
                real max_sp = *std::max_element(override_spacings.begin(), override_spacings.end());
                real variation = (max_sp - min_sp) / min_sp;

                if (variation > 0.01) {
                    std::cout << "    Override [" << ov.start * 1e9 << "-" << ov.end * 1e9
                              << "] nm: variation " << variation * 100 << "% (target: uniform)\n";
                }
            }
        }
    }

    // ===== Helper: Sample refractive index at position =====
    real sample_n_at_position(
        const std::vector<AutoMeshStructureBounds>& structures,
        int axis,
        real pos,
        real background_n
    ) {
        for (const auto& s : structures) {
            real s_min, s_max;
            switch (axis) {
                case 0: s_min = s.x_min; s_max = s.x_max; break;
                case 1: s_min = s.y_min; s_max = s.y_max; break;
                case 2: s_min = s.z_min; s_max = s.z_max; break;
                default: continue;
            }

            if (pos >= s_min && pos < s_max) {
                return s.n_material;
            }
        }
        return background_n;
    }

    // ===== Print physical mesh statistics =====
    void print_physical_mesh_stats(
        const std::vector<real>& x_lines,
        const std::vector<real>& y_lines,
        const std::vector<real>& z_lines
    ) {
        std::cout << "\n[AutoMesh] Physical mesh statistics:\n";

        auto print_1d = [](const std::vector<real>& lines, const char* name) {
            if (lines.size() < 2) return;
            std::vector<real> dx(lines.size() - 1);
            for (size_t i = 0; i < dx.size(); ++i) {
                dx[i] = lines[i + 1] - lines[i];
            }
            real min_dx = *std::min_element(dx.begin(), dx.end());
            real max_dx = *std::max_element(dx.begin(), dx.end());
            real domain = lines.back() - lines.front();

            std::cout << "  " << name << ": " << dx.size() << " cells, domain="
                      << domain * 1e9 << "nm, dx=[" << min_dx * 1e9 << ", "
                      << max_dx * 1e9 << "] nm\n";
        };

        print_1d(x_lines, "X");
        print_1d(y_lines, "Y");
        print_1d(z_lines, "Z");
    }

    // ===== Verify override boundaries exist in grid =====
    void verify_override_boundaries(
        const std::vector<real>& x_lines,
        const std::vector<real>& y_lines,
        const std::vector<real>& z_lines
    ) {
        if (!config.mesh_override_enabled) return;

        std::cout << "\n[AutoMesh] Override boundary verification (physical mesh):\n";

        auto find_nearest = [](const std::vector<real>& lines, real target) -> std::pair<size_t, real> {
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
        };

        for (size_t ri = 0; ri < config.override_regions.size(); ++ri) {
            const auto& r = config.override_regions[ri];
            if (!r.enabled) continue;

            std::cout << "  Region " << ri << " [" << r.x0*1e9 << "-" << r.x1*1e9
                      << "] x [" << r.y0*1e9 << "-" << r.y1*1e9 << "] nm:\n";

            // Check X boundaries
            auto [x0_idx, x0_err] = find_nearest(x_lines, r.x0);
            auto [x1_idx, x1_err] = find_nearest(x_lines, r.x1);
            std::cout << "    X: x0=" << r.x0*1e9 << "nm -> idx " << x0_idx
                      << " (err=" << x0_err*1e9 << "nm" << (x0_err < 1e-12 ? " OK" : " MISS") << "), "
                      << "x1=" << r.x1*1e9 << "nm -> idx " << x1_idx
                      << " (err=" << x1_err*1e9 << "nm" << (x1_err < 1e-12 ? " OK" : " MISS") << ")\n";

            // Check Y boundaries
            auto [y0_idx, y0_err] = find_nearest(y_lines, r.y0);
            auto [y1_idx, y1_err] = find_nearest(y_lines, r.y1);
            std::cout << "    Y: y0=" << r.y0*1e9 << "nm -> idx " << y0_idx
                      << " (err=" << y0_err*1e9 << "nm" << (y0_err < 1e-12 ? " OK" : " MISS") << "), "
                      << "y1=" << r.y1*1e9 << "nm -> idx " << y1_idx
                      << " (err=" << y1_err*1e9 << "nm" << (y1_err < 1e-12 ? " OK" : " MISS") << ")\n";
        }
    }

    // ===== Print final mesh statistics =====
    void print_final_mesh_stats(
        const GridSpacing& grid,
        size_t npml,
        size_t phys_nx, size_t phys_ny, size_t phys_nz
    ) {
        std::cout << "\n[AutoMesh] Final total mesh statistics:\n";

        auto print_axis = [&](const std::vector<real>& dx, const std::vector<real>& bounds,
                             const char* name, size_t phys_n) {
            real min_dx = *std::min_element(dx.begin(), dx.end());
            real max_dx = *std::max_element(dx.begin(), dx.end());
            real total = bounds.back() - bounds.front();

            real max_ratio = 1.0;
            for (size_t i = 1; i < dx.size(); ++i) {
                real ratio = std::max(dx[i] / dx[i-1], dx[i-1] / dx[i]);
                max_ratio = std::max(max_ratio, ratio);
            }

            std::cout << "  " << name << ": " << dx.size() << " total cells ("
                      << phys_n << " physical + " << 2*npml << " PML)\n";
            std::cout << "      Domain: " << total * 1e9 << " nm, dx=["
                      << min_dx * 1e9 << ", " << max_dx * 1e9 << "] nm\n";
            std::cout << "      Max grading: " << max_ratio;
            if (max_ratio > config.max_grading_ratio * 1.01) {
                std::cout << " (EXCEEDS " << config.max_grading_ratio << ")";
            }
            std::cout << "\n";
        };

        print_axis(grid.dx, grid.x_bounds, "X", phys_nx);
        print_axis(grid.dy, grid.y_bounds, "Y", phys_ny);
        print_axis(grid.dz, grid.z_bounds, "Z", phys_nz);

        size_t total_cells = grid.dx.size() * grid.dy.size() * grid.dz.size();
        std::cout << "  Total cells: " << total_cells << "\n";

        // Print physical domain boundaries in total domain
        if (npml > 0 && npml < grid.x_bounds.size()) {
            std::cout << "\n[AutoMesh] Physical domain in total coordinates:\n";
            std::cout << "  X: [" << grid.x_bounds[npml]*1e9 << ", "
                      << grid.x_bounds[grid.x_bounds.size()-1-npml]*1e9 << "] nm\n";
            std::cout << "  Y: [" << grid.y_bounds[npml]*1e9 << ", "
                      << grid.y_bounds[grid.y_bounds.size()-1-npml]*1e9 << "] nm\n";
            std::cout << "  Z: [" << grid.z_bounds[npml]*1e9 << ", "
                      << grid.z_bounds[grid.z_bounds.size()-1-npml]*1e9 << "] nm\n";
        }

        // Print override region uniformity in total mesh
        if (config.mesh_override_enabled && !config.override_regions.empty() && npml > 0) {
            std::cout << "\n[AutoMesh] Override regions in total mesh (shifted by PML):\n";

            for (size_t ri = 0; ri < config.override_regions.size(); ++ri) {
                const auto& r = config.override_regions[ri];
                if (!r.enabled) continue;

                // Override in physical coords -> shift by PML offset for total coords
                real pml_offset_x = grid.x_bounds[npml];
                real pml_offset_y = grid.y_bounds[npml];

                real total_x0 = pml_offset_x + r.x0;
                real total_x1 = pml_offset_x + r.x1;

                // Check uniformity
                std::vector<real> inside_dx;
                for (size_t i = 0; i < grid.dx.size(); ++i) {
                    real center = 0.5 * (grid.x_bounds[i] + grid.x_bounds[i + 1]);
                    if (center > total_x0 && center < total_x1) {
                        inside_dx.push_back(grid.dx[i]);
                    }
                }

                if (!inside_dx.empty()) {
                    real mean_dx = std::accumulate(inside_dx.begin(), inside_dx.end(), 0.0) / inside_dx.size();
                    real min_dx = *std::min_element(inside_dx.begin(), inside_dx.end());
                    real max_dx = *std::max_element(inside_dx.begin(), inside_dx.end());
                    real var_pct = (max_dx - min_dx) / mean_dx * 100;

                    std::cout << "  Region " << ri << " (phys [" << r.x0*1e9 << "-" << r.x1*1e9
                              << "] nm -> total [" << total_x0*1e9 << "-" << total_x1*1e9 << "] nm):\n";
                    std::cout << "    " << inside_dx.size() << " cells, target=" << r.dx*1e9
                              << "nm, actual=" << mean_dx*1e9 << "nm, var=" << var_pct << "%\n";
                }
            }
        }
    }
};

// ========== Helper: Convert StructureBounds to AutoMeshStructureBounds ==========
inline std::vector<AutoMeshStructureBounds> convert_to_auto_mesh_bounds(
    const std::vector<StructureBounds>& bounds
) {
    std::vector<AutoMeshStructureBounds> result;
    result.reserve(bounds.size());
    for (const auto& b : bounds) {
        AutoMeshStructureBounds v2;
        v2.x_min = b.x_min;
        v2.x_max = b.x_max;
        v2.y_min = b.y_min;
        v2.y_max = b.y_max;
        v2.z_min = b.z_min;
        v2.z_max = b.z_max;
        v2.n_material = b.n_material;
        result.push_back(v2);
    }
    return result;
}

// ========== Convenience function ==========
inline GridSpacing generate_auto_mesh(
    real domain_x, real domain_y, real domain_z,
    size_t npml,
    const std::vector<AutoMeshStructureBounds>& structures,
    real lambda_min,
    int mesh_accuracy = 2,
    real background_n = 1.0,
    real max_grading_ratio = 1.41421356237
) {
    AutoMeshConfig config;
    config.lambda_min = lambda_min;
    config.mesh_accuracy = mesh_accuracy;
    config.max_grading_ratio = max_grading_ratio;
    config.dx_min = lambda_min / 100;
    config.dx_max = lambda_min / 4;

    AutoMeshGenerator generator(config);
    return generator.generate(domain_x, domain_y, domain_z, npml, structures, background_n);
}
