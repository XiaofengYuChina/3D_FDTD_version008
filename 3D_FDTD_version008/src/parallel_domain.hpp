// parallel_domain.hpp — Domain Decomposition Parallel Module for 3D FDTD
//
// This module provides domain decomposition parallelization for the FDTD simulation.
// It splits the computational domain into multiple subdomains and handles:
// - Domain partitioning along a chosen axis (X, Y, or Z)
// - Halo (ghost cell) data exchange between adjacent subdomains
// - Parallel field updates with proper synchronization
//
// Compatible with: PML boundaries, non-uniform meshes, two-level system
//
// Usage:
// 1. Set PARALLEL_ENABLED = true in user_config.hpp
// 2. Set NUM_DOMAINS to desired number of parallel regions
// 3. Choose decomposition axis (X=0, Y=1, Z=2)

#pragma once

#include <vector>
#include <cstddef>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <functional>
#include <atomic>

#include "omp_config.hpp"
#include "global_function.hpp"

namespace Parallel {

// ============================================================================
//                          CONFIGURATION
// ============================================================================

// Decomposition axis enum
enum class DecompAxis { X = 0, Y = 1, Z = 2 };

// Parallel configuration structure
struct ParallelConfig {
    bool enabled = false;           // Enable parallel domain decomposition
    int num_domains = 1;            // Number of domains to split into
    DecompAxis axis = DecompAxis::Z; // Axis along which to decompose
    int halo_width = 1;             // Halo (ghost cell) width for data exchange

    // Validate and adjust configuration
    void validate(size_t NxT, size_t NyT, size_t NzT, size_t npml) {
        if (!enabled || num_domains <= 1) {
            enabled = false;
            num_domains = 1;
            return;
        }

        // Get the size along decomposition axis
        size_t axis_size = (axis == DecompAxis::X) ? NxT :
                          (axis == DecompAxis::Y) ? NyT : NzT;

        // Ensure we have enough cells per domain (minimum 4 + 2*halo + 2*npml per domain)
        size_t min_cells_per_domain = 4 + 2 * halo_width;
        size_t max_domains = (axis_size > 2 * npml) ?
                            (axis_size - 2 * npml) / min_cells_per_domain : 1;

        if (num_domains > (int)max_domains) {
            std::cout << "[Parallel] Warning: Reducing num_domains from " << num_domains
                      << " to " << max_domains << " due to grid size constraints\n";
            num_domains = (int)max_domains;
        }

        if (num_domains <= 1) {
            enabled = false;
            num_domains = 1;
        }
    }
};

// ============================================================================
//                          SUBDOMAIN STRUCTURE
// ============================================================================

// Represents a single subdomain in the decomposition
struct Subdomain {
    int id;                         // Domain ID (0 to num_domains-1)

    // Global index range [start, end) along decomposition axis
    size_t global_start;            // Start index in global array (inclusive)
    size_t global_end;              // End index in global array (exclusive)

    // Local dimensions including halo
    size_t local_Nx, local_Ny, local_Nz;
    size_t local_Ntot;              // Total local cells

    // Halo information
    int left_neighbor;              // ID of left neighbor (-1 if none, boundary)
    int right_neighbor;             // ID of right neighbor (-1 if none, boundary)
    size_t halo_width;

    // Offsets for converting local to global indices
    size_t offset_i, offset_j, offset_k;

    // Is this domain at a global boundary (needs PML)?
    bool is_left_boundary;
    bool is_right_boundary;

    // Local field arrays (E, H, J)
    std::vector<real> Ex, Ey, Ez;
    std::vector<real> Hx, Hy, Hz;
    std::vector<real> Jx, Jy, Jz;

    // Local coefficient arrays (copied from global for this subdomain)
    std::vector<real> aEx, bEx, aEy, bEy, aEz, bEz;
    std::vector<real> aHx, bHx, aHy, bHy, aHz, bHz;

    // Local grid spacing (subset of global spacing)
    std::vector<real> dx_local, dy_local, dz_local;

    // Allocate local arrays
    void allocate() {
        local_Ntot = local_Nx * local_Ny * local_Nz;

        Ex.assign(local_Ntot, 0);
        Ey.assign(local_Ntot, 0);
        Ez.assign(local_Ntot, 0);
        Hx.assign(local_Ntot, 0);
        Hy.assign(local_Ntot, 0);
        Hz.assign(local_Ntot, 0);
        Jx.assign(local_Ntot, 0);
        Jy.assign(local_Ntot, 0);
        Jz.assign(local_Ntot, 0);

        aEx.assign(local_Ntot, 1);
        bEx.assign(local_Ntot, 0);
        aEy.assign(local_Ntot, 1);
        bEy.assign(local_Ntot, 0);
        aEz.assign(local_Ntot, 1);
        bEz.assign(local_Ntot, 0);
        aHx.assign(local_Ntot, 1);
        bHx.assign(local_Ntot, 0);
        aHy.assign(local_Ntot, 1);
        bHy.assign(local_Ntot, 0);
        aHz.assign(local_Ntot, 1);
        bHz.assign(local_Ntot, 0);
    }

    // Convert local index to 1D index
    inline size_t local_idx(size_t i, size_t j, size_t k) const {
        return (i * local_Ny + j) * local_Nz + k;
    }

    // Convert local index to global index
    inline void local_to_global(size_t li, size_t lj, size_t lk,
                                size_t& gi, size_t& gj, size_t& gk) const {
        gi = li + offset_i;
        gj = lj + offset_j;
        gk = lk + offset_k;
    }
};

// ============================================================================
//                     DOMAIN DECOMPOSITION MANAGER
// ============================================================================

class DomainDecomposition {
public:
    ParallelConfig config;
    size_t NxT, NyT, NzT, Ntot;
    size_t npml;

    std::vector<Subdomain> subdomains;

    // Synchronization primitives
    std::unique_ptr<std::mutex[]> exchange_mutexes;
    int num_exchange_mutexes = 0;
    std::atomic<int> barrier_count{0};
    std::mutex barrier_mutex;
    std::condition_variable barrier_cv;

    DomainDecomposition() = default;

    // Initialize domain decomposition
    void initialize(const ParallelConfig& cfg,
                   size_t Nx, size_t Ny, size_t Nz,
                   size_t pml_thickness,
                   const GridSpacing& grid_spacing,
                   const MaterialGrids& mats) {

        config = cfg;
        NxT = Nx;
        NyT = Ny;
        NzT = Nz;
        Ntot = NxT * NyT * NzT;
        npml = pml_thickness;

        // Validate configuration
        config.validate(NxT, NyT, NzT, npml);

        if (!config.enabled) {
            std::cout << "[Parallel] Domain decomposition disabled (single domain mode)\n";
            return;
        }

        std::cout << "[Parallel] Initializing domain decomposition:\n";
        std::cout << "  Number of domains: " << config.num_domains << "\n";
        std::cout << "  Decomposition axis: " << (config.axis == DecompAxis::X ? "X" :
                                                  config.axis == DecompAxis::Y ? "Y" : "Z") << "\n";
        std::cout << "  Halo width: " << config.halo_width << "\n";

        // Create subdomains
        create_subdomains(grid_spacing, mats);

        // Initialize synchronization
        num_exchange_mutexes = config.num_domains;
        exchange_mutexes = std::make_unique<std::mutex[]>(num_exchange_mutexes);

        std::cout << "[Parallel] Domain decomposition initialized successfully\n";
    }

    // Check if parallel mode is active
    bool is_parallel() const { return config.enabled && config.num_domains > 1; }

    // Get number of domains
    int num_domains() const { return config.num_domains; }

    // Scatter global fields to subdomains
    void scatter_fields(const std::vector<real>& global_Ex,
                       const std::vector<real>& global_Ey,
                       const std::vector<real>& global_Ez,
                       const std::vector<real>& global_Hx,
                       const std::vector<real>& global_Hy,
                       const std::vector<real>& global_Hz) {
        if (!is_parallel()) return;

        for (auto& dom : subdomains) {
            scatter_to_subdomain(dom, global_Ex, global_Ey, global_Ez,
                                global_Hx, global_Hy, global_Hz);
        }
    }

    // Gather subdomain fields back to global arrays
    void gather_fields(std::vector<real>& global_Ex,
                      std::vector<real>& global_Ey,
                      std::vector<real>& global_Ez,
                      std::vector<real>& global_Hx,
                      std::vector<real>& global_Hy,
                      std::vector<real>& global_Hz) {
        if (!is_parallel()) return;

        for (const auto& dom : subdomains) {
            gather_from_subdomain(dom, global_Ex, global_Ey, global_Ez,
                                 global_Hx, global_Hy, global_Hz);
        }
    }

    // Exchange halo data between adjacent subdomains
    // Call after H-field update and before E-field update
    void exchange_halos_H() {
        if (!is_parallel()) return;
        exchange_halos_internal(true);
    }

    // Call after E-field update and before H-field update
    void exchange_halos_E() {
        if (!is_parallel()) return;
        exchange_halos_internal(false);
    }

    // Parallel barrier - wait for all domains to complete
    void barrier() {
        if (!is_parallel()) return;

        std::unique_lock<std::mutex> lock(barrier_mutex);
        int expected = barrier_count.fetch_add(1) + 1;

        if (expected == config.num_domains) {
            barrier_count.store(0);
            barrier_cv.notify_all();
        } else {
            barrier_cv.wait(lock, [this]() { return barrier_count.load() == 0; });
        }
    }

private:
    // Create subdomain partitions
    void create_subdomains(const GridSpacing& grid_spacing,
                          const MaterialGrids& mats) {
        subdomains.clear();
        subdomains.resize(config.num_domains);

        // Get axis size
        size_t axis_size = (config.axis == DecompAxis::X) ? NxT :
                          (config.axis == DecompAxis::Y) ? NyT : NzT;

        // Calculate interior size (excluding PML)
        size_t interior_size = axis_size - 2 * npml;

        // Distribute interior cells among domains
        size_t base_size = interior_size / config.num_domains;
        size_t remainder = interior_size % config.num_domains;

        size_t current_start = 0;

        for (int d = 0; d < config.num_domains; ++d) {
            Subdomain& dom = subdomains[d];
            dom.id = d;
            dom.halo_width = config.halo_width;

            // Calculate this domain's portion of interior cells
            size_t domain_interior_cells = base_size + (d < (int)remainder ? 1 : 0);

            // Global start and end indices
            if (d == 0) {
                // First domain includes left PML
                dom.global_start = 0;
                dom.is_left_boundary = true;
            } else {
                // Start after previous domain (with overlap for halo)
                dom.global_start = npml + current_start - config.halo_width;
                dom.is_left_boundary = false;
            }

            current_start += domain_interior_cells;

            if (d == config.num_domains - 1) {
                // Last domain includes right PML
                dom.global_end = axis_size;
                dom.is_right_boundary = true;
            } else {
                // End with some extra for halo
                dom.global_end = npml + current_start + config.halo_width;
                dom.is_right_boundary = false;
            }

            // Set neighbors
            dom.left_neighbor = (d > 0) ? d - 1 : -1;
            dom.right_neighbor = (d < config.num_domains - 1) ? d + 1 : -1;

            // Calculate local dimensions based on decomposition axis
            switch (config.axis) {
                case DecompAxis::X:
                    dom.local_Nx = dom.global_end - dom.global_start;
                    dom.local_Ny = NyT;
                    dom.local_Nz = NzT;
                    dom.offset_i = dom.global_start;
                    dom.offset_j = 0;
                    dom.offset_k = 0;
                    break;
                case DecompAxis::Y:
                    dom.local_Nx = NxT;
                    dom.local_Ny = dom.global_end - dom.global_start;
                    dom.local_Nz = NzT;
                    dom.offset_i = 0;
                    dom.offset_j = dom.global_start;
                    dom.offset_k = 0;
                    break;
                case DecompAxis::Z:
                    dom.local_Nx = NxT;
                    dom.local_Ny = NyT;
                    dom.local_Nz = dom.global_end - dom.global_start;
                    dom.offset_i = 0;
                    dom.offset_j = 0;
                    dom.offset_k = dom.global_start;
                    break;
            }

            // Allocate local arrays
            dom.allocate();

            // Copy grid spacing for this subdomain
            copy_grid_spacing(dom, grid_spacing);

            // Copy material coefficients for this subdomain
            copy_material_coeffs(dom, mats);

            std::cout << "  Domain " << d << ": global [" << dom.global_start << ", "
                      << dom.global_end << "), local size " << dom.local_Nx << "x"
                      << dom.local_Ny << "x" << dom.local_Nz << "\n";
        }
    }

    // Copy grid spacing to subdomain
    void copy_grid_spacing(Subdomain& dom, const GridSpacing& grid_spacing) {
        // Copy relevant portion of spacing arrays
        switch (config.axis) {
            case DecompAxis::X:
                dom.dx_local.assign(grid_spacing.dx.begin() + dom.global_start,
                                   grid_spacing.dx.begin() + dom.global_end);
                dom.dy_local = grid_spacing.dy;
                dom.dz_local = grid_spacing.dz;
                break;
            case DecompAxis::Y:
                dom.dx_local = grid_spacing.dx;
                dom.dy_local.assign(grid_spacing.dy.begin() + dom.global_start,
                                   grid_spacing.dy.begin() + dom.global_end);
                dom.dz_local = grid_spacing.dz;
                break;
            case DecompAxis::Z:
                dom.dx_local = grid_spacing.dx;
                dom.dy_local = grid_spacing.dy;
                dom.dz_local.assign(grid_spacing.dz.begin() + dom.global_start,
                                   grid_spacing.dz.begin() + dom.global_end);
                break;
        }
    }

    // Copy material coefficients to subdomain
    void copy_material_coeffs(Subdomain& dom, const MaterialGrids& mats) {
        for (size_t li = 0; li < dom.local_Nx; ++li) {
            for (size_t lj = 0; lj < dom.local_Ny; ++lj) {
                for (size_t lk = 0; lk < dom.local_Nz; ++lk) {
                    size_t gi, gj, gk;
                    dom.local_to_global(li, lj, lk, gi, gj, gk);

                    size_t lid = dom.local_idx(li, lj, lk);
                    size_t gid = idx3(gi, gj, gk, NyT, NzT);

                    dom.aEx[lid] = mats.aEx[gid];
                    dom.bEx[lid] = mats.bEx[gid];
                    dom.aEy[lid] = mats.aEy[gid];
                    dom.bEy[lid] = mats.bEy[gid];
                    dom.aEz[lid] = mats.aEz[gid];
                    dom.bEz[lid] = mats.bEz[gid];
                    dom.aHx[lid] = mats.aHx[gid];
                    dom.bHx[lid] = mats.bHx[gid];
                    dom.aHy[lid] = mats.aHy[gid];
                    dom.bHy[lid] = mats.bHy[gid];
                    dom.aHz[lid] = mats.aHz[gid];
                    dom.bHz[lid] = mats.bHz[gid];
                }
            }
        }
    }

    // Scatter global fields to a single subdomain
    void scatter_to_subdomain(Subdomain& dom,
                             const std::vector<real>& global_Ex,
                             const std::vector<real>& global_Ey,
                             const std::vector<real>& global_Ez,
                             const std::vector<real>& global_Hx,
                             const std::vector<real>& global_Hy,
                             const std::vector<real>& global_Hz) {

        for (size_t li = 0; li < dom.local_Nx; ++li) {
            for (size_t lj = 0; lj < dom.local_Ny; ++lj) {
                for (size_t lk = 0; lk < dom.local_Nz; ++lk) {
                    size_t gi, gj, gk;
                    dom.local_to_global(li, lj, lk, gi, gj, gk);

                    size_t lid = dom.local_idx(li, lj, lk);
                    size_t gid = idx3(gi, gj, gk, NyT, NzT);

                    dom.Ex[lid] = global_Ex[gid];
                    dom.Ey[lid] = global_Ey[gid];
                    dom.Ez[lid] = global_Ez[gid];
                    dom.Hx[lid] = global_Hx[gid];
                    dom.Hy[lid] = global_Hy[gid];
                    dom.Hz[lid] = global_Hz[gid];
                }
            }
        }
    }

    // Gather fields from a single subdomain (only interior, not halo)
    void gather_from_subdomain(const Subdomain& dom,
                              std::vector<real>& global_Ex,
                              std::vector<real>& global_Ey,
                              std::vector<real>& global_Ez,
                              std::vector<real>& global_Hx,
                              std::vector<real>& global_Hy,
                              std::vector<real>& global_Hz) {

        // Determine interior range (exclude halo cells)
        size_t i_start = dom.is_left_boundary ? 0 : config.halo_width;
        size_t i_end = dom.local_Nx - (dom.is_right_boundary ? 0 : config.halo_width);
        size_t j_start = 0, j_end = dom.local_Ny;
        size_t k_start = 0, k_end = dom.local_Nz;

        // Adjust based on decomposition axis
        switch (config.axis) {
            case DecompAxis::X:
                // Already set above
                break;
            case DecompAxis::Y:
                i_start = 0; i_end = dom.local_Nx;
                j_start = dom.is_left_boundary ? 0 : config.halo_width;
                j_end = dom.local_Ny - (dom.is_right_boundary ? 0 : config.halo_width);
                break;
            case DecompAxis::Z:
                i_start = 0; i_end = dom.local_Nx;
                k_start = dom.is_left_boundary ? 0 : config.halo_width;
                k_end = dom.local_Nz - (dom.is_right_boundary ? 0 : config.halo_width);
                break;
        }

        for (size_t li = i_start; li < i_end; ++li) {
            for (size_t lj = j_start; lj < j_end; ++lj) {
                for (size_t lk = k_start; lk < k_end; ++lk) {
                    size_t gi, gj, gk;
                    dom.local_to_global(li, lj, lk, gi, gj, gk);

                    size_t lid = dom.local_idx(li, lj, lk);
                    size_t gid = idx3(gi, gj, gk, NyT, NzT);

                    global_Ex[gid] = dom.Ex[lid];
                    global_Ey[gid] = dom.Ey[lid];
                    global_Ez[gid] = dom.Ez[lid];
                    global_Hx[gid] = dom.Hx[lid];
                    global_Hy[gid] = dom.Hy[lid];
                    global_Hz[gid] = dom.Hz[lid];
                }
            }
        }
    }

    // Internal halo exchange implementation
    void exchange_halos_internal(bool after_H) {
        // For each pair of adjacent domains, exchange halo data
        for (int d = 0; d < config.num_domains - 1; ++d) {
            Subdomain& left = subdomains[d];
            Subdomain& right = subdomains[d + 1];

            // Exchange data based on decomposition axis
            switch (config.axis) {
                case DecompAxis::X:
                    exchange_halos_X(left, right, after_H);
                    break;
                case DecompAxis::Y:
                    exchange_halos_Y(left, right, after_H);
                    break;
                case DecompAxis::Z:
                    exchange_halos_Z(left, right, after_H);
                    break;
            }
        }
    }

    // Exchange halos along X axis
    void exchange_halos_X(Subdomain& left, Subdomain& right, bool after_H) {
        size_t hw = config.halo_width;

        // Copy from left's right edge to right's left halo
        // Copy from right's left edge to left's right halo
        for (size_t j = 0; j < NyT; ++j) {
            for (size_t k = 0; k < NzT; ++k) {
                for (size_t h = 0; h < hw; ++h) {
                    // Left domain's right interior -> Right domain's left halo
                    size_t left_src = left.local_idx(left.local_Nx - 2*hw + h, j, k);
                    size_t right_dst = right.local_idx(h, j, k);

                    // Right domain's left interior -> Left domain's right halo
                    size_t right_src = right.local_idx(hw + h, j, k);
                    size_t left_dst = left.local_idx(left.local_Nx - hw + h, j, k);

                    if (after_H) {
                        // After H update: exchange H fields
                        right.Hx[right_dst] = left.Hx[left_src];
                        right.Hy[right_dst] = left.Hy[left_src];
                        right.Hz[right_dst] = left.Hz[left_src];

                        left.Hx[left_dst] = right.Hx[right_src];
                        left.Hy[left_dst] = right.Hy[right_src];
                        left.Hz[left_dst] = right.Hz[right_src];
                    } else {
                        // After E update: exchange E fields
                        right.Ex[right_dst] = left.Ex[left_src];
                        right.Ey[right_dst] = left.Ey[left_src];
                        right.Ez[right_dst] = left.Ez[left_src];

                        left.Ex[left_dst] = right.Ex[right_src];
                        left.Ey[left_dst] = right.Ey[right_src];
                        left.Ez[left_dst] = right.Ez[right_src];
                    }
                }
            }
        }
    }

    // Exchange halos along Y axis
    void exchange_halos_Y(Subdomain& left, Subdomain& right, bool after_H) {
        size_t hw = config.halo_width;

        for (size_t i = 0; i < NxT; ++i) {
            for (size_t k = 0; k < NzT; ++k) {
                for (size_t h = 0; h < hw; ++h) {
                    size_t left_src = left.local_idx(i, left.local_Ny - 2*hw + h, k);
                    size_t right_dst = right.local_idx(i, h, k);

                    size_t right_src = right.local_idx(i, hw + h, k);
                    size_t left_dst = left.local_idx(i, left.local_Ny - hw + h, k);

                    if (after_H) {
                        right.Hx[right_dst] = left.Hx[left_src];
                        right.Hy[right_dst] = left.Hy[left_src];
                        right.Hz[right_dst] = left.Hz[left_src];

                        left.Hx[left_dst] = right.Hx[right_src];
                        left.Hy[left_dst] = right.Hy[right_src];
                        left.Hz[left_dst] = right.Hz[right_src];
                    } else {
                        right.Ex[right_dst] = left.Ex[left_src];
                        right.Ey[right_dst] = left.Ey[left_src];
                        right.Ez[right_dst] = left.Ez[left_src];

                        left.Ex[left_dst] = right.Ex[right_src];
                        left.Ey[left_dst] = right.Ey[right_src];
                        left.Ez[left_dst] = right.Ez[right_src];
                    }
                }
            }
        }
    }

    // Exchange halos along Z axis
    void exchange_halos_Z(Subdomain& left, Subdomain& right, bool after_H) {
        size_t hw = config.halo_width;

        for (size_t i = 0; i < NxT; ++i) {
            for (size_t j = 0; j < NyT; ++j) {
                for (size_t h = 0; h < hw; ++h) {
                    size_t left_src = left.local_idx(i, j, left.local_Nz - 2*hw + h);
                    size_t right_dst = right.local_idx(i, j, h);

                    size_t right_src = right.local_idx(i, j, hw + h);
                    size_t left_dst = left.local_idx(i, j, left.local_Nz - hw + h);

                    if (after_H) {
                        right.Hx[right_dst] = left.Hx[left_src];
                        right.Hy[right_dst] = left.Hy[left_src];
                        right.Hz[right_dst] = left.Hz[left_src];

                        left.Hx[left_dst] = right.Hx[right_src];
                        left.Hy[left_dst] = right.Hy[right_src];
                        left.Hz[left_dst] = right.Hz[right_src];
                    } else {
                        right.Ex[right_dst] = left.Ex[left_src];
                        right.Ey[right_dst] = left.Ey[left_src];
                        right.Ez[right_dst] = left.Ez[left_src];

                        left.Ex[left_dst] = right.Ex[right_src];
                        left.Ey[left_dst] = right.Ey[right_src];
                        left.Ez[left_dst] = right.Ez[right_src];
                    }
                }
            }
        }
    }
};

// ============================================================================
//                     PARALLEL FDTD STEPPER FUNCTIONS
// ============================================================================

// Update H field for a single subdomain
template<typename Real>
inline void fdtd_update_H_subdomain(Subdomain& dom) {
    using namespace fdtd_math;

    const size_t NxL = dom.local_Nx;
    const size_t NyL = dom.local_Ny;
    const size_t NzL = dom.local_Nz;

    const size_t sI = NyL * NzL;
    const size_t sJ = NzL;
    const size_t sK = 1;

    const Real* __restrict dx_arr = dom.dx_local.data();
    const Real* __restrict dy_arr = dom.dy_local.data();
    const Real* __restrict dz_arr = dom.dz_local.data();

    Real* __restrict Ex = dom.Ex.data();
    Real* __restrict Ey = dom.Ey.data();
    Real* __restrict Ez = dom.Ez.data();
    Real* __restrict Hx = dom.Hx.data();
    Real* __restrict Hy = dom.Hy.data();
    Real* __restrict Hz = dom.Hz.data();

    const Real* __restrict aHx = dom.aHx.data();
    const Real* __restrict bHx = dom.bHx.data();
    const Real* __restrict aHy = dom.aHy.data();
    const Real* __restrict bHy = dom.bHy.data();
    const Real* __restrict aHz = dom.aHz.data();
    const Real* __restrict bHz = dom.bHz.data();

#if FDTD_OMP_ENABLED
#pragma omp parallel for
#endif
    for (size_t i = 0; i < NxL - 1; ++i) {
        for (size_t j = 0; j < NyL - 1; ++j) {
            for (size_t k = 0; k < NzL - 1; ++k) {
                size_t id = dom.local_idx(i, j, k);

                const Real inv_dx = Real(1) / dx_arr[i];
                const Real inv_dy = Real(1) / dy_arr[j];
                const Real inv_dz = Real(1) / dz_arr[k];

                Real curlEx = diff_y(Ez, id, sJ, inv_dy) - diff_z(Ey, id, sK, inv_dz);
                Real curlEy = diff_z(Ex, id, sK, inv_dz) - diff_x(Ez, id, sI, inv_dx);
                Real curlEz = diff_x(Ey, id, sI, inv_dx) - diff_y(Ex, id, sJ, inv_dy);

                Hx[id] = aHx[id] * Hx[id] - bHx[id] * curlEx;
                Hy[id] = aHy[id] * Hy[id] - bHy[id] * curlEy;
                Hz[id] = aHz[id] * Hz[id] - bHz[id] * curlEz;
            }
        }
    }
}

// Update E field for a single subdomain
template<typename Real>
inline void fdtd_update_E_subdomain(Subdomain& dom) {
    using namespace fdtd_math;

    const size_t NxL = dom.local_Nx;
    const size_t NyL = dom.local_Ny;
    const size_t NzL = dom.local_Nz;

    const size_t sI = NyL * NzL;
    const size_t sJ = NzL;
    const size_t sK = 1;

    const Real* __restrict dx_arr = dom.dx_local.data();
    const Real* __restrict dy_arr = dom.dy_local.data();
    const Real* __restrict dz_arr = dom.dz_local.data();

    Real* __restrict Ex = dom.Ex.data();
    Real* __restrict Ey = dom.Ey.data();
    Real* __restrict Ez = dom.Ez.data();
    const Real* __restrict Hx = dom.Hx.data();
    const Real* __restrict Hy = dom.Hy.data();
    const Real* __restrict Hz = dom.Hz.data();
    const Real* __restrict Jx = dom.Jx.data();
    const Real* __restrict Jy = dom.Jy.data();
    const Real* __restrict Jz = dom.Jz.data();

    const Real* __restrict aEx = dom.aEx.data();
    const Real* __restrict bEx = dom.bEx.data();
    const Real* __restrict aEy = dom.aEy.data();
    const Real* __restrict bEy = dom.bEy.data();
    const Real* __restrict aEz = dom.aEz.data();
    const Real* __restrict bEz = dom.bEz.data();

#if FDTD_OMP_ENABLED
#pragma omp parallel for
#endif
    for (size_t i = 1; i < NxL; ++i) {
        for (size_t j = 1; j < NyL; ++j) {
            for (size_t k = 1; k < NzL; ++k) {
                size_t id = dom.local_idx(i, j, k);

                const Real inv_dx = Real(1) / dx_arr[i];
                const Real inv_dy = Real(1) / dy_arr[j];
                const Real inv_dz = Real(1) / dz_arr[k];

                Real curlHx = diff_ym(Hz, id, sJ, inv_dy) - diff_zm(Hy, id, sK, inv_dz);
                Real curlHy = diff_zm(Hx, id, sK, inv_dz) - diff_xm(Hz, id, sI, inv_dx);
                Real curlHz = diff_xm(Hy, id, sI, inv_dx) - diff_ym(Hx, id, sJ, inv_dy);

                Ex[id] = aEx[id] * Ex[id] + bEx[id] * (curlHx - Jx[id]);
                Ey[id] = aEy[id] * Ey[id] + bEy[id] * (curlHy - Jy[id]);
                Ez[id] = aEz[id] * Ez[id] + bEz[id] * (curlHz - Jz[id]);
            }
        }
    }
}

// ============================================================================
//                     PARALLEL EXECUTION HELPER
// ============================================================================

// Execute a function in parallel across all subdomains
template<typename Func>
void parallel_for_domains(DomainDecomposition& decomp, Func&& func) {
    if (!decomp.is_parallel()) {
        return;
    }

#if FDTD_OMP_ENABLED
    #pragma omp parallel for
#endif
    for (int d = 0; d < decomp.num_domains(); ++d) {
        func(decomp.subdomains[d]);
    }
}

// Execute FDTD time step in parallel
// This is the main parallel stepping function
inline void parallel_fdtd_step_H(DomainDecomposition& decomp) {
    if (!decomp.is_parallel()) return;

    // Update H fields in all subdomains
    parallel_for_domains(decomp, [](Subdomain& dom) {
        fdtd_update_H_subdomain<real>(dom);
    });

    // Exchange H-field halos
    decomp.exchange_halos_H();
}

inline void parallel_fdtd_step_E(DomainDecomposition& decomp) {
    if (!decomp.is_parallel()) return;

    // Update E fields in all subdomains
    parallel_for_domains(decomp, [](Subdomain& dom) {
        fdtd_update_E_subdomain<real>(dom);
    });

    // Exchange E-field halos
    decomp.exchange_halos_E();
}

// Clear J fields in all subdomains
inline void parallel_clear_J(DomainDecomposition& decomp) {
    if (!decomp.is_parallel()) return;

    parallel_for_domains(decomp, [](Subdomain& dom) {
        std::fill(dom.Jx.begin(), dom.Jx.end(), 0.0);
        std::fill(dom.Jy.begin(), dom.Jy.end(), 0.0);
        std::fill(dom.Jz.begin(), dom.Jz.end(), 0.0);
    });
}

// Inject source current into appropriate subdomain
inline void parallel_inject_source(DomainDecomposition& decomp,
                                   size_t global_i, size_t global_j, size_t global_k,
                                   real Jz_value) {
    if (!decomp.is_parallel()) return;

    // Find which subdomain contains this source point
    for (auto& dom : decomp.subdomains) {
        size_t gi, gj, gk;
        bool found = false;

        // Check if source is in this subdomain
        switch (decomp.config.axis) {
            case DecompAxis::X:
                if (global_i >= dom.offset_i && global_i < dom.offset_i + dom.local_Nx) {
                    size_t li = global_i - dom.offset_i;
                    size_t lid = dom.local_idx(li, global_j, global_k);
                    dom.Jz[lid] += Jz_value;
                    found = true;
                }
                break;
            case DecompAxis::Y:
                if (global_j >= dom.offset_j && global_j < dom.offset_j + dom.local_Ny) {
                    size_t lj = global_j - dom.offset_j;
                    size_t lid = dom.local_idx(global_i, lj, global_k);
                    dom.Jz[lid] += Jz_value;
                    found = true;
                }
                break;
            case DecompAxis::Z:
                if (global_k >= dom.offset_k && global_k < dom.offset_k + dom.local_Nz) {
                    size_t lk = global_k - dom.offset_k;
                    size_t lid = dom.local_idx(global_i, global_j, lk);
                    dom.Jz[lid] += Jz_value;
                    found = true;
                }
                break;
        }

        if (found) break;
    }
}

} // namespace Parallel
