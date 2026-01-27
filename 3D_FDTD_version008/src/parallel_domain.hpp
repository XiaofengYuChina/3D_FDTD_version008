// parallel_domain.hpp - Domain Decomposition Parallel Module for 3D FDTD
//
// Provides domain decomposition parallelization for the FDTD simulation.
// Splits the computational domain into multiple subdomains and handles:
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

enum class DecompAxis { X = 0, Y = 1, Z = 2 };

struct ParallelConfig {
    bool enabled = false;
    int num_domains = 1;
    DecompAxis axis = DecompAxis::Z;
    int halo_width = 1;

    void validate(size_t NxT, size_t NyT, size_t NzT, size_t npml) {
        if (!enabled || num_domains <= 1) {
            enabled = false;
            num_domains = 1;
            return;
        }

        size_t axis_size = (axis == DecompAxis::X) ? NxT :
                          (axis == DecompAxis::Y) ? NyT : NzT;

        // Minimum 4 + 2*halo cells per domain in interior region
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

struct Subdomain {
    int id;
    size_t global_start, global_end;
    size_t local_Nx, local_Ny, local_Nz, local_Ntot;
    int left_neighbor, right_neighbor;
    size_t halo_width;
    size_t offset_i, offset_j, offset_k;
    bool is_left_boundary, is_right_boundary;

    std::vector<real> Ex, Ey, Ez;
    std::vector<real> Hx, Hy, Hz;
    std::vector<real> Jx, Jy, Jz;
    std::vector<real> aEx, bEx, aEy, bEy, aEz, bEz;
    std::vector<real> aHx, bHx, aHy, bHy, aHz, bHz;
    std::vector<real> dx_local, dy_local, dz_local;

    void allocate() {
        local_Ntot = local_Nx * local_Ny * local_Nz;

        Ex.assign(local_Ntot, 0); Ey.assign(local_Ntot, 0); Ez.assign(local_Ntot, 0);
        Hx.assign(local_Ntot, 0); Hy.assign(local_Ntot, 0); Hz.assign(local_Ntot, 0);
        Jx.assign(local_Ntot, 0); Jy.assign(local_Ntot, 0); Jz.assign(local_Ntot, 0);

        aEx.assign(local_Ntot, 1); bEx.assign(local_Ntot, 0);
        aEy.assign(local_Ntot, 1); bEy.assign(local_Ntot, 0);
        aEz.assign(local_Ntot, 1); bEz.assign(local_Ntot, 0);
        aHx.assign(local_Ntot, 1); bHx.assign(local_Ntot, 0);
        aHy.assign(local_Ntot, 1); bHy.assign(local_Ntot, 0);
        aHz.assign(local_Ntot, 1); bHz.assign(local_Ntot, 0);
    }

    inline size_t local_idx(size_t i, size_t j, size_t k) const {
        return (i * local_Ny + j) * local_Nz + k;
    }

    inline void local_to_global(size_t li, size_t lj, size_t lk,
                                size_t& gi, size_t& gj, size_t& gk) const {
        gi = li + offset_i;
        gj = lj + offset_j;
        gk = lk + offset_k;
    }
};

class DomainDecomposition {
public:
    ParallelConfig config;
    size_t NxT, NyT, NzT, Ntot;
    size_t npml;
    std::vector<Subdomain> subdomains;

    std::unique_ptr<std::mutex[]> exchange_mutexes;
    int num_exchange_mutexes = 0;
    std::atomic<int> barrier_count{0};
    std::mutex barrier_mutex;
    std::condition_variable barrier_cv;

    DomainDecomposition() = default;

    void initialize(const ParallelConfig& cfg,
                   size_t Nx, size_t Ny, size_t Nz,
                   size_t pml_thickness,
                   const GridSpacing& grid_spacing,
                   const MaterialGrids& mats) {
        config = cfg;
        NxT = Nx; NyT = Ny; NzT = Nz;
        Ntot = NxT * NyT * NzT;
        npml = pml_thickness;

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

        create_subdomains(grid_spacing, mats);

        num_exchange_mutexes = config.num_domains;
        exchange_mutexes = std::make_unique<std::mutex[]>(num_exchange_mutexes);

        std::cout << "[Parallel] Domain decomposition initialized successfully\n";
    }

    bool is_parallel() const { return config.enabled && config.num_domains > 1; }
    int num_domains() const { return config.num_domains; }

    void scatter_fields(const std::vector<real>& global_Ex,
                       const std::vector<real>& global_Ey,
                       const std::vector<real>& global_Ez,
                       const std::vector<real>& global_Hx,
                       const std::vector<real>& global_Hy,
                       const std::vector<real>& global_Hz) {
        if (!is_parallel()) return;
        for (auto& dom : subdomains)
            scatter_to_subdomain(dom, global_Ex, global_Ey, global_Ez,
                                global_Hx, global_Hy, global_Hz);
    }

    void gather_fields(std::vector<real>& global_Ex,
                      std::vector<real>& global_Ey,
                      std::vector<real>& global_Ez,
                      std::vector<real>& global_Hx,
                      std::vector<real>& global_Hy,
                      std::vector<real>& global_Hz) {
        if (!is_parallel()) return;
        for (const auto& dom : subdomains)
            gather_from_subdomain(dom, global_Ex, global_Ey, global_Ez,
                                 global_Hx, global_Hy, global_Hz);
    }

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
    void create_subdomains(const GridSpacing& grid_spacing, const MaterialGrids& mats) {
        subdomains.clear();
        subdomains.resize(config.num_domains);

        size_t axis_size = (config.axis == DecompAxis::X) ? NxT :
                          (config.axis == DecompAxis::Y) ? NyT : NzT;

        size_t interior_size = axis_size - 2 * npml;
        size_t base_size = interior_size / config.num_domains;
        size_t remainder = interior_size % config.num_domains;
        size_t current_start = 0;

        for (int d = 0; d < config.num_domains; ++d) {
            Subdomain& dom = subdomains[d];
            dom.id = d;
            dom.halo_width = config.halo_width;

            size_t domain_interior_cells = base_size + (d < (int)remainder ? 1 : 0);

            if (d == 0) {
                dom.global_start = 0;
                dom.is_left_boundary = true;
            } else {
                dom.global_start = npml + current_start - config.halo_width;
                dom.is_left_boundary = false;
            }

            current_start += domain_interior_cells;

            if (d == config.num_domains - 1) {
                dom.global_end = axis_size;
                dom.is_right_boundary = true;
            } else {
                dom.global_end = npml + current_start + config.halo_width;
                dom.is_right_boundary = false;
            }

            dom.left_neighbor = (d > 0) ? d - 1 : -1;
            dom.right_neighbor = (d < config.num_domains - 1) ? d + 1 : -1;

            switch (config.axis) {
                case DecompAxis::X:
                    dom.local_Nx = dom.global_end - dom.global_start;
                    dom.local_Ny = NyT; dom.local_Nz = NzT;
                    dom.offset_i = dom.global_start;
                    dom.offset_j = 0; dom.offset_k = 0;
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
                    dom.local_Nx = NxT; dom.local_Ny = NyT;
                    dom.local_Nz = dom.global_end - dom.global_start;
                    dom.offset_i = 0; dom.offset_j = 0;
                    dom.offset_k = dom.global_start;
                    break;
            }

            dom.allocate();
            copy_grid_spacing(dom, grid_spacing);
            copy_material_coeffs(dom, mats);

            std::cout << "  Domain " << d << ": global [" << dom.global_start << ", "
                      << dom.global_end << "), local size " << dom.local_Nx << "x"
                      << dom.local_Ny << "x" << dom.local_Nz << "\n";
        }
    }

    void copy_grid_spacing(Subdomain& dom, const GridSpacing& grid_spacing) {
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

    void copy_material_coeffs(Subdomain& dom, const MaterialGrids& mats) {
        for (size_t li = 0; li < dom.local_Nx; ++li) {
            for (size_t lj = 0; lj < dom.local_Ny; ++lj) {
                for (size_t lk = 0; lk < dom.local_Nz; ++lk) {
                    size_t gi, gj, gk;
                    dom.local_to_global(li, lj, lk, gi, gj, gk);

                    size_t lid = dom.local_idx(li, lj, lk);
                    size_t gid = idx3(gi, gj, gk, NyT, NzT);

                    dom.aEx[lid] = mats.aEx[gid]; dom.bEx[lid] = mats.bEx[gid];
                    dom.aEy[lid] = mats.aEy[gid]; dom.bEy[lid] = mats.bEy[gid];
                    dom.aEz[lid] = mats.aEz[gid]; dom.bEz[lid] = mats.bEz[gid];
                    dom.aHx[lid] = mats.aHx[gid]; dom.bHx[lid] = mats.bHx[gid];
                    dom.aHy[lid] = mats.aHy[gid]; dom.bHy[lid] = mats.bHy[gid];
                    dom.aHz[lid] = mats.aHz[gid]; dom.bHz[lid] = mats.bHz[gid];
                }
            }
        }
    }

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

                    dom.Ex[lid] = global_Ex[gid]; dom.Ey[lid] = global_Ey[gid]; dom.Ez[lid] = global_Ez[gid];
                    dom.Hx[lid] = global_Hx[gid]; dom.Hy[lid] = global_Hy[gid]; dom.Hz[lid] = global_Hz[gid];
                }
            }
        }
    }

    // Gathers only interior cells, not halo regions
    void gather_from_subdomain(const Subdomain& dom,
                              std::vector<real>& global_Ex,
                              std::vector<real>& global_Ey,
                              std::vector<real>& global_Ez,
                              std::vector<real>& global_Hx,
                              std::vector<real>& global_Hy,
                              std::vector<real>& global_Hz) {
        size_t i_start = dom.is_left_boundary ? 0 : config.halo_width;
        size_t i_end = dom.local_Nx - (dom.is_right_boundary ? 0 : config.halo_width);
        size_t j_start = 0, j_end = dom.local_Ny;
        size_t k_start = 0, k_end = dom.local_Nz;

        switch (config.axis) {
            case DecompAxis::X:
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

                    global_Ex[gid] = dom.Ex[lid]; global_Ey[gid] = dom.Ey[lid]; global_Ez[gid] = dom.Ez[lid];
                    global_Hx[gid] = dom.Hx[lid]; global_Hy[gid] = dom.Hy[lid]; global_Hz[gid] = dom.Hz[lid];
                }
            }
        }
    }

    void exchange_halos_internal(bool after_H) {
        for (int d = 0; d < config.num_domains - 1; ++d) {
            Subdomain& left = subdomains[d];
            Subdomain& right = subdomains[d + 1];

            switch (config.axis) {
                case DecompAxis::X: exchange_halos_X(left, right, after_H); break;
                case DecompAxis::Y: exchange_halos_Y(left, right, after_H); break;
                case DecompAxis::Z: exchange_halos_Z(left, right, after_H); break;
            }
        }
    }

    void exchange_halos_X(Subdomain& left, Subdomain& right, bool after_H) {
        size_t hw = config.halo_width;

        for (size_t j = 0; j < NyT; ++j) {
            for (size_t k = 0; k < NzT; ++k) {
                for (size_t h = 0; h < hw; ++h) {
                    size_t left_src = left.local_idx(left.local_Nx - 2*hw + h, j, k);
                    size_t right_dst = right.local_idx(h, j, k);
                    size_t right_src = right.local_idx(hw + h, j, k);
                    size_t left_dst = left.local_idx(left.local_Nx - hw + h, j, k);

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

template<typename Real>
inline void fdtd_update_H_subdomain(Subdomain& dom) {
    using namespace fdtd_math;

    const size_t NxL = dom.local_Nx, NyL = dom.local_Ny, NzL = dom.local_Nz;
    const size_t sI = NyL * NzL, sJ = NzL, sK = 1;

    const Real* __restrict dx_arr = dom.dx_local.data();
    const Real* __restrict dy_arr = dom.dy_local.data();
    const Real* __restrict dz_arr = dom.dz_local.data();

    Real* __restrict Ex = dom.Ex.data(); Real* __restrict Ey = dom.Ey.data(); Real* __restrict Ez = dom.Ez.data();
    Real* __restrict Hx = dom.Hx.data(); Real* __restrict Hy = dom.Hy.data(); Real* __restrict Hz = dom.Hz.data();

    const Real* __restrict aHx = dom.aHx.data(); const Real* __restrict bHx = dom.bHx.data();
    const Real* __restrict aHy = dom.aHy.data(); const Real* __restrict bHy = dom.bHy.data();
    const Real* __restrict aHz = dom.aHz.data(); const Real* __restrict bHz = dom.bHz.data();

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

template<typename Real>
inline void fdtd_update_E_subdomain(Subdomain& dom) {
    using namespace fdtd_math;

    const size_t NxL = dom.local_Nx, NyL = dom.local_Ny, NzL = dom.local_Nz;
    const size_t sI = NyL * NzL, sJ = NzL, sK = 1;

    const Real* __restrict dx_arr = dom.dx_local.data();
    const Real* __restrict dy_arr = dom.dy_local.data();
    const Real* __restrict dz_arr = dom.dz_local.data();

    Real* __restrict Ex = dom.Ex.data(); Real* __restrict Ey = dom.Ey.data(); Real* __restrict Ez = dom.Ez.data();
    const Real* __restrict Hx = dom.Hx.data(); const Real* __restrict Hy = dom.Hy.data(); const Real* __restrict Hz = dom.Hz.data();
    const Real* __restrict Jx = dom.Jx.data(); const Real* __restrict Jy = dom.Jy.data(); const Real* __restrict Jz = dom.Jz.data();

    const Real* __restrict aEx = dom.aEx.data(); const Real* __restrict bEx = dom.bEx.data();
    const Real* __restrict aEy = dom.aEy.data(); const Real* __restrict bEy = dom.bEy.data();
    const Real* __restrict aEz = dom.aEz.data(); const Real* __restrict bEz = dom.bEz.data();

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

template<typename Func>
void parallel_for_domains(DomainDecomposition& decomp, Func&& func) {
    if (!decomp.is_parallel()) return;

#if FDTD_OMP_ENABLED
    #pragma omp parallel for
#endif
    for (int d = 0; d < decomp.num_domains(); ++d)
        func(decomp.subdomains[d]);
}

inline void parallel_fdtd_step_H(DomainDecomposition& decomp) {
    if (!decomp.is_parallel()) return;
    parallel_for_domains(decomp, [](Subdomain& dom) { fdtd_update_H_subdomain<real>(dom); });
    decomp.exchange_halos_H();
}

inline void parallel_fdtd_step_E(DomainDecomposition& decomp) {
    if (!decomp.is_parallel()) return;
    parallel_for_domains(decomp, [](Subdomain& dom) { fdtd_update_E_subdomain<real>(dom); });
    decomp.exchange_halos_E();
}

inline void parallel_clear_J(DomainDecomposition& decomp) {
    if (!decomp.is_parallel()) return;
    parallel_for_domains(decomp, [](Subdomain& dom) {
        std::fill(dom.Jx.begin(), dom.Jx.end(), 0.0);
        std::fill(dom.Jy.begin(), dom.Jy.end(), 0.0);
        std::fill(dom.Jz.begin(), dom.Jz.end(), 0.0);
    });
}

inline void parallel_inject_source(DomainDecomposition& decomp,
                                   size_t global_i, size_t global_j, size_t global_k,
                                   real Jz_value) {
    if (!decomp.is_parallel()) return;

    for (auto& dom : decomp.subdomains) {
        bool found = false;

        switch (decomp.config.axis) {
            case DecompAxis::X:
                if (global_i >= dom.offset_i && global_i < dom.offset_i + dom.local_Nx) {
                    size_t li = global_i - dom.offset_i;
                    dom.Jz[dom.local_idx(li, global_j, global_k)] += Jz_value;
                    found = true;
                }
                break;
            case DecompAxis::Y:
                if (global_j >= dom.offset_j && global_j < dom.offset_j + dom.local_Ny) {
                    size_t lj = global_j - dom.offset_j;
                    dom.Jz[dom.local_idx(global_i, lj, global_k)] += Jz_value;
                    found = true;
                }
                break;
            case DecompAxis::Z:
                if (global_k >= dom.offset_k && global_k < dom.offset_k + dom.local_Nz) {
                    size_t lk = global_k - dom.offset_k;
                    dom.Jz[dom.local_idx(global_i, global_j, lk)] += Jz_value;
                    found = true;
                }
                break;
        }

        if (found) break;
    }
}

} // namespace Parallel
