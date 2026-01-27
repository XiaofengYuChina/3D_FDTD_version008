// parallel_domain_v2.hpp - Domain Decomposition Parallel Module for 3D FDTD
//
// TRUE domain decomposition: each subdomain is independent with its own field arrays.
// Only halo regions are exchanged between neighbors (O(N^2) per step).
// Compatible with: Non-uniform meshes, two-level system gain medium

#pragma once

#include <vector>
#include <cstddef>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <functional>

#include "omp_config.hpp"
#include "global_function.hpp"
#include "two_level_system.hpp"

namespace ParallelV2 {

enum class DecompAxis { X = 0, Y = 1, Z = 2 };

struct ParallelConfig {
    bool enabled = false;
    int num_domains = 1;
    DecompAxis axis = DecompAxis::Z;
    int halo_width = 1;
    int npml = 8;

    void validate(size_t NxT, size_t NyT, size_t NzT) {
        if (!enabled || num_domains <= 1) {
            enabled = false;
            num_domains = 1;
            return;
        }

        size_t axis_size = (axis == DecompAxis::X) ? NxT :
                          (axis == DecompAxis::Y) ? NyT : NzT;

        size_t interior_size = axis_size - 2 * npml;
        size_t min_cells_per_domain = 4 + 2 * halo_width;
        size_t max_domains = interior_size / min_cells_per_domain;

        if (num_domains > (int)max_domains) {
            std::cout << "[ParallelV2] Reducing domains from " << num_domains
                      << " to " << max_domains << "\n";
            num_domains = (int)max_domains;
        }

        if (num_domains <= 1) {
            enabled = false;
            num_domains = 1;
        }
    }
};

// CPML auxiliary field data for a single subdomain
struct SubdomainCPML {
    // 1D profiles
    std::vector<real> kx, ky, kz;
    std::vector<real> inv_kx, inv_ky, inv_kz;  // Pre-computed inverse kappa
    std::vector<real> bEx, cEx, bEy, cEy, bEz, cEz;
    std::vector<real> bHx, cHx, bHy, cHy, bHz, cHz;

    // 3D psi arrays
    std::vector<real> psi_Ex_y, psi_Ex_z;
    std::vector<real> psi_Ey_z, psi_Ey_x;
    std::vector<real> psi_Ez_x, psi_Ez_y;
    std::vector<real> psi_Hx_y, psi_Hx_z;
    std::vector<real> psi_Hy_z, psi_Hy_x;
    std::vector<real> psi_Hz_x, psi_Hz_y;

    // PML boundaries along decomposition axis
    bool has_left_pml = false;
    bool has_right_pml = false;

    real dt, eps0, mu0, c0;
    int npml;
    real m = 3.0;
    real kappa_max = 8.0;
    real alpha0 = 0.05;

    void initialize(size_t local_Nx, size_t local_Ny, size_t local_Nz,
                   const std::vector<real>& dx_local,
                   const std::vector<real>& dy_local,
                   const std::vector<real>& dz_local,
                   size_t global_offset_along_axis,
                   size_t global_axis_size,
                   DecompAxis axis,
                   int pml_thickness,
                   real dt_, real eps0_, real mu0_, real c0_) {
        dt = dt_; eps0 = eps0_; mu0 = mu0_; c0 = c0_;
        npml = pml_thickness;

        size_t Ntot = local_Nx * local_Ny * local_Nz;

        kx.assign(local_Nx, 1.0);
        ky.assign(local_Ny, 1.0);
        kz.assign(local_Nz, 1.0);

        bEx.assign(local_Nx, 1.0); cEx.assign(local_Nx, 0.0);
        bEy.assign(local_Ny, 1.0); cEy.assign(local_Ny, 0.0);
        bEz.assign(local_Nz, 1.0); cEz.assign(local_Nz, 0.0);
        bHx.assign(local_Nx, 1.0); cHx.assign(local_Nx, 0.0);
        bHy.assign(local_Ny, 1.0); cHy.assign(local_Ny, 0.0);
        bHz.assign(local_Nz, 1.0); cHz.assign(local_Nz, 0.0);

        psi_Ex_y.assign(Ntot, 0.0); psi_Ex_z.assign(Ntot, 0.0);
        psi_Ey_z.assign(Ntot, 0.0); psi_Ey_x.assign(Ntot, 0.0);
        psi_Ez_x.assign(Ntot, 0.0); psi_Ez_y.assign(Ntot, 0.0);
        psi_Hx_y.assign(Ntot, 0.0); psi_Hx_z.assign(Ntot, 0.0);
        psi_Hy_z.assign(Ntot, 0.0); psi_Hy_x.assign(Ntot, 0.0);
        psi_Hz_x.assign(Ntot, 0.0); psi_Hz_y.assign(Ntot, 0.0);

        has_left_pml = (global_offset_along_axis == 0);
        has_right_pml = (global_offset_along_axis +
            (axis == DecompAxis::X ? local_Nx : axis == DecompAxis::Y ? local_Ny : local_Nz)
            >= global_axis_size);

        build_profile_1d(local_Nx, dx_local,
                        axis == DecompAxis::X ? has_left_pml : true,
                        axis == DecompAxis::X ? has_right_pml : true,
                        axis == DecompAxis::X ? global_offset_along_axis : 0,
                        axis == DecompAxis::X ? global_axis_size : local_Nx,
                        kx, bEx, cEx, bHx, cHx, true);

        build_profile_1d(local_Ny, dy_local,
                        axis == DecompAxis::Y ? has_left_pml : true,
                        axis == DecompAxis::Y ? has_right_pml : true,
                        axis == DecompAxis::Y ? global_offset_along_axis : 0,
                        axis == DecompAxis::Y ? global_axis_size : local_Ny,
                        ky, bEy, cEy, bHy, cHy, true);

        build_profile_1d(local_Nz, dz_local,
                        axis == DecompAxis::Z ? has_left_pml : true,
                        axis == DecompAxis::Z ? has_right_pml : true,
                        axis == DecompAxis::Z ? global_offset_along_axis : 0,
                        axis == DecompAxis::Z ? global_axis_size : local_Nz,
                        kz, bEz, cEz, bHz, cHz, true);

        inv_kx.resize(local_Nx);
        inv_ky.resize(local_Ny);
        inv_kz.resize(local_Nz);
        for (size_t i = 0; i < local_Nx; ++i) inv_kx[i] = 1.0 / kx[i];
        for (size_t j = 0; j < local_Ny; ++j) inv_ky[j] = 1.0 / ky[j];
        for (size_t k = 0; k < local_Nz; ++k) inv_kz[k] = 1.0 / kz[k];
    }

private:
    void build_profile_1d(size_t local_N, const std::vector<real>& ds_local,
                         bool has_low_pml, bool has_high_pml,
                         size_t global_offset, size_t global_N,
                         std::vector<real>& k,
                         std::vector<real>& bE, std::vector<real>& cE,
                         std::vector<real>& bH, std::vector<real>& cH,
                         bool /*is_decomp_axis*/) {

        // Calculate average spacing in PML region for sigma_max
        real ds_avg = 0;
        int count = 0;
        for (size_t i = 0; i < local_N && i < (size_t)npml; ++i) {
            ds_avg += ds_local[i];
            count++;
        }
        if (count > 0) ds_avg /= count;
        else ds_avg = ds_local[0];

        real sigma_max = ((m + 1.0) * std::log(1e10) * eps0 * c0) / (2.0 * npml * ds_avg);

        for (size_t n = 0; n < local_N; ++n) {
            size_t global_n = global_offset + n;
            int dist_from_low = (int)global_n;
            int dist_from_high = (int)(global_N - 1 - global_n);

            bool in_low_pml = has_low_pml && (dist_from_low < npml);
            bool in_high_pml = has_high_pml && (dist_from_high < npml);

            if (!in_low_pml && !in_high_pml) {
                k[n] = 1.0;
                bE[n] = 1.0; cE[n] = 0.0;
                bH[n] = 1.0; cH[n] = 0.0;
                continue;
            }

            // Use the closer PML boundary
            int d = in_low_pml ? dist_from_low : dist_from_high;

            real r = real(npml - (d + 0.5)) / real(npml);
            real rm = std::pow(r, m);
            real kap = 1.0 + (kappa_max - 1.0) * rm;
            real sig = sigma_max * kap * rm;
            real alf = alpha0 * (c0 / ds_local[n]) * (1.0 - r);

            k[n] = kap;

            real sig_e = sig / eps0;
            real be = std::exp(-(sig_e / kap + alf) * dt);
            be = std::max<real>(be, 1e-30);
            real ce = 0.0;
            if (std::abs(sig_e) > 1e-30 || std::abs(alf) > 1e-30) {
                ce = (sig_e * (be - 1.0)) / (kap * (sig_e + kap * alf));
            }
            bE[n] = be;
            cE[n] = ce;

            real sig_h = sig / mu0;
            real bh = std::exp(-(sig_h / kap + alf) * dt);
            bh = std::max<real>(bh, 1e-30);
            real ch = 0.0;
            if (std::abs(sig_h) > 1e-30 || std::abs(alf) > 1e-30) {
                ch = (sig_h * (bh - 1.0)) / (kap * (sig_h + kap * alf));
            }
            bH[n] = bh;
            cH[n] = ch;
        }
    }
};

struct Subdomain {
    int id;

    // Global position info
    size_t global_start, global_end;
    size_t offset_i, offset_j, offset_k;

    // Local dimensions
    size_t local_Nx, local_Ny, local_Nz;
    size_t local_Ntot;

    // Neighbor info (-1 if at boundary)
    int left_neighbor, right_neighbor;
    size_t halo_width;

    bool is_left_boundary;
    bool is_right_boundary;

    // Local field arrays
    std::vector<real> Ex, Ey, Ez;
    std::vector<real> Hx, Hy, Hz;
    std::vector<real> Jx, Jy, Jz;

    // Local coefficient arrays
    std::vector<real> aEx, bEx, aEy, bEy, aEz, bEz;
    std::vector<real> aHx, bHx, aHy, bHy, aHz, bHz;

    // Local grid spacing and pre-computed inverses
    std::vector<real> dx_local, dy_local, dz_local;
    std::vector<real> inv_dx_local, inv_dy_local, inv_dz_local;

    SubdomainCPML cpml;

    // TLS (Two-Level System) local data - NO halo exchange needed
    bool tls_enabled = false;
    std::vector<real> tls_Ng, tls_Nu, tls_Ng0, tls_Ndip;
    std::vector<real> tls_Pz, tls_Pz_prev, tls_dPz_dt, tls_Ez_old;

    void allocate_tls() {
        tls_Ng.assign(local_Ntot, 0);
        tls_Nu.assign(local_Ntot, 0);
        tls_Ng0.assign(local_Ntot, 0);
        tls_Ndip.assign(local_Ntot, 0);
        tls_Pz.assign(local_Ntot, 0);
        tls_Pz_prev.assign(local_Ntot, 0);
        tls_dPz_dt.assign(local_Ntot, 0);
        tls_Ez_old.assign(local_Ntot, 0);
    }

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

    bool contains_global(size_t gi, size_t gj, size_t gk, DecompAxis axis) const {
        switch (axis) {
            case DecompAxis::X: return gi >= offset_i && gi < offset_i + local_Nx;
            case DecompAxis::Y: return gj >= offset_j && gj < offset_j + local_Ny;
            case DecompAxis::Z: return gk >= offset_k && gk < offset_k + local_Nz;
        }
        return false;
    }
};

class DomainDecomposition {
public:
    ParallelConfig config;
    size_t NxT, NyT, NzT, Ntot;
    size_t npml;
    std::vector<Subdomain> subdomains;
    real dt, eps0, mu0, c0;

    DomainDecomposition() = default;

    void initialize(const ParallelConfig& cfg,
                   size_t Nx, size_t Ny, size_t Nz,
                   size_t pml_thickness,
                   const GridSpacing& grid_spacing,
                   const MaterialGrids& mats,
                   real dt_, real eps0_, real mu0_, real c0_) {

        config = cfg;
        NxT = Nx; NyT = Ny; NzT = Nz;
        Ntot = NxT * NyT * NzT;
        npml = pml_thickness;
        dt = dt_; eps0 = eps0_; mu0 = mu0_; c0 = c0_;

        config.npml = pml_thickness;
        config.validate(NxT, NyT, NzT);

        if (!config.enabled) {
            std::cout << "[ParallelV2] Domain decomposition disabled\n";
            return;
        }

        std::cout << "\n[ParallelV2] Initializing TRUE domain decomposition:\n";
        std::cout << "  Domains: " << config.num_domains << "\n";
        std::cout << "  Axis: " << (config.axis == DecompAxis::X ? "X" :
                                    config.axis == DecompAxis::Y ? "Y" : "Z") << "\n";
        std::cout << "  Halo width: " << config.halo_width << "\n";
        std::cout << "  Each subdomain has independent PML (no scatter/gather needed)\n";

        create_subdomains(grid_spacing, mats);

        std::cout << "[ParallelV2] Initialization complete\n\n";
    }

    bool is_parallel() const { return config.enabled && config.num_domains > 1; }
    int num_domains() const { return config.num_domains; }

    void initialize_tls(const TwoLevelState& global_tls) {
        if (!is_parallel()) return;

        std::cout << "[ParallelV2] Initializing TLS for " << config.num_domains << " subdomains...\n";

        for (auto& dom : subdomains) {
            dom.tls_enabled = true;
            dom.allocate_tls();

            for (size_t li = 0; li < dom.local_Nx; ++li) {
                for (size_t lj = 0; lj < dom.local_Ny; ++lj) {
                    for (size_t lk = 0; lk < dom.local_Nz; ++lk) {
                        size_t gi, gj, gk;
                        dom.local_to_global(li, lj, lk, gi, gj, gk);

                        size_t lid = dom.local_idx(li, lj, lk);
                        size_t gid = idx3(gi, gj, gk, NyT, NzT);

                        dom.tls_Ng[lid] = global_tls.Ng[gid];
                        dom.tls_Nu[lid] = global_tls.Nu[gid];
                        dom.tls_Ng0[lid] = global_tls.Ng0[gid];
                        dom.tls_Ndip[lid] = global_tls.Ndip[gid];
                        dom.tls_Pz[lid] = global_tls.Pz[gid];
                        dom.tls_Pz_prev[lid] = global_tls.Pz_prev[gid];
                        dom.tls_dPz_dt[lid] = global_tls.dPz_dt[gid];
                        dom.tls_Ez_old[lid] = global_tls.Ez_old[gid];
                    }
                }
            }
        }

        std::cout << "[ParallelV2] TLS initialization complete\n";
    }

    // Gather TLS state back to global arrays (for diagnostics/output)
    void gather_tls_for_output(TwoLevelState& global_tls) {
        if (!is_parallel()) return;

        for (const auto& dom : subdomains) {
            if (!dom.tls_enabled) continue;

            size_t i_start, i_end, j_start, j_end, k_start, k_end;
            get_interior_range(dom, i_start, i_end, j_start, j_end, k_start, k_end);

            for (size_t li = i_start; li < i_end; ++li) {
                for (size_t lj = j_start; lj < j_end; ++lj) {
                    for (size_t lk = k_start; lk < k_end; ++lk) {
                        size_t gi, gj, gk;
                        dom.local_to_global(li, lj, lk, gi, gj, gk);

                        size_t lid = dom.local_idx(li, lj, lk);
                        size_t gid = idx3(gi, gj, gk, NyT, NzT);

                        global_tls.Ng[gid] = dom.tls_Ng[lid];
                        global_tls.Nu[gid] = dom.tls_Nu[lid];
                        global_tls.Pz[gid] = dom.tls_Pz[lid];
                        global_tls.Pz_prev[gid] = dom.tls_Pz_prev[lid];
                        global_tls.dPz_dt[gid] = dom.tls_dPz_dt[lid];
                    }
                }
            }
        }
    }

    void exchange_halos() {
        if (!is_parallel()) return;

        for (int d = 0; d < config.num_domains - 1; ++d) {
            exchange_halo_pair(subdomains[d], subdomains[d + 1]);
        }
    }

    // Gather fields to global arrays (only for output/diagnostics)
    void gather_fields_for_output(std::vector<real>& global_Ex,
                                  std::vector<real>& global_Ey,
                                  std::vector<real>& global_Ez,
                                  std::vector<real>& global_Hx,
                                  std::vector<real>& global_Hy,
                                  std::vector<real>& global_Hz) {
        if (!is_parallel()) return;

        for (const auto& dom : subdomains) {
            size_t i_start, i_end, j_start, j_end, k_start, k_end;
            get_interior_range(dom, i_start, i_end, j_start, j_end, k_start, k_end);

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
    }

    // Initial scatter - only called once at simulation start
    void initial_scatter(const std::vector<real>& global_Ex,
                        const std::vector<real>& global_Ey,
                        const std::vector<real>& global_Ez,
                        const std::vector<real>& global_Hx,
                        const std::vector<real>& global_Hy,
                        const std::vector<real>& global_Hz) {
        if (!is_parallel()) return;

        for (auto& dom : subdomains) {
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
    }

    // Get interior range excluding halo cells
    void get_interior_range(const Subdomain& dom,
                           size_t& i_start, size_t& i_end,
                           size_t& j_start, size_t& j_end,
                           size_t& k_start, size_t& k_end) const {
        i_start = 0; i_end = dom.local_Nx;
        j_start = 0; j_end = dom.local_Ny;
        k_start = 0; k_end = dom.local_Nz;

        switch (config.axis) {
            case DecompAxis::X:
                if (!dom.is_left_boundary) i_start = config.halo_width;
                if (!dom.is_right_boundary) i_end = dom.local_Nx - config.halo_width;
                break;
            case DecompAxis::Y:
                if (!dom.is_left_boundary) j_start = config.halo_width;
                if (!dom.is_right_boundary) j_end = dom.local_Ny - config.halo_width;
                break;
            case DecompAxis::Z:
                if (!dom.is_left_boundary) k_start = config.halo_width;
                if (!dom.is_right_boundary) k_end = dom.local_Nz - config.halo_width;
                break;
        }
    }

private:
    void create_subdomains(const GridSpacing& grid_spacing, const MaterialGrids& mats) {
        subdomains.clear();
        subdomains.resize(config.num_domains);

        size_t axis_size = (config.axis == DecompAxis::X) ? NxT :
                          (config.axis == DecompAxis::Y) ? NyT : NzT;

        size_t interior_start = npml;
        size_t interior_end = axis_size - npml;
        size_t interior_size = interior_end - interior_start;

        size_t base_size = interior_size / config.num_domains;
        size_t remainder = interior_size % config.num_domains;
        size_t current_pos = 0;

        for (int d = 0; d < config.num_domains; ++d) {
            Subdomain& dom = subdomains[d];
            dom.id = d;
            dom.halo_width = config.halo_width;

            size_t dom_interior = base_size + (d < (int)remainder ? 1 : 0);

            if (d == 0) {
                dom.global_start = 0;
                dom.is_left_boundary = true;
            } else {
                dom.global_start = interior_start + current_pos - config.halo_width;
                dom.is_left_boundary = false;
            }

            current_pos += dom_interior;

            if (d == config.num_domains - 1) {
                dom.global_end = axis_size;
                dom.is_right_boundary = true;
            } else {
                dom.global_end = interior_start + current_pos + config.halo_width;
                dom.is_right_boundary = false;
            }

            dom.left_neighbor = (d > 0) ? d - 1 : -1;
            dom.right_neighbor = (d < config.num_domains - 1) ? d + 1 : -1;

            set_local_dimensions(dom, grid_spacing);
            dom.allocate();
            copy_grid_spacing(dom, grid_spacing);
            copy_material_coeffs(dom, mats);

            dom.cpml.initialize(
                dom.local_Nx, dom.local_Ny, dom.local_Nz,
                dom.dx_local, dom.dy_local, dom.dz_local,
                (config.axis == DecompAxis::X) ? dom.global_start :
                (config.axis == DecompAxis::Y) ? dom.global_start : dom.global_start,
                axis_size, config.axis, npml, dt, eps0, mu0, c0
            );

            std::cout << "  Domain " << d << ": [" << dom.global_start << ", "
                      << dom.global_end << "), size " << dom.local_Nx << "x"
                      << dom.local_Ny << "x" << dom.local_Nz;
            if (dom.is_left_boundary) std::cout << " [LEFT PML]";
            if (dom.is_right_boundary) std::cout << " [RIGHT PML]";
            std::cout << "\n";
        }
    }

    void set_local_dimensions(Subdomain& dom, const GridSpacing& /*gs*/) {
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
    }

    void copy_grid_spacing(Subdomain& dom, const GridSpacing& gs) {
        switch (config.axis) {
            case DecompAxis::X:
                dom.dx_local.assign(gs.dx.begin() + dom.global_start, gs.dx.begin() + dom.global_end);
                dom.dy_local = gs.dy;
                dom.dz_local = gs.dz;
                dom.inv_dx_local.assign(gs.inv_dx.begin() + dom.global_start, gs.inv_dx.begin() + dom.global_end);
                dom.inv_dy_local = gs.inv_dy;
                dom.inv_dz_local = gs.inv_dz;
                break;
            case DecompAxis::Y:
                dom.dx_local = gs.dx;
                dom.dy_local.assign(gs.dy.begin() + dom.global_start, gs.dy.begin() + dom.global_end);
                dom.dz_local = gs.dz;
                dom.inv_dx_local = gs.inv_dx;
                dom.inv_dy_local.assign(gs.inv_dy.begin() + dom.global_start, gs.inv_dy.begin() + dom.global_end);
                dom.inv_dz_local = gs.inv_dz;
                break;
            case DecompAxis::Z:
                dom.dx_local = gs.dx;
                dom.dy_local = gs.dy;
                dom.dz_local.assign(gs.dz.begin() + dom.global_start, gs.dz.begin() + dom.global_end);
                dom.inv_dx_local = gs.inv_dx;
                dom.inv_dy_local = gs.inv_dy;
                dom.inv_dz_local.assign(gs.inv_dz.begin() + dom.global_start, gs.inv_dz.begin() + dom.global_end);
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

    void exchange_halo_pair(Subdomain& left, Subdomain& right) {
        size_t hw = config.halo_width;

        switch (config.axis) {
            case DecompAxis::X: exchange_halo_X(left, right, hw); break;
            case DecompAxis::Y: exchange_halo_Y(left, right, hw); break;
            case DecompAxis::Z: exchange_halo_Z(left, right, hw); break;
        }
    }

    void exchange_halo_X(Subdomain& left, Subdomain& right, size_t hw) {
        for (size_t j = 0; j < NyT; ++j) {
            for (size_t k = 0; k < NzT; ++k) {
                for (size_t h = 0; h < hw; ++h) {
                    size_t left_src = left.local_idx(left.local_Nx - 2*hw + h, j, k);
                    size_t right_dst = right.local_idx(h, j, k);
                    size_t right_src = right.local_idx(hw + h, j, k);
                    size_t left_dst = left.local_idx(left.local_Nx - hw + h, j, k);

                    std::swap(right.Ex[right_dst], left.Ex[left_src]); right.Ex[right_dst] = left.Ex[left_src];
                    std::swap(right.Ey[right_dst], left.Ey[left_src]); right.Ey[right_dst] = left.Ey[left_src];
                    std::swap(right.Ez[right_dst], left.Ez[left_src]); right.Ez[right_dst] = left.Ez[left_src];
                    std::swap(right.Hx[right_dst], left.Hx[left_src]); right.Hx[right_dst] = left.Hx[left_src];
                    std::swap(right.Hy[right_dst], left.Hy[left_src]); right.Hy[right_dst] = left.Hy[left_src];
                    std::swap(right.Hz[right_dst], left.Hz[left_src]); right.Hz[right_dst] = left.Hz[left_src];

                    left.Ex[left_dst] = right.Ex[right_src];
                    left.Ey[left_dst] = right.Ey[right_src];
                    left.Ez[left_dst] = right.Ez[right_src];
                    left.Hx[left_dst] = right.Hx[right_src];
                    left.Hy[left_dst] = right.Hy[right_src];
                    left.Hz[left_dst] = right.Hz[right_src];
                }
            }
        }
    }

    void exchange_halo_Y(Subdomain& left, Subdomain& right, size_t hw) {
        for (size_t i = 0; i < NxT; ++i) {
            for (size_t k = 0; k < NzT; ++k) {
                for (size_t h = 0; h < hw; ++h) {
                    size_t left_src = left.local_idx(i, left.local_Ny - 2*hw + h, k);
                    size_t right_dst = right.local_idx(i, h, k);
                    size_t right_src = right.local_idx(i, hw + h, k);
                    size_t left_dst = left.local_idx(i, left.local_Ny - hw + h, k);

                    right.Ex[right_dst] = left.Ex[left_src];
                    right.Ey[right_dst] = left.Ey[left_src];
                    right.Ez[right_dst] = left.Ez[left_src];
                    right.Hx[right_dst] = left.Hx[left_src];
                    right.Hy[right_dst] = left.Hy[left_src];
                    right.Hz[right_dst] = left.Hz[left_src];

                    left.Ex[left_dst] = right.Ex[right_src];
                    left.Ey[left_dst] = right.Ey[right_src];
                    left.Ez[left_dst] = right.Ez[right_src];
                    left.Hx[left_dst] = right.Hx[right_src];
                    left.Hy[left_dst] = right.Hy[right_src];
                    left.Hz[left_dst] = right.Hz[right_src];
                }
            }
        }
    }

    void exchange_halo_Z(Subdomain& left, Subdomain& right, size_t hw) {
        for (size_t i = 0; i < NxT; ++i) {
            for (size_t j = 0; j < NyT; ++j) {
                for (size_t h = 0; h < hw; ++h) {
                    size_t left_src = left.local_idx(i, j, left.local_Nz - 2*hw + h);
                    size_t right_dst = right.local_idx(i, j, h);
                    size_t right_src = right.local_idx(i, j, hw + h);
                    size_t left_dst = left.local_idx(i, j, left.local_Nz - hw + h);

                    right.Ex[right_dst] = left.Ex[left_src];
                    right.Ey[right_dst] = left.Ey[left_src];
                    right.Ez[right_dst] = left.Ez[left_src];
                    right.Hx[right_dst] = left.Hx[left_src];
                    right.Hy[right_dst] = left.Hy[left_src];
                    right.Hz[right_dst] = left.Hz[left_src];

                    left.Ex[left_dst] = right.Ex[right_src];
                    left.Ey[left_dst] = right.Ey[right_src];
                    left.Ez[left_dst] = right.Ez[right_src];
                    left.Hx[left_dst] = right.Hx[right_src];
                    left.Hy[left_dst] = right.Hy[right_src];
                    left.Hz[left_dst] = right.Hz[right_src];
                }
            }
        }
    }
};

// H-field update for subdomain with CPML
inline void fdtd_update_H_subdomain(Subdomain& dom) {
    using namespace fdtd_math;

    const size_t NxL = dom.local_Nx;
    const size_t NyL = dom.local_Ny;
    const size_t NzL = dom.local_Nz;
    const size_t sI = NyL * NzL;
    const size_t sJ = NzL;
    const size_t sK = 1;

    const real* inv_dx = dom.inv_dx_local.data();
    const real* inv_dy = dom.inv_dy_local.data();
    const real* inv_dz = dom.inv_dz_local.data();

    real* Ex = dom.Ex.data(); real* Ey = dom.Ey.data(); real* Ez = dom.Ez.data();
    real* Hx = dom.Hx.data(); real* Hy = dom.Hy.data(); real* Hz = dom.Hz.data();

    const real* aHx = dom.aHx.data(); const real* bHx = dom.bHx.data();
    const real* aHy = dom.aHy.data(); const real* bHy = dom.bHy.data();
    const real* aHz = dom.aHz.data(); const real* bHz = dom.bHz.data();

    SubdomainCPML& cpml = dom.cpml;
    const real* inv_kx = cpml.inv_kx.data();
    const real* inv_ky = cpml.inv_ky.data();
    const real* inv_kz = cpml.inv_kz.data();
    const real* bHx_pml = cpml.bHx.data(); const real* cHx_pml = cpml.cHx.data();
    const real* bHy_pml = cpml.bHy.data(); const real* cHy_pml = cpml.cHy.data();
    const real* bHz_pml = cpml.bHz.data(); const real* cHz_pml = cpml.cHz.data();

    real* psi_Hx_y = cpml.psi_Hx_y.data(); real* psi_Hx_z = cpml.psi_Hx_z.data();
    real* psi_Hy_z = cpml.psi_Hy_z.data(); real* psi_Hy_x = cpml.psi_Hy_x.data();
    real* psi_Hz_x = cpml.psi_Hz_x.data(); real* psi_Hz_y = cpml.psi_Hz_y.data();

#if FDTD_OMP_ENABLED
#pragma omp parallel for
#endif
    for (int i = 0; i < (int)(NxL - 1); ++i) {
        for (int j = 0; j < (int)(NyL - 1); ++j) {
            for (int k = 0; k < (int)(NzL - 1); ++k) {
                size_t id = ((size_t)i * NyL + (size_t)j) * NzL + (size_t)k;

                const real inv_dx_i = inv_dx[i];
                const real inv_dy_j = inv_dy[j];
                const real inv_dz_k = inv_dz[k];

                real dEz_dy = (Ez[id + sJ] - Ez[id]) * inv_dy_j;
                real dEy_dz = (Ey[id + sK] - Ey[id]) * inv_dz_k;
                real dEx_dz = (Ex[id + sK] - Ex[id]) * inv_dz_k;
                real dEz_dx = (Ez[id + sI] - Ez[id]) * inv_dx_i;
                real dEy_dx = (Ey[id + sI] - Ey[id]) * inv_dx_i;
                real dEx_dy = (Ex[id + sJ] - Ex[id]) * inv_dy_j;

                // Hx: dEz/dy - dEy/dz
                real corr_y = dEz_dy * inv_ky[j];
                psi_Hx_y[id] = bHy_pml[j] * psi_Hx_y[id] + cHy_pml[j] * dEz_dy;
                corr_y += psi_Hx_y[id];

                real corr_z = dEy_dz * inv_kz[k];
                psi_Hx_z[id] = bHz_pml[k] * psi_Hx_z[id] + cHz_pml[k] * dEy_dz;
                corr_z += psi_Hx_z[id];

                real curlEx = corr_y - corr_z;

                // Hy: dEx/dz - dEz/dx
                corr_z = dEx_dz * inv_kz[k];
                psi_Hy_z[id] = bHz_pml[k] * psi_Hy_z[id] + cHz_pml[k] * dEx_dz;
                corr_z += psi_Hy_z[id];

                real corr_x = dEz_dx * inv_kx[i];
                psi_Hy_x[id] = bHx_pml[i] * psi_Hy_x[id] + cHx_pml[i] * dEz_dx;
                corr_x += psi_Hy_x[id];

                real curlEy = corr_z - corr_x;

                // Hz: dEy/dx - dEx/dy
                corr_x = dEy_dx * inv_kx[i];
                psi_Hz_x[id] = bHx_pml[i] * psi_Hz_x[id] + cHx_pml[i] * dEy_dx;
                corr_x += psi_Hz_x[id];

                corr_y = dEx_dy * inv_ky[j];
                psi_Hz_y[id] = bHy_pml[j] * psi_Hz_y[id] + cHy_pml[j] * dEx_dy;
                corr_y += psi_Hz_y[id];

                real curlEz = corr_x - corr_y;

                Hx[id] = aHx[id] * Hx[id] - bHx[id] * curlEx;
                Hy[id] = aHy[id] * Hy[id] - bHy[id] * curlEy;
                Hz[id] = aHz[id] * Hz[id] - bHz[id] * curlEz;
            }
        }
    }
}

// E-field update for subdomain with CPML
inline void fdtd_update_E_subdomain(Subdomain& dom) {
    using namespace fdtd_math;

    const size_t NxL = dom.local_Nx;
    const size_t NyL = dom.local_Ny;
    const size_t NzL = dom.local_Nz;
    const size_t sI = NyL * NzL;
    const size_t sJ = NzL;
    const size_t sK = 1;

    const real* inv_dx = dom.inv_dx_local.data();
    const real* inv_dy = dom.inv_dy_local.data();
    const real* inv_dz = dom.inv_dz_local.data();

    real* Ex = dom.Ex.data(); real* Ey = dom.Ey.data(); real* Ez = dom.Ez.data();
    const real* Hx = dom.Hx.data(); const real* Hy = dom.Hy.data(); const real* Hz = dom.Hz.data();
    const real* Jx = dom.Jx.data(); const real* Jy = dom.Jy.data(); const real* Jz = dom.Jz.data();

    const real* aEx = dom.aEx.data(); const real* bEx = dom.bEx.data();
    const real* aEy = dom.aEy.data(); const real* bEy = dom.bEy.data();
    const real* aEz = dom.aEz.data(); const real* bEz = dom.bEz.data();

    SubdomainCPML& cpml = dom.cpml;
    const real* inv_kx = cpml.inv_kx.data();
    const real* inv_ky = cpml.inv_ky.data();
    const real* inv_kz = cpml.inv_kz.data();
    const real* bEx_pml = cpml.bEx.data(); const real* cEx_pml = cpml.cEx.data();
    const real* bEy_pml = cpml.bEy.data(); const real* cEy_pml = cpml.cEy.data();
    const real* bEz_pml = cpml.bEz.data(); const real* cEz_pml = cpml.cEz.data();

    real* psi_Ex_y = cpml.psi_Ex_y.data(); real* psi_Ex_z = cpml.psi_Ex_z.data();
    real* psi_Ey_z = cpml.psi_Ey_z.data(); real* psi_Ey_x = cpml.psi_Ey_x.data();
    real* psi_Ez_x = cpml.psi_Ez_x.data(); real* psi_Ez_y = cpml.psi_Ez_y.data();

#if FDTD_OMP_ENABLED
#pragma omp parallel for
#endif
    for (int i = 1; i < (int)NxL; ++i) {
        for (int j = 1; j < (int)NyL; ++j) {
            for (int k = 1; k < (int)NzL; ++k) {
                size_t id = ((size_t)i * NyL + (size_t)j) * NzL + (size_t)k;

                const real inv_dx_i = inv_dx[i];
                const real inv_dy_j = inv_dy[j];
                const real inv_dz_k = inv_dz[k];

                real dHz_dy = (Hz[id] - Hz[id - sJ]) * inv_dy_j;
                real dHy_dz = (Hy[id] - Hy[id - sK]) * inv_dz_k;
                real dHx_dz = (Hx[id] - Hx[id - sK]) * inv_dz_k;
                real dHz_dx = (Hz[id] - Hz[id - sI]) * inv_dx_i;
                real dHy_dx = (Hy[id] - Hy[id - sI]) * inv_dx_i;
                real dHx_dy = (Hx[id] - Hx[id - sJ]) * inv_dy_j;

                // Ex: dHz/dy - dHy/dz
                real corr_y = dHz_dy * inv_ky[j];
                psi_Ex_y[id] = bEy_pml[j] * psi_Ex_y[id] + cEy_pml[j] * dHz_dy;
                corr_y += psi_Ex_y[id];

                real corr_z = dHy_dz * inv_kz[k];
                psi_Ex_z[id] = bEz_pml[k] * psi_Ex_z[id] + cEz_pml[k] * dHy_dz;
                corr_z += psi_Ex_z[id];

                real curlHx = corr_y - corr_z;

                // Ey: dHx/dz - dHz/dx
                corr_z = dHx_dz * inv_kz[k];
                psi_Ey_z[id] = bEz_pml[k] * psi_Ey_z[id] + cEz_pml[k] * dHx_dz;
                corr_z += psi_Ey_z[id];

                real corr_x = dHz_dx * inv_kx[i];
                psi_Ey_x[id] = bEx_pml[i] * psi_Ey_x[id] + cEx_pml[i] * dHz_dx;
                corr_x += psi_Ey_x[id];

                real curlHy = corr_z - corr_x;

                // Ez: dHy/dx - dHx/dy
                corr_x = dHy_dx * inv_kx[i];
                psi_Ez_x[id] = bEx_pml[i] * psi_Ez_x[id] + cEx_pml[i] * dHy_dx;
                corr_x += psi_Ez_x[id];

                corr_y = dHx_dy * inv_ky[j];
                psi_Ez_y[id] = bEy_pml[j] * psi_Ez_y[id] + cEy_pml[j] * dHx_dy;
                corr_y += psi_Ez_y[id];

                real curlHz = corr_x - corr_y;

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
    omp_set_nested(0);
    #pragma omp parallel for
#endif
    for (int d = 0; d < decomp.num_domains(); ++d) {
        func(decomp.subdomains[d]);
    }
}

inline void parallel_update_H(DomainDecomposition& decomp) {
    if (!decomp.is_parallel()) return;

    parallel_for_domains(decomp, [](Subdomain& dom) {
        fdtd_update_H_subdomain(dom);
    });

    decomp.exchange_halos();
}

inline void parallel_update_E(DomainDecomposition& decomp) {
    if (!decomp.is_parallel()) return;

    parallel_for_domains(decomp, [](Subdomain& dom) {
        fdtd_update_E_subdomain(dom);
    });

    decomp.exchange_halos();
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
        if (dom.contains_global(global_i, global_j, global_k, decomp.config.axis)) {
            size_t li = global_i - dom.offset_i;
            size_t lj = global_j - dom.offset_j;
            size_t lk = global_k - dom.offset_k;
            size_t lid = dom.local_idx(li, lj, lk);
            dom.Jz[lid] += Jz_value;
            break;
        }
    }
}

// TLS polarization update using three-point recursion
inline void update_polarization_subdomain(Subdomain& dom, const TwoLevelParams& params) {
    if (!dom.tls_enabled || !params.coefficients_initialized) return;

    const size_t NxL = dom.local_Nx;
    const size_t NyL = dom.local_Ny;
    const size_t NzL = dom.local_Nz;

    const real kapa_coeff = params.kapa_coeff;
    const real pa2 = params.pa2;
    const real pa3 = params.pa3;
    const real inv_dt = params.inv_dt;

    real* Ez = dom.Ez.data();
    real* Ndip = dom.tls_Ndip.data();
    real* Ng = dom.tls_Ng.data();
    real* Nu = dom.tls_Nu.data();
    real* Ng0 = dom.tls_Ng0.data();
    real* Pz = dom.tls_Pz.data();
    real* Pz_prev = dom.tls_Pz_prev.data();
    real* dPz_dt = dom.tls_dPz_dt.data();

#if FDTD_OMP_ENABLED
#pragma omp parallel for
#endif
    for (int i = 0; i < (int)NxL; ++i) {
        for (int j = 0; j < (int)NyL; ++j) {
            for (int k = 0; k < (int)NzL; ++k) {
                size_t id = ((size_t)i * NyL + (size_t)j) * NzL + (size_t)k;

                if (Ndip[id] <= 0) continue;

                real Ng_frac = (Ng0[id] > 0) ? (Ng[id] - Nu[id]) / Ng0[id] : 0.0;
                real driving_factor = Ndip[id] * Ng_frac * Ez[id];

                // Three-point recursion: P^{n+1} = kapa*F + pa2*P^n + pa3*P^{n-1}
                real Pz_new = kapa_coeff * driving_factor + pa2 * Pz[id] + pa3 * Pz_prev[id];

                real Pz_old = Pz[id];
                real dPz_dt_half = (Pz_new - Pz_old) * inv_dt;

                Pz_prev[id] = Pz_old;
                Pz[id] = Pz_new;
                dPz_dt[id] = dPz_dt_half;
            }
        }
    }
}

inline void update_populations_subdomain(Subdomain& dom, real dt, const TwoLevelParams& params) {
    if (!dom.tls_enabled) return;

    const size_t NxL = dom.local_Nx;
    const size_t NyL = dom.local_Ny;
    const size_t NzL = dom.local_Nz;

    const real tau = params.tau;
    const real hbar = params.hbar;
    const real omega_a = params.omega_a;
    const real inv_dt = params.inv_dt;
    const real inv_tau = 1.0 / tau;
    const real inv_hbar_omega = 1.0 / (hbar * omega_a);

    real* Ez = dom.Ez.data();
    real* Ez_old = dom.tls_Ez_old.data();
    real* Ndip = dom.tls_Ndip.data();
    real* Ng = dom.tls_Ng.data();
    real* Nu = dom.tls_Nu.data();
    real* Ng0 = dom.tls_Ng0.data();
    real* Pz = dom.tls_Pz.data();
    real* Pz_prev = dom.tls_Pz_prev.data();

#if FDTD_OMP_ENABLED
#pragma omp parallel for
#endif
    for (int i = 0; i < (int)NxL; ++i) {
        for (int j = 0; j < (int)NyL; ++j) {
            for (int k = 0; k < (int)NzL; ++k) {
                size_t id = ((size_t)i * NyL + (size_t)j) * NzL + (size_t)k;

                if (Ndip[id] <= 0) continue;

                real Nu_curr = Nu[id];
                real Ntotal = Ng0[id];
                if (Ntotal <= 0) continue;

                // Stimulated term using energy-conserving rate form
                real E_avg = 0.5 * (Ez[id] + Ez_old[id]);
                real delta_P = Pz[id] - Pz_prev[id];
                real dP_dt = delta_P * inv_dt;
                real stim_rate = E_avg * dP_dt * inv_hbar_omega;

                // Rate equations: stim_rate > 0 = absorption, < 0 = stimulated emission
                real dNu_dt = -inv_tau * Nu_curr + stim_rate;
                real dNg_dt = +inv_tau * Nu_curr - stim_rate;

                Nu[id] += dt * dNu_dt;
                Ng[id] += dt * dNg_dt;

                Nu[id] = std::max(0.0, Nu[id]);
                Ng[id] = std::max(0.0, Ng[id]);
            }
        }
    }
}

// E-field update with TLS contribution
inline void fdtd_update_E_with_tls_subdomain(Subdomain& dom) {
    using namespace fdtd_math;

    const size_t NxL = dom.local_Nx;
    const size_t NyL = dom.local_Ny;
    const size_t NzL = dom.local_Nz;
    const size_t sI = NyL * NzL;
    const size_t sJ = NzL;
    const size_t sK = 1;

    const real* inv_dx = dom.inv_dx_local.data();
    const real* inv_dy = dom.inv_dy_local.data();
    const real* inv_dz = dom.inv_dz_local.data();

    real* Ex = dom.Ex.data(); real* Ey = dom.Ey.data(); real* Ez = dom.Ez.data();
    const real* Hx = dom.Hx.data(); const real* Hy = dom.Hy.data(); const real* Hz = dom.Hz.data();
    const real* Jx = dom.Jx.data(); const real* Jy = dom.Jy.data(); const real* Jz = dom.Jz.data();

    const real* aEx = dom.aEx.data(); const real* bEx = dom.bEx.data();
    const real* aEy = dom.aEy.data(); const real* bEy = dom.bEy.data();
    const real* aEz = dom.aEz.data(); const real* bEz = dom.bEz.data();

    real* Ez_old = dom.tls_Ez_old.data();
    const real* dPz_dt = dom.tls_dPz_dt.data();
    const bool tls_enabled = dom.tls_enabled;

    SubdomainCPML& cpml = dom.cpml;
    const real* inv_kx = cpml.inv_kx.data();
    const real* inv_ky = cpml.inv_ky.data();
    const real* inv_kz = cpml.inv_kz.data();
    const real* bEx_pml = cpml.bEx.data(); const real* cEx_pml = cpml.cEx.data();
    const real* bEy_pml = cpml.bEy.data(); const real* cEy_pml = cpml.cEy.data();
    const real* bEz_pml = cpml.bEz.data(); const real* cEz_pml = cpml.cEz.data();

    real* psi_Ex_y = cpml.psi_Ex_y.data(); real* psi_Ex_z = cpml.psi_Ex_z.data();
    real* psi_Ey_z = cpml.psi_Ey_z.data(); real* psi_Ey_x = cpml.psi_Ey_x.data();
    real* psi_Ez_x = cpml.psi_Ez_x.data(); real* psi_Ez_y = cpml.psi_Ez_y.data();

#if FDTD_OMP_ENABLED
#pragma omp parallel for
#endif
    for (int i = 1; i < (int)NxL; ++i) {
        for (int j = 1; j < (int)NyL; ++j) {
            for (int k = 1; k < (int)NzL; ++k) {
                size_t id = ((size_t)i * NyL + (size_t)j) * NzL + (size_t)k;

                const real inv_dx_i = inv_dx[i];
                const real inv_dy_j = inv_dy[j];
                const real inv_dz_k = inv_dz[k];

                real dHz_dy = (Hz[id] - Hz[id - sJ]) * inv_dy_j;
                real dHy_dz = (Hy[id] - Hy[id - sK]) * inv_dz_k;
                real dHx_dz = (Hx[id] - Hx[id - sK]) * inv_dz_k;
                real dHz_dx = (Hz[id] - Hz[id - sI]) * inv_dx_i;
                real dHy_dx = (Hy[id] - Hy[id - sI]) * inv_dx_i;
                real dHx_dy = (Hx[id] - Hx[id - sJ]) * inv_dy_j;

                // Ex
                real corr_y = dHz_dy * inv_ky[j];
                psi_Ex_y[id] = bEy_pml[j] * psi_Ex_y[id] + cEy_pml[j] * dHz_dy;
                corr_y += psi_Ex_y[id];

                real corr_z = dHy_dz * inv_kz[k];
                psi_Ex_z[id] = bEz_pml[k] * psi_Ex_z[id] + cEz_pml[k] * dHy_dz;
                corr_z += psi_Ex_z[id];

                real curlHx = corr_y - corr_z;

                // Ey
                corr_z = dHx_dz * inv_kz[k];
                psi_Ey_z[id] = bEz_pml[k] * psi_Ey_z[id] + cEz_pml[k] * dHx_dz;
                corr_z += psi_Ey_z[id];

                real corr_x = dHz_dx * inv_kx[i];
                psi_Ey_x[id] = bEx_pml[i] * psi_Ey_x[id] + cEx_pml[i] * dHz_dx;
                corr_x += psi_Ey_x[id];

                real curlHy = corr_z - corr_x;

                // Ez
                corr_x = dHy_dx * inv_kx[i];
                psi_Ez_x[id] = bEx_pml[i] * psi_Ez_x[id] + cEx_pml[i] * dHy_dx;
                corr_x += psi_Ez_x[id];

                corr_y = dHx_dy * inv_ky[j];
                psi_Ez_y[id] = bEy_pml[j] * psi_Ez_y[id] + cEy_pml[j] * dHx_dy;
                corr_y += psi_Ez_y[id];

                real curlHz = corr_x - corr_y;

                Ex[id] = aEx[id] * Ex[id] + bEx[id] * (curlHx - Jx[id]);
                Ey[id] = aEy[id] * Ey[id] + bEy[id] * (curlHy - Jy[id]);

                // Store Ez_old before update (for stimulated term calculation)
                if (tls_enabled) {
                    Ez_old[id] = Ez[id];
                }

                // Ez with TLS polarization source: -dPz/dt
                real Pz_source = tls_enabled ? dPz_dt[id] : 0.0;
                Ez[id] = aEz[id] * Ez[id] + bEz[id] * (curlHz - Jz[id] - Pz_source);
            }
        }
    }
}

inline void parallel_update_E_with_tls(DomainDecomposition& decomp, real dt, const TwoLevelParams& params) {
    if (!decomp.is_parallel()) return;

    parallel_for_domains(decomp, [&params](Subdomain& dom) {
        update_polarization_subdomain(dom, params);
    });

    parallel_for_domains(decomp, [](Subdomain& dom) {
        fdtd_update_E_with_tls_subdomain(dom);
    });

    decomp.exchange_halos();

    parallel_for_domains(decomp, [dt, &params](Subdomain& dom) {
        update_populations_subdomain(dom, dt, params);
    });
}

// Compute integrated upper state population (density * dV -> atom count)
inline real parallel_compute_integrated_Nu(const DomainDecomposition& decomp) {
    if (!decomp.is_parallel()) return 0.0;

    real total = 0.0;
    for (const auto& dom : decomp.subdomains) {
        if (!dom.tls_enabled) continue;

        size_t i_start, i_end, j_start, j_end, k_start, k_end;
        decomp.get_interior_range(dom, i_start, i_end, j_start, j_end, k_start, k_end);

        for (size_t i = i_start; i < i_end; ++i) {
            for (size_t j = j_start; j < j_end; ++j) {
                for (size_t k = k_start; k < k_end; ++k) {
                    size_t id = dom.local_idx(i, j, k);
                    if (dom.tls_Ndip[id] <= 0) continue;

                    real dV = dom.dx_local[i] * dom.dy_local[j] * dom.dz_local[k];
                    total += dom.tls_Nu[id] * dV;
                }
            }
        }
    }
    return total;
}

// Compute integrated ground state population (density * dV -> atom count)
inline real parallel_compute_integrated_Ng(const DomainDecomposition& decomp) {
    if (!decomp.is_parallel()) return 0.0;

    real total = 0.0;
    for (const auto& dom : decomp.subdomains) {
        if (!dom.tls_enabled) continue;

        size_t i_start, i_end, j_start, j_end, k_start, k_end;
        decomp.get_interior_range(dom, i_start, i_end, j_start, j_end, k_start, k_end);

        for (size_t i = i_start; i < i_end; ++i) {
            for (size_t j = j_start; j < j_end; ++j) {
                for (size_t k = k_start; k < k_end; ++k) {
                    size_t id = dom.local_idx(i, j, k);
                    if (dom.tls_Ndip[id] <= 0) continue;

                    real dV = dom.dx_local[i] * dom.dy_local[j] * dom.dz_local[k];
                    total += dom.tls_Ng[id] * dV;
                }
            }
        }
    }
    return total;
}

// Compute integrated inversion: (Nu - Ng) in atoms
inline real parallel_compute_integrated_inversion(const DomainDecomposition& decomp) {
    return parallel_compute_integrated_Nu(decomp) - parallel_compute_integrated_Ng(decomp);
}

} // namespace ParallelV2
