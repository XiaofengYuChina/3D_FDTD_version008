// structure_material.hpp - Material structure package

#pragma once

#include <vector>
#include <memory>
#include <functional>
#include <cmath>
#include <string>
#include <tuple>
#include <unordered_map>

#include "global_function.hpp"
#include "boundary.hpp"

struct Material {
    real eps_r{ 1.0 };    // Relative permittivity
    real mu_r{ 1.0 };     // Relative permeability
    real sigma{ 0.0 };    // Electrical conductivity (S/m)
    real sigma_m{ 0.0 };  // Magnetic loss (rarely used)
};

// Create material from refractive index n (simplified: loss via sigma, not k)
inline Material make_nk(real n, real k = 0.0, real mu_r = 1.0, real sigma = 0.0) {
    Material m;
    m.eps_r = n * n;
    m.mu_r = mu_r;
    m.sigma = sigma;
    return m;
}

struct AABB {
    real x0, x1, y0, y1, z0, z1;

    AABB merge(const AABB& other) const {
        return {
            std::min(x0, other.x0), std::max(x1, other.x1),
            std::min(y0, other.y0), std::max(y1, other.y1),
            std::min(z0, other.z0), std::max(z1, other.z1)
        };
    }
};

struct Shape {
    virtual ~Shape() = default;
    virtual bool contains(real x, real y, real z) const = 0;
    virtual AABB bounding_box() const = 0;
};

struct Box : public Shape {
    AABB bb;
    explicit Box(AABB a) : bb(a) {}
    bool contains(real x, real y, real z) const override {
        return (x >= bb.x0 && x < bb.x1 &&
                y >= bb.y0 && y < bb.y1 &&
                z >= bb.z0 && z < bb.z1);
    }
    AABB bounding_box() const override { return bb; }
};

struct Sphere : public Shape {
    real cx, cy, cz, r;
    Sphere(real cx_, real cy_, real cz_, real r_) : cx(cx_), cy(cy_), cz(cz_), r(r_) {}
    bool contains(real x, real y, real z) const override {
        real dx = x - cx, dy = y - cy, dz = z - cz;
        return (dx * dx + dy * dy + dz * dz) <= r * r;
    }
    AABB bounding_box() const override {
        return {cx - r, cx + r, cy - r, cy + r, cz - r, cz + r};
    }
};

struct CylinderZ : public Shape {
    real cx, cy, r, z0, z1;
    CylinderZ(real cx_, real cy_, real r_, real z0_, real z1_)
        : cx(cx_), cy(cy_), r(r_), z0(z0_), z1(z1_) {}
    bool contains(real x, real y, real z) const override {
        if (z < z0 || z >= z1) return false;
        real dx = x - cx, dy = y - cy;
        return (dx * dx + dy * dy) <= r * r;
    }
    AABB bounding_box() const override {
        return {cx - r, cx + r, cy - r, cy + r, z0, z1};
    }
};

struct StructureTLSConfig {
    bool enabled = false;
    real lambda0 = 1500e-9;         // Transition wavelength (m)
    real gamma = 7e12;              // Polarization damping rate (1/s)
    real tau = 1e-12;               // Upper level lifetime (s)
    real N0 = 1e25;                 // Dipole density (atoms/m^3)
    real inversion_fraction = 1.0;  // Initial population inversion (0 to 1)
};

struct StructureItem {
    std::shared_ptr<Shape> shape;
    Material mat;
    StructureTLSConfig tls_config;
};

struct StructureScene {
    size_t NxT{}, NyT{}, NzT{};
    size_t Nx_core{}, Ny_core{}, Nz_core{};  // Physical region dimensions
    size_t npml_cells{};                      // PML thickness (cells)

    GridSpacing grid_spacing;
    Material bg{};

    std::vector<StructureItem> items;

    inline std::tuple<real, real, real> cell_center(size_t i, size_t j, size_t k) const {
        real x = 0.5 * (grid_spacing.x_bounds[i] + grid_spacing.x_bounds[i+1]);
        real y = 0.5 * (grid_spacing.y_bounds[j] + grid_spacing.y_bounds[j+1]);
        real z = 0.5 * (grid_spacing.z_bounds[k] + grid_spacing.z_bounds[k+1]);
        return {x, y, z};
    }

    inline std::tuple<size_t, size_t, size_t> core_center_ijk() const {
        size_t ic = npml_cells + Nx_core / 2;
        size_t jc = npml_cells + Ny_core / 2;
        size_t kc = npml_cells + Nz_core / 2;
        return {ic, jc, kc};
    }

    inline bool is_in_pml(size_t i, size_t j, size_t k) const {
        const size_t Nx = NxT, Ny = NyT, Nz = NzT;
        const size_t np = npml_cells;
        const size_t di = std::min(i, Nx - 1 - i);
        const size_t dj = std::min(j, Ny - 1 - j);
        const size_t dk = std::min(k, Nz - 1 - k);
        return (di < np) || (dj < np) || (dk < np);
    }

    void add_box(AABB bb, const Material& m, const StructureTLSConfig& tls = {}) {
        items.push_back({std::make_shared<Box>(bb), m, tls});
    }

    void add_sphere(real cx, real cy, real cz, real r, const Material& m, const StructureTLSConfig& tls = {}) {
        items.push_back({std::make_shared<Sphere>(cx, cy, cz, r), m, tls});
    }

    void add_cylinder_z(real cx, real cy, real r, real z0, real z1, const Material& m, const StructureTLSConfig& tls = {}) {
        items.push_back({std::make_shared<CylinderZ>(cx, cy, r, z0, z1), m, tls});
    }

    std::vector<const StructureItem*> get_tls_structures() const {
        std::vector<const StructureItem*> tls_items;
        for (const auto& item : items) {
            if (item.tls_config.enabled) {
                tls_items.push_back(&item);
            }
        }
        return tls_items;
    }

    bool has_any_tls() const {
        for (const auto& item : items) {
            if (item.tls_config.enabled) return true;
        }
        return false;
    }

    struct StructBounds {
        real x_min, x_max;
        real y_min, y_max;
        real z_min, z_max;
        real n_material;
    };

    std::vector<StructBounds> get_structure_bounds() const {
        std::vector<StructBounds> bounds;
        for (const auto& item : items) {
            AABB bb = item.shape->bounding_box();
            real n = std::sqrt(item.mat.eps_r * item.mat.mu_r);
            bounds.push_back({bb.x0, bb.x1, bb.y0, bb.y1, bb.z0, bb.z1, n});
        }
        return bounds;
    }

    AABB get_total_bounds() const {
        if (items.empty()) {
            return {0, 0, 0, 0, 0, 0};
        }
        AABB total = items[0].shape->bounding_box();
        for (size_t i = 1; i < items.size(); ++i) {
            total = total.merge(items[i].shape->bounding_box());
        }
        return total;
    }

    // Add cube using core region cell indices (auto-adjusts for PML offset)
    void add_cube_by_core_cells(size_t ic_core, size_t jc_core, size_t kc_core,
                                 size_t side_cells, const Material& m) {
        const size_t ic = npml_cells + ic_core;
        const size_t jc = npml_cells + jc_core;
        const size_t kc = npml_cells + kc_core;
        const size_t i0 = ic - side_cells / 2, i1 = i0 + side_cells;
        const size_t j0 = jc - side_cells / 2, j1 = j0 + side_cells;
        const size_t k0 = kc - side_cells / 2, k1 = k0 + side_cells;

        AABB bb{
            grid_spacing.x_bounds[i0], grid_spacing.x_bounds[i1],
            grid_spacing.y_bounds[j0], grid_spacing.y_bounds[j1],
            grid_spacing.z_bounds[k0], grid_spacing.z_bounds[k1]
        };
        add_box(bb, m);
    }

    // subpixel_samples: volume fraction samples (1=single point, 8=2x2x2, 27=3x3x3...)
    void bake(MaterialGrids& mg, real dt, real eps0_arg, real mu0_arg, int subpixel_samples = 1) const {
        mg.allocate(NxT, NyT, NzT);

        auto sample_mat = [&](real x, real y, real z) -> Material {
            Material m = bg;  // Later additions override
            for (const auto& it : items) {
                if (it.shape->contains(x, y, z)) m = it.mat;
            }
            return m;
        };

        auto for_each_sub = [&](auto&& F) {
            int s = std::round(std::cbrt((double)subpixel_samples));
            if (s < 1) s = 1;
            for (int a = 0; a < s; ++a)
                for (int b = 0; b < s; ++b)
                    for (int c = 0; c < s; ++c)
                        F((a + 0.5) / s, (b + 0.5) / s, (c + 0.5) / s);
        };

        auto avg_eps_mu_sigma_at = [&](size_t i, size_t j, size_t k,
                                        real ux, real uy, real uz) -> std::tuple<real, real, real, real> {
            real dx_local = grid_spacing.dx[i];
            real dy_local = grid_spacing.dy[j];
            real dz_local = grid_spacing.dz[k];
            real x0 = grid_spacing.x_bounds[i];
            real y0 = grid_spacing.y_bounds[j];
            real z0 = grid_spacing.z_bounds[k];

            real eps_sum = 0, mu_sum = 0, sg_sum = 0, sgm_sum = 0;
            int cnt = 0;

            for_each_sub([&](real rx, real ry, real rz) {
                real xs = x0 + (ux + rx) * dx_local;
                real ys = y0 + (uy + ry) * dy_local;
                real zs = z0 + (uz + rz) * dz_local;
                auto m = sample_mat(xs, ys, zs);
                eps_sum += m.eps_r * eps0_arg;
                mu_sum += m.mu_r * mu0_arg;
                sg_sum += m.sigma;
                sgm_sum += m.sigma_m;
                ++cnt;
            });
            return {eps_sum / cnt, mu_sum / cnt, sg_sum / cnt, sgm_sum / cnt};
        };

        // Yee grid stagger offsets (ux,uy,uz):
        // Ex: (0, 0.5, 0.5), Ey: (0.5, 0, 0.5), Ez: (0.5, 0.5, 0)
        // Hx: (0.5, 0, 0),   Hy: (0, 0.5, 0),   Hz: (0, 0, 0.5)

        for (size_t i = 0; i < NxT; ++i)
            for (size_t j = 0; j < NyT; ++j)
                for (size_t k = 0; k < NzT; ++k) {
                    size_t id = idx3(i, j, k, NyT, NzT);

                    auto do_E = [&](real ux, real uy, real uz, real& aE, real& bE) {
                        auto [eps, mu, sigma, sigm] = avg_eps_mu_sigma_at(i, j, k, ux, uy, uz);
                        real c1 = sigma * dt / (2 * eps);
                        aE = (1 - c1) / (1 + c1);
                        bE = (dt / eps) / (1 + c1);
                    };

                    auto do_H = [&](real ux, real uy, real uz, real& aH, real& bH) {
                        auto [eps, mu, sigma, sigm] = avg_eps_mu_sigma_at(i, j, k, ux, uy, uz);
                        real c2 = sigm * dt / (2 * mu);
                        aH = (1 - c2) / (1 + c2);
                        bH = (dt / mu) / (1 + c2);
                    };

                    do_E(0.0, 0.5, 0.5, mg.aEx[id], mg.bEx[id]);
                    do_E(0.5, 0.0, 0.5, mg.aEy[id], mg.bEy[id]);
                    do_E(0.5, 0.5, 0.0, mg.aEz[id], mg.bEz[id]);

                    do_H(0.5, 0.0, 0.0, mg.aHx[id], mg.bHx[id]);
                    do_H(0.0, 0.5, 0.0, mg.aHy[id], mg.bHy[id]);
                    do_H(0.0, 0.0, 0.5, mg.aHz[id], mg.bHz[id]);
                }
    }
};

inline void make_n_grid_from_scene(const StructureScene& scene, std::vector<real>& n_grid) {
    const size_t NxT = scene.NxT, NyT = scene.NyT, NzT = scene.NzT;
    n_grid.assign(NxT * NyT * NzT, 1.0);

    auto contains_mat = [&](real x, real y, real z) -> Material {
        Material m = scene.bg;
        for (const auto& it : scene.items) {
            if (it.shape->contains(x, y, z)) m = it.mat;
        }
        return m;
    };

    for (size_t i = 0; i < NxT; ++i)
        for (size_t j = 0; j < NyT; ++j)
            for (size_t k = 0; k < NzT; ++k) {
                auto [xc, yc, zc] = scene.cell_center(i, j, k);
                const Material m = contains_mat(xc, yc, zc);
                n_grid[idx3(i, j, k, NyT, NzT)] = std::sqrt(m.eps_r * m.mu_r);
            }
}
