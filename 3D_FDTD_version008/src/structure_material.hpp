// structure_material.hpp —— Material structure package

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


// ========== Material definition ==========
struct Material {
    real eps_r{ 1.0 };   // Relative permittivity
    real mu_r{ 1.0 };    // Relative permeability
    real sigma{ 0.0 };   // Electrical conductivity S/m
    real sigma_m{ 0.0 }; // Magnetic permeability dissipation, rarely used
};

// Common helper: n,k -> eps_r (assuming mu_r=1)
inline Material make_nk(real n, real k = 0.0, real mu_r = 1.0, real sigma = 0.0) {
    // Complex refractive index n_c = n + i k, eps_r = n_c^2 = (n^2-k^2) + i(2nk)
    // Here putting loss in sigma is more intuitive; for strict treatment, use dispersion model.
    Material m;
    m.eps_r = n * n; // Simplified: ignore k -> put in sigma more appropriate (frequency dependent)
    m.mu_r = mu_r;
    m.sigma = sigma;
    return m;
}

// ========== Shapes ==========
struct AABB {
    real x0, x1, y0, y1, z0, z1;

    // Merge two AABBs
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
    // Returns "whether point (x,y,z) is inside this shape"
    virtual bool contains(real x, real y, real z) const = 0;
    // Returns the axis-aligned bounding box of this shape
    virtual AABB bounding_box() const = 0;
};

// ========== Box ==========
struct Box : public Shape {
    AABB bb;
    explicit Box(AABB a) :bb(a) {}
    bool contains(real x, real y, real z) const override {
        return (x >= bb.x0 && x < bb.x1 &&
                y >= bb.y0 && y < bb.y1 &&
                z >= bb.z0 && z < bb.z1);
    }
    AABB bounding_box() const override { return bb; }
};

// ========== Sphere ==========
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

// ========== Cylinder (z direction) ==========
struct CylinderZ : public Shape {
    real cx, cy, r, z0, z1;     // Along z axis
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

// ========== TLS Configuration for Structure ==========
// Stores TLS parameters bound to a specific structure
struct StructureTLSConfig {
    bool enabled = false;           // Is TLS enabled for this structure?
    real lambda0 = 1500e-9;         // Transition wavelength (m)
    real gamma = 7e12;              // Polarization damping rate (1/s)
    real tau = 1e-12;               // Upper level lifetime (s)
    real N0 = 1e25;                 // Dipole density (atoms/m³)
    real inversion_fraction = 1.0;  // Initial population inversion (0 to 1)
};

// Structure item = geometry + material + optional TLS
struct StructureItem {
    std::shared_ptr<Shape> shape;
    Material mat;
    StructureTLSConfig tls_config;  // TLS configuration for this structure
};

// ========== Scene container ==========
struct StructureScene {

    // Grid (including PML region)
    size_t NxT{}, NyT{}, NzT{};
    size_t Nx_core{}, Ny_core{}, Nz_core{};     // Physical region dimensions only
    size_t npml_cells{};                        // PML thickness (cells)

    // DELETED - for uniform
    // real x0{ 0 }, y0{ 0 }, z0{ 0 };
    // real dx{ 1 }, dy{ 1 }, dz{ 1 };

    GridSpacing grid_spacing;

    Material bg{};                              // Background material

    std::vector<StructureItem> items;

    // —— Utilities: coordinate/index conversion —— //
    inline std::tuple<real, real, real> cell_center(size_t i, size_t j, size_t k) const {       //MODIFIED
        //return { x0 + (i + 0.5) * dx, y0 + (j + 0.5) * dy, z0 + (k + 0.5) * dz };     //DELETED
        real x = 0.5 * (grid_spacing.x_bounds[i] + grid_spacing.x_bounds[i+1]);
        real y = 0.5 * (grid_spacing.y_bounds[j] + grid_spacing.y_bounds[j+1]);
        real z = 0.5 * (grid_spacing.z_bounds[k] + grid_spacing.z_bounds[k+1]);
        return {x, y, z};
    }

    inline std::tuple<size_t, size_t, size_t> core_center_ijk() const {
        size_t ic = npml_cells + Nx_core / 2;
        size_t jc = npml_cells + Ny_core / 2;
        size_t kc = npml_cells + Nz_core / 2;
        return { ic,jc,kc };
    }

    inline bool is_in_pml(size_t i, size_t j, size_t k) const {
        const size_t Nx = NxT, Ny = NyT, Nz = NzT;
        const size_t np = npml_cells;
        const size_t di = std::min(i, Nx - 1 - i);
        const size_t dj = std::min(j, Ny - 1 - j);
        const size_t dk = std::min(k, Nz - 1 - k);
        return (di < np) || (dj < np) || (dk < np);
    }

    // —— Add shapes with physical coordinates —— //
    // Each add_* method now accepts optional TLS configuration
    void add_box(AABB bb, const Material& m, const StructureTLSConfig& tls = {}) {
        items.push_back({ std::make_shared<Box>(bb), m, tls });
    }
    void add_sphere(real cx, real cy, real cz, real r, const Material& m, const StructureTLSConfig& tls = {}) {
        items.push_back({ std::make_shared<Sphere>(cx,cy,cz,r), m, tls });
    }
    void add_cylinder_z(real cx, real cy, real r, real z0, real z1, const Material& m, const StructureTLSConfig& tls = {}) {
        items.push_back({ std::make_shared<CylinderZ>(cx,cy,r,z0,z1), m, tls });
    }

    // —— Get structures with TLS enabled —— //
    // Returns a vector of items that have TLS enabled
    std::vector<const StructureItem*> get_tls_structures() const {
        std::vector<const StructureItem*> tls_items;
        for (const auto& item : items) {
            if (item.tls_config.enabled) {
                tls_items.push_back(&item);
            }
        }
        return tls_items;
    }

    // Check if any structure has TLS enabled
    bool has_any_tls() const {
        for (const auto& item : items) {
            if (item.tls_config.enabled) return true;
        }
        return false;
    }

    // —— Extract structure bounds for auto mesh generator —— //
    // Returns vector of {x_min, x_max, y_min, y_max, z_min, z_max, n_material}
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

    // —— Get overall domain bounds from all structures —— //
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

    // —— Add cube with "core region integer cell indices" (convenience function you want) ——
    // Pass core region center (ic_core,jc_core,kc_core) and side length (in cells); automatically accounts for npml offset
    void add_cube_by_core_cells(size_t ic_core, size_t jc_core, size_t kc_core,
        size_t side_cells, const Material& m) {
        const size_t ic = npml_cells + ic_core;
        const size_t jc = npml_cells + jc_core;
        const size_t kc = npml_cells + kc_core;
        const size_t i0 = ic - side_cells / 2, i1 = i0 + side_cells;
        const size_t j0 = jc - side_cells / 2, j1 = j0 + side_cells;
        const size_t k0 = kc - side_cells / 2, k1 = k0 + side_cells;
        
        // Use cumulative bounds
        AABB bb{
            grid_spacing.x_bounds[i0], grid_spacing.x_bounds[i1],
            grid_spacing.y_bounds[j0], grid_spacing.y_bounds[j1],
            grid_spacing.z_bounds[k0], grid_spacing.z_bounds[k1]
        };
        add_box(bb, m);
    }

    // ========= Baking: generate aE/bE/aH/bH =========
    // subpixel_samples: number of volume fraction samples (1=single point, 8=2x2x2, 27=3x3x3...)
    void bake(MaterialGrids& mg, real dt, real eps0_arg, real mu0_arg, int subpixel_samples = 1) const {
        mg.allocate(NxT, NyT, NzT);

        auto sample_mat = [&](real x, real y, real z)->Material {
            Material m = bg; // Later additions override
            for (const auto& it : items) {
                if (it.shape->contains(x, y, z)) m = it.mat;
            }
            return m;
            };

        // Sample point generation (uniform grid within cube)
        auto for_each_sub = [&](auto&& F) {
            int s = std::round(std::cbrt((double)subpixel_samples));
            if (s < 1) s = 1;
            for (int a = 0; a < s; ++a)
                for (int b = 0; b < s; ++b)
                    for (int c = 0; c < s; ++c)
                        F((a + 0.5) / s, (b + 0.5) / s, (c + 0.5) / s);
            };

        auto avg_eps_mu_sigma_at = [&](size_t i, size_t j, size_t k,
            real ux, real uy, real uz)->std::tuple<real, real, real, real> {
                // Uniform sampling in small voxel at (x0,y0,z0) + u*(dx,dy,dz)
                // ux,uy,uz ∈ {0,0.5,1} etc., for Yee stagger (see below)

                // Get local cell dimensions
                real dx_local = grid_spacing.dx[i];
                real dy_local = grid_spacing.dy[j];
                real dz_local = grid_spacing.dz[k];
                
                // Cell corner position
                real x0 = grid_spacing.x_bounds[i];
                real y0 = grid_spacing.y_bounds[j];
                real z0 = grid_spacing.z_bounds[k];

                real eps_sum = 0, mu_sum = 0, sg_sum = 0, sgm_sum = 0;
                int cnt = 0;

                for_each_sub([&](real rx, real ry, real rz) {    
                    real xs = x0 + (ux + rx * (1.0)) * dx_local;    //MODIFIED    didnt have the * (1.0) on thing?
                    real ys = y0 + (uy + ry * (1.0)) * dy_local;    //MODIFIED
                    real zs = z0 + (uz + rz * (1.0)) * dz_local;    //MODIFIED
                    auto m = sample_mat(xs, ys, zs);
                    eps_sum += m.eps_r * eps0_arg;
                    mu_sum += m.mu_r * mu0_arg;
                    sg_sum += m.sigma;
                    sgm_sum += m.sigma_m;
                    ++cnt;
                });
                return { eps_sum / cnt, mu_sum / cnt, sg_sum / cnt, sgm_sum / cnt };
            };

        // —— Strict positions of Yee grid points:
        // Ex(i,j,k): (i*dx, (j+0.5)dy, (k+0.5)dz) -> (ux,uy,uz) = (0,0.5,0.5)
        // Ey(i,j,k): ((i+0.5)dx, 0.5*dy, (k+0.5)dz) -> (0.5,0,0.5)
        // Ez(i,j,k): ((i+0.5)dx, (j+0.5)dy, 0.5*dz) -> (0.5,0.5,0)
        // Hx(i,j,k): ((i+0.5)dx, j*dy, k*dz)     -> (0.5,0,0)
        // Hy(i,j,k): (i*dx, (j+0.5)dy, k*dz)     -> (0,0.5,0)
        // Hz(i,j,k): (i*dx, j*dy, (k+0.5)dz)     -> (0,0,0.5)

        for (size_t i = 0; i < NxT; ++i)
            for (size_t j = 0; j < NyT; ++j)
                for (size_t k = 0; k < NzT; ++k) {
                    size_t id = idx3(i, j, k, NyT, NzT);
                    //real x_cell = x0 + i * dx, y_cell = y0 + j * dy, z_cell = z0 + k * dz;    DELTED? maybe? it didnt say to delete just didnt show up

                    auto do_E = [&](real ux, real uy, real uz, real& aE, real& bE) {
                        auto [eps, mu, sigma, sigm] = avg_eps_mu_sigma_at(i, j, k, ux, uy, uz);    //MODIFIED - i j k used to be x_cell y_cell z_cell
                        real c1 = sigma * dt / (2 * eps);
                        aE = (1 - c1) / (1 + c1);
                        bE = (dt / eps) / (1 + c1);
                    };
                    auto do_H = [&](real ux, real uy, real uz, real& aH, real& bH) {
                        auto [eps, mu, sigma, sigm] = avg_eps_mu_sigma_at(i, j, k, ux, uy, uz);    //MODIFIED - i j k used to be x_cell y_cell z_cell
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

// ========= Generate refractive index 3D grid from StructureScene (cell-center sampling) =========
inline void make_n_grid_from_scene(const StructureScene& scene,
    std::vector<real>& n_grid) {
    const size_t NxT = scene.NxT, NyT = scene.NyT, NzT = scene.NzT;
    n_grid.assign(NxT * NyT * NzT, 1.0);

    auto contains_mat = [&](real x, real y, real z)->Material {
        Material m = scene.bg;
        for (const auto& it : scene.items) {
            if (it.shape->contains(x, y, z)) m = it.mat;
        }
        return m;
    };

    for (size_t i = 0; i < NxT; ++i)
        for (size_t j = 0; j < NyT; ++j)
            for (size_t k = 0; k < NzT; ++k) {
                // const real xc = scene.x0 + (i + 0.5) * scene.dx;
                // const real yc = scene.y0 + (j + 0.5) * scene.dy;
                // const real zc = scene.z0 + (k + 0.5) * scene.dz;
                // CHANGED to: use cell_center helper
                auto [xc, yc, zc] = scene.cell_center(i, j, k);

                const Material m = contains_mat(xc, yc, zc);
                const real eps_r = m.eps_r;
                const real mu_r = m.mu_r;  // Usually 1
                n_grid[idx3(i, j, k, NyT, NzT)] = std::sqrt(eps_r * mu_r);
            }
}



