// global_function.hpp — Global utility functions and common definitions

#pragma once
#include <cstddef>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numbers>

using real = double;

namespace PhysConst {
    inline constexpr real PI   = std::numbers::pi_v<real>;
    inline constexpr real C0   = 299792458.0;              // Speed of light (m/s)
    inline constexpr real EPS0 = 8.854187817e-12;          // Vacuum permittivity (F/m)
    inline constexpr real MU0  = 4.0 * PI * 1e-7;          // Vacuum permeability (H/m)
    inline constexpr real HBAR = 1.054571817e-34;          // Reduced Planck constant (J·s)
    inline constexpr real Z0   = 376.730313668;            // Vacuum impedance (Ohm)
}

constexpr inline std::size_t idx3(std::size_t i, std::size_t j, std::size_t k,
    std::size_t Ny, std::size_t Nz) noexcept {
    return (i * Ny + j) * Nz + k;
}

// tau = sqrt(2 ln2) / (pi * df_fwhm) for Gaussian envelope
constexpr inline real tau_from_df_fwhm(real df_fwhm) {
    if (df_fwhm <= 0.0) return 0.0;
    return std::sqrt(2.0 * std::log(2.0)) / (3.14159265358979323846 * df_fwhm);
}

#if defined(_MSC_VER)
#define FDTD_ALWAYS_INLINE __forceinline
#elif defined(__GNUC__) || defined(__clang__)
#define FDTD_ALWAYS_INLINE inline __attribute__((always_inline))
#else
#define FDTD_ALWAYS_INLINE inline
#endif

namespace fdtd_math {

template<typename T>
FDTD_ALWAYS_INLINE T diff_x(const T* __restrict F, std::size_t id, std::size_t sI, T inv_dx) {
    return (F[id + sI] - F[id]) * inv_dx;
}
template<typename T>
FDTD_ALWAYS_INLINE T diff_y(const T* __restrict F, std::size_t id, std::size_t sJ, T inv_dy) {
    return (F[id + sJ] - F[id]) * inv_dy;
}
template<typename T>
FDTD_ALWAYS_INLINE T diff_z(const T* __restrict F, std::size_t id, std::size_t sK, T inv_dz) {
    return (F[id + sK] - F[id]) * inv_dz;
}

template<typename T>
FDTD_ALWAYS_INLINE T diff_xm(const T* __restrict F, std::size_t id, std::size_t sI, T inv_dx) {
    return (F[id] - F[id - sI]) * inv_dx;
}
template<typename T>
FDTD_ALWAYS_INLINE T diff_ym(const T* __restrict F, std::size_t id, std::size_t sJ, T inv_dy) {
    return (F[id] - F[id - sJ]) * inv_dy;
}
template<typename T>
FDTD_ALWAYS_INLINE T diff_zm(const T* __restrict F, std::size_t id, std::size_t sK, T inv_dz) {
    return (F[id] - F[id - sK]) * inv_dz;
}

}

struct MaterialGrids {
    size_t NxT{}, NyT{}, NzT{};

    std::vector<real> aEx, bEx;
    std::vector<real> aEy, bEy;
    std::vector<real> aEz, bEz;
    std::vector<real> aHx, bHx;
    std::vector<real> aHy, bHy;
    std::vector<real> aHz, bHz;

    void allocate(size_t Nx, size_t Ny, size_t Nz) {
        NxT = Nx; NyT = Ny; NzT = Nz;
        const size_t N = NxT * NyT * NzT;
        aEx.assign(N, 1); bEx.assign(N, 0);
        aEy.assign(N, 1); bEy.assign(N, 0);
        aEz.assign(N, 1); bEz.assign(N, 0);
        aHx.assign(N, 1); bHx.assign(N, 0);
        aHy.assign(N, 1); bHy.assign(N, 0);
        aHz.assign(N, 1); bHz.assign(N, 0);
    }
};

struct GridSpacing {
    std::vector<real> dx, dy, dz;
    std::vector<real> inv_dx, inv_dy, inv_dz;
    std::vector<real> x_bounds, y_bounds, z_bounds;
    size_t npml = 0;

    void compute_inverses() {
        inv_dx.resize(dx.size());
        inv_dy.resize(dy.size());
        inv_dz.resize(dz.size());
        for (size_t i = 0; i < dx.size(); ++i) inv_dx[i] = 1.0 / dx[i];
        for (size_t j = 0; j < dy.size(); ++j) inv_dy[j] = 1.0 / dy[j];
        for (size_t k = 0; k < dz.size(); ++k) inv_dz[k] = 1.0 / dz[k];
    }

    real get_x(size_t i) const { return x_bounds[i]; }
    real get_y(size_t j) const { return y_bounds[j]; }
    real get_z(size_t k) const { return z_bounds[k]; }

    real pml_offset_x() const { return (npml < x_bounds.size()) ? x_bounds[npml] : 0.0; }
    real pml_offset_y() const { return (npml < y_bounds.size()) ? y_bounds[npml] : 0.0; }
    real pml_offset_z() const { return (npml < z_bounds.size()) ? z_bounds[npml] : 0.0; }

    real to_physical_x(real total_x) const { return total_x - pml_offset_x(); }
    real to_physical_y(real total_y) const { return total_y - pml_offset_y(); }
    real to_physical_z(real total_z) const { return total_z - pml_offset_z(); }

    real to_total_x(real phys_x) const { return phys_x + pml_offset_x(); }
    real to_total_y(real phys_y) const { return phys_y + pml_offset_y(); }
    real to_total_z(real phys_z) const { return phys_z + pml_offset_z(); }

    real cell_center_physical_x(size_t i) const {
        if (i + 1 >= x_bounds.size()) return to_physical_x(x_bounds.back());
        return to_physical_x(0.5 * (x_bounds[i] + x_bounds[i + 1]));
    }
    real cell_center_physical_y(size_t j) const {
        if (j + 1 >= y_bounds.size()) return to_physical_y(y_bounds.back());
        return to_physical_y(0.5 * (y_bounds[j] + y_bounds[j + 1]));
    }
    real cell_center_physical_z(size_t k) const {
        if (k + 1 >= z_bounds.size()) return to_physical_z(z_bounds.back());
        return to_physical_z(0.5 * (z_bounds[k] + z_bounds[k + 1]));
    }

    // Uses binary search to find the cell containing the given position
    size_t physical_to_index_x(real phys_x) const {
        real total_x = to_total_x(phys_x);
        if (total_x <= x_bounds.front()) return 0;
        if (total_x >= x_bounds.back()) return x_bounds.size() - 2;
        auto it = std::upper_bound(x_bounds.begin(), x_bounds.end(), total_x);
        return (it - x_bounds.begin()) - 1;
    }
    size_t physical_to_index_y(real phys_y) const {
        real total_y = to_total_y(phys_y);
        if (total_y <= y_bounds.front()) return 0;
        if (total_y >= y_bounds.back()) return y_bounds.size() - 2;
        auto it = std::upper_bound(y_bounds.begin(), y_bounds.end(), total_y);
        return (it - y_bounds.begin()) - 1;
    }
    size_t physical_to_index_z(real phys_z) const {
        real total_z = to_total_z(phys_z);
        if (total_z <= z_bounds.front()) return 0;
        if (total_z >= z_bounds.back()) return z_bounds.size() - 2;
        auto it = std::upper_bound(z_bounds.begin(), z_bounds.end(), total_z);
        return (it - z_bounds.begin()) - 1;
    }
};

inline void build_cumulative_bounds(const std::vector<real>& spacing, std::vector<real>& bounds) {
    size_t N = spacing.size();
    bounds.resize(N + 1);
    bounds[0] = 0.0;
    for (size_t i = 0; i < N; ++i) {
        bounds[i + 1] = bounds[i] + spacing[i];
    }
}

namespace NumericUtils {

template<typename T>
inline bool has_nan_or_inf(const std::vector<T>& v) {
    for (const auto& val : v) {
        if (std::isnan(val) || std::isinf(val)) return true;
    }
    return false;
}

template<typename T>
inline T max_abs(const std::vector<T>& v) {
    T max_val = 0;
    for (const auto& val : v) {
        max_val = std::max(max_val, std::abs(val));
    }
    return max_val;
}

template<typename T>
inline T min_value(const std::vector<T>& v) {
    if (v.empty()) return T{};
    return *std::min_element(v.begin(), v.end());
}

template<typename T>
inline T max_value(const std::vector<T>& v) {
    if (v.empty()) return T{};
    return *std::max_element(v.begin(), v.end());
}

inline real safe_div(real numerator, real denominator, real threshold = 1e-30) {
    return (std::abs(denominator) > threshold) ? (numerator / denominator) : 0.0;
}

template<typename T>
inline T clamp(T val, T lo, T hi) {
    return std::max(lo, std::min(val, hi));
}

inline real lerp(real a, real b, real t) {
    return a + t * (b - a);
}

}

namespace CoordUtils {

inline size_t pos_to_index(real pos, const std::vector<real>& bounds) {
    if (pos <= bounds.front()) return 0;
    if (pos >= bounds.back()) return bounds.size() - 2;
    auto it = std::upper_bound(bounds.begin(), bounds.end(), pos);
    return (it - bounds.begin()) - 1;
}

inline real index_to_center(size_t idx, const std::vector<real>& bounds) {
    if (idx + 1 >= bounds.size()) return bounds.back();
    return 0.5 * (bounds[idx] + bounds[idx + 1]);
}

inline bool is_in_pml(size_t i, size_t j, size_t k,
                      size_t NxT, size_t NyT, size_t NzT, size_t npml) {
    size_t di = std::min(i, NxT - 1 - i);
    size_t dj = std::min(j, NyT - 1 - j);
    size_t dk = std::min(k, NzT - 1 - k);
    return (di < npml) || (dj < npml) || (dk < npml);
}

inline size_t safe_sub(size_t a, size_t b) {
    return (a >= b) ? (a - b) : 0;
}

}

namespace UnitConv {

inline real lambda_to_freq(real lambda, real c0 = PhysConst::C0) { return c0 / lambda; }
inline real freq_to_lambda(real freq, real c0 = PhysConst::C0) { return c0 / freq; }
inline real freq_to_omega(real freq) { return 2.0 * PhysConst::PI * freq; }
inline real lambda_to_omega(real lambda, real c0 = PhysConst::C0) { return 2.0 * PhysConst::PI * c0 / lambda; }

inline constexpr real nm_to_m(real nm) { return nm * 1e-9; }
inline constexpr real m_to_nm(real m) { return m * 1e9; }
inline constexpr real fs_to_s(real fs) { return fs * 1e-15; }
inline constexpr real s_to_fs(real s) { return s * 1e15; }

}
