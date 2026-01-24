// boundary.hpp —— Boundary conditions

#pragma once

#include <vector>
#include <cmath>
#include <cstddef>
#include <array>
#include <string>
#include <memory>
#include <algorithm>
#include <iostream>

#include "global_function.hpp"

// === Only two types ===
enum class BcType { PEC, CPML_RC };

// === CPML common configuration ===
struct CPMLConfig {
    int   npml = 8;            // Thickness (cells)
    real  m = 3.0;             // σ/κ polynomial order
    real  Rerr = 1e-10;        // Target reflection
    real  alpha0 = 0.05;       // α0 ~ (c0/Δs)*0~0.2
    real  kappa_max = 8.0;     // κ_max 6~12
    bool  alpha_linear = true; // α linear distribution
};

struct BoundaryParams {
    BcType type = BcType::CPML_RC;
    CPMLConfig cpml;
};

// === Boundary interface ===
struct IBoundary {
    virtual ~IBoundary() = default;

    // PML correction called after main core update
    virtual void apply_after_H(
        std::vector<real>& Ex, std::vector<real>& Ey, std::vector<real>& Ez,
        std::vector<real>& Hx, std::vector<real>& Hy, std::vector<real>& Hz) {
        (void)Ex; (void)Ey; (void)Ez; (void)Hx; (void)Hy; (void)Hz;
    }

    virtual void apply_after_E(
        std::vector<real>& Ex, std::vector<real>& Ey, std::vector<real>& Ez,
        std::vector<real>& Hx, std::vector<real>& Hy, std::vector<real>& Hz) = 0;

    // To make PML increment scaling consistent with main core, bind MaterialGrids' bE/bH (optional)
    virtual void bind_material_coeffs(
        const std::vector<real>*, const std::vector<real>*, const std::vector<real>*,
        const std::vector<real>*, const std::vector<real>*, const std::vector<real>*) {
    }

    // Dimensions
    virtual size_t NxT() const = 0;
    virtual size_t NyT() const = 0;
    virtual size_t NzT() const = 0;
    virtual size_t npml() const = 0;
};

// === PEC boundary ===
struct PECBoundary final : public IBoundary {
    size_t Nx_, Ny_, Nz_;
    explicit PECBoundary(size_t Nx, size_t Ny, size_t Nz) : Nx_(Nx), Ny_(Ny), Nz_(Nz) {}

    void apply_after_E(
        std::vector<real>& Ex, std::vector<real>& Ey, std::vector<real>& Ez,
        std::vector<real>& Hx, std::vector<real>& Hy, std::vector<real>& Hz) override {
        (void)Hx; (void)Hy; (void)Hz;
        for (size_t j = 0; j < Ny_; ++j)
            for (size_t k = 0; k < Nz_; ++k) {
            Ex[idx3(0, j, k, Ny_, Nz_)] = Ey[idx3(0, j, k, Ny_, Nz_)] = Ez[idx3(0, j, k, Ny_, Nz_)] = 0.0;
            Ex[idx3(Nx_ - 1, j, k, Ny_, Nz_)] = Ey[idx3(Nx_ - 1, j, k, Ny_, Nz_)] = Ez[idx3(Nx_ - 1, j, k, Ny_, Nz_)] = 0.0;
        }
        for (size_t i = 0; i < Nx_; ++i)
            for (size_t k = 0; k < Nz_; ++k) {
            Ex[idx3(i, 0, k, Ny_, Nz_)] = Ey[idx3(i, 0, k, Ny_, Nz_)] = Ez[idx3(i, 0, k, Ny_, Nz_)] = 0.0;
            Ex[idx3(i, Ny_ - 1, k, Ny_, Nz_)] = Ey[idx3(i, Ny_ - 1, k, Ny_, Nz_)] = Ez[idx3(i, Ny_ - 1, k, Ny_, Nz_)] = 0.0;
        }
        for (size_t i = 0; i < Nx_; ++i)
            for (size_t j = 0; j < Ny_; ++j) {
            Ex[idx3(i, j, 0, Ny_, Nz_)] = Ey[idx3(i, j, 0, Ny_, Nz_)] = Ez[idx3(i, j, 0, Ny_, Nz_)] = 0.0;
            Ex[idx3(i, j, Nz_ - 1, Ny_, Nz_)] = Ey[idx3(i, j, Nz_ - 1, Ny_, Nz_)] = Ez[idx3(i, j, Nz_ - 1, Ny_, Nz_)] = 0.0;
        }
    }
    size_t NxT() const override { return Nx_; }
    size_t NyT() const override { return Ny_; }
    size_t NzT() const override { return Nz_; }
    size_t npml() const override { return 0; }
};

// === RC-CPML boundary ===
struct CPMLBoundary final : public IBoundary {

    // Core & total dimensions
    size_t NxC_, NyC_, NzC_;
    int    npml_;
    size_t NxT_, NyT_, NzT_;

    // Dimensions/constants
    real /*dx_, dy_, dz_,*/ dt_, eps0_, mu0_, c0_;
    // Spacing arrays
    std::vector<real> dx_array_, dy_array_, dz_array_;
    
    // 1D profiles
    std::vector<real> kx_, ky_, kz_;
    std::vector<real> sx_, sy_, sz_;
    std::vector<real> ax_, ay_, az_;

    // Single-pole coefficients and ψ (6 derivatives each for E/H, total 12 ψ)
    std::vector<real> bEx_, cEx_, bEy_, cEy_, bEz_, cEz_;
    std::vector<real> bHx_, cHx_, bHy_, cHy_, bHz_, cHz_;

    std::vector<real> psi_Ex_y_, psi_Ex_z_;
    std::vector<real> psi_Ey_z_, psi_Ey_x_;
    std::vector<real> psi_Ez_x_, psi_Ez_y_;
    std::vector<real> psi_Hx_y_, psi_Hx_z_;
    std::vector<real> psi_Hy_z_, psi_Hy_x_;
    std::vector<real> psi_Hz_x_, psi_Hz_y_;

    // Bind main core bE/bH (nullable, falls back to dt/eps0 or dt/mu0 if null)
    const std::vector<real>* bEx_mg_{ nullptr };
    const std::vector<real>* bEy_mg_{ nullptr };
    const std::vector<real>* bEz_mg_{ nullptr };
    const std::vector<real>* bHx_mg_{ nullptr };
    const std::vector<real>* bHy_mg_{ nullptr };
    const std::vector<real>* bHz_mg_{ nullptr };

    CPMLBoundary(size_t NxCore, size_t NyCore, size_t NzCore,
        // real dx, real dy, real dz, 
        const GridSpacing& grid_spacing,  // CHANGED: Pass full spacing
        real dt, real eps0, real mu0, real c0,
        const CPMLConfig& cfg)

        : NxC_(NxCore), NyC_(NyCore), NzC_(NzCore),
        npml_(cfg.npml),
        NxT_(NxCore + 2 * cfg.npml),
        NyT_(NyCore + 2 * cfg.npml),
        NzT_(NzCore + 2 * cfg.npml),
        // dx_(dx), dy_(dy), dz_(dz), deleted
        dt_(dt), eps0_(eps0), mu0_(mu0), c0_(c0)
    {

        // Store spacing
        dx_array_ = grid_spacing.dx;
        dy_array_ = grid_spacing.dy;
        dz_array_ = grid_spacing.dz;

        const size_t N = NxT_ * NyT_ * NzT_;
        kx_.assign(NxT_, 1.0); ky_.assign(NyT_, 1.0); kz_.assign(NzT_, 1.0);
        sx_.assign(NxT_, 0.0); sy_.assign(NyT_, 0.0); sz_.assign(NzT_, 0.0);
        ax_.assign(NxT_, 0.0); ay_.assign(NyT_, 0.0); az_.assign(NzT_, 0.0);

        bEx_.assign(NxT_, 1.0); cEx_.assign(NxT_, 0.0);
        bEy_.assign(NyT_, 1.0); cEy_.assign(NyT_, 0.0);
        bEz_.assign(NzT_, 1.0); cEz_.assign(NzT_, 0.0);
        bHx_.assign(NxT_, 1.0); cHx_.assign(NxT_, 0.0);
        bHy_.assign(NyT_, 1.0); cHy_.assign(NyT_, 0.0);
        bHz_.assign(NzT_, 1.0); cHz_.assign(NzT_, 0.0);

        psi_Ex_y_.assign(N, 0.0); psi_Ex_z_.assign(N, 0.0);
        psi_Ey_z_.assign(N, 0.0); psi_Ey_x_.assign(N, 0.0);
        psi_Ez_x_.assign(N, 0.0); psi_Ez_y_.assign(N, 0.0);
        psi_Hx_y_.assign(N, 0.0); psi_Hx_z_.assign(N, 0.0);
        psi_Hy_z_.assign(N, 0.0); psi_Hy_x_.assign(N, 0.0);
        psi_Hz_x_.assign(N, 0.0); psi_Hz_y_.assign(N, 0.0);

        build_profiles_and_coeffs(cfg);

    }

    void bind_material_coeffs(
        const std::vector<real>* bEx, const std::vector<real>* bEy, const std::vector<real>* bEz,
        const std::vector<real>* bHx, const std::vector<real>* bHy, const std::vector<real>* bHz) override {
        bEx_mg_ = bEx; bEy_mg_ = bEy; bEz_mg_ = bEz;
        bHx_mg_ = bHx; bHy_mg_ = bHy; bHz_mg_ = bHz;
    }

    size_t NxT() const override { return NxT_; }
    size_t NyT() const override { return NyT_; }
    size_t NzT() const override { return NzT_; }
    size_t npml() const override { return size_t(npml_); }

    // === Correction after H ===
    void apply_after_H(
        std::vector<real>& Ex, std::vector<real>& Ey, std::vector<real>& Ez,
        std::vector<real>& Hx, std::vector<real>& Hy, std::vector<real>& Hz) override
    {
        const size_t NxT = NxT_, NyT = NyT_, NzT = NzT_;

        auto in_pml_x = [&](int i) { return std::min(i, int(NxT) - 1 - i) < npml_; };
        auto in_pml_y = [&](int j) { return std::min(j, int(NyT) - 1 - j) < npml_; };
        auto in_pml_z = [&](int k) { return std::min(k, int(NzT) - 1 - k) < npml_; };

        using namespace fdtd_math;
        const size_t sI = NyT * NzT;
        const size_t sJ = NzT;
        const size_t sK = 1;

        const real* Exp = Ex.data();
        const real* Eyp = Ey.data();
        const real* Ezp = Ez.data();
        const real* Hxp = Hx.data();
        const real* Hyp = Hy.data();
        const real* Hzp = Hz.data();

        for (size_t i = 0; i < NxT - 1; ++i)
            for (size_t j = 0; j < NyT - 1; ++j)
                for (size_t k = 0; k < NzT - 1; ++k) {

                    size_t id = idx3(i, j, k, NyT, NzT);

                    const real inv_dx = 1.0 / dx_array_[i];
                    const real inv_dy = 1.0 / dy_array_[j];
                    const real inv_dz = 1.0 / dz_array_[k];

                    // Hx: dEz/dy and dEy/dz
                    if (in_pml_y(j)) {
                        real d = diff_y(Ezp, id, sJ, inv_dy);
                        real corr = d / ky_[j];
                        real& p = psi_Hx_y_[id]; p = bHy_[j] * p + cHy_[j] * d; corr += p;
                        const real bh = bHx_mg_ ? (*bHx_mg_)[id] : (dt_ / mu0_);
                        Hx[id] -= bh * (corr - d);
                    }
                    if (in_pml_z(k)) {
                        real d = diff_z(Eyp, id, sK, inv_dz);
                        real corr = d / kz_[k];
                        real& p = psi_Hx_z_[id]; p = bHz_[k] * p + cHz_[k] * d; corr += p;
                        const real bh = bHx_mg_ ? (*bHx_mg_)[id] : (dt_ / mu0_);
                        Hx[id] -= bh * (-(corr - d));
                    }

                    // Hy: dEx/dz and dEz/dx
                    if (in_pml_z(k)) {
                        real d = diff_z(Exp, id, sK, inv_dz);
                        real corr = d / kz_[k];
                        real& p = psi_Hy_z_[id]; p = bHz_[k] * p + cHz_[k] * d; corr += p;
                        const real bh = bHy_mg_ ? (*bHy_mg_)[id] : (dt_ / mu0_);
                        Hy[id] -= bh * (corr - d);
                    }
                    if (in_pml_x(i)) {
                        real d = diff_x(Ezp, id, sI, inv_dx);
                        real corr = d / kx_[i];
                        real& p = psi_Hy_x_[id]; p = bHx_[i] * p + cHx_[i] * d; corr += p;
                        const real bh = bHy_mg_ ? (*bHy_mg_)[id] : (dt_ / mu0_);
                        Hy[id] -= bh * (-(corr - d));
                    }

                    // Hz: dEy/dx and dEx/dy
                    if (in_pml_x(i)) {
                        real d = diff_x(Eyp, id, sI, inv_dx);
                        real corr = d / kx_[i];
                        real& p = psi_Hz_x_[id]; p = bHx_[i] * p + cHx_[i] * d; corr += p;
                        const real bh = bHz_mg_ ? (*bHz_mg_)[id] : (dt_ / mu0_);
                        Hz[id] -= bh * (corr - d);
                    }
                    if (in_pml_y(j)) {
                        real d = diff_y(Exp, id, sJ, inv_dy);
                        real corr = d / ky_[j];
                        real& p = psi_Hz_y_[id]; p = bHy_[j] * p + cHy_[j] * d; corr += p;
                        const real bh = bHz_mg_ ? (*bHz_mg_)[id] : (dt_ / mu0_);
                        Hz[id] -= bh * (-(corr - d));
                    }
                }
    }

    // === Correction after E ===
    void apply_after_E(
        std::vector<real>& Ex, std::vector<real>& Ey, std::vector<real>& Ez,
        std::vector<real>& Hx, std::vector<real>& Hy, std::vector<real>& Hz) override
    {
        const size_t NxT = NxT_, NyT = NyT_, NzT = NzT_;
        auto in_pml_x = [&](int i) { return std::min(i, int(NxT) - 1 - i) < npml_; };
        auto in_pml_y = [&](int j) { return std::min(j, int(NyT) - 1 - j) < npml_; };
        auto in_pml_z = [&](int k) { return std::min(k, int(NzT) - 1 - k) < npml_; };

        using namespace fdtd_math;
        const size_t sI = NyT * NzT;
        const size_t sJ = NzT;
        const size_t sK = 1;

        const real* Exp = Ex.data();
        const real* Eyp = Ey.data();
        const real* Ezp = Ez.data();
        const real* Hxp = Hx.data();
        const real* Hyp = Hy.data();
        const real* Hzp = Hz.data();

        for (size_t i = 1; i < NxT; ++i)
            for (size_t j = 1; j < NyT; ++j)
                for (size_t k = 1; k < NzT; ++k) {
                    size_t id = idx3(i, j, k, NyT, NzT);

                    const real inv_dx = 1.0 / dx_array_[i];
                    const real inv_dy = 1.0 / dy_array_[j];
                    const real inv_dz = 1.0 / dz_array_[k];

                    // Ex: dHz/dy, dHy/dz
                    if (in_pml_y(j)) {
                        real d = diff_ym(Hzp, id, sJ, inv_dy);
                        real corr = d / ky_[j];
                        real& p = psi_Ex_y_[id]; p = bEy_[j] * p + cEy_[j] * d; corr += p;
                        const real be = bEx_mg_ ? (*bEx_mg_)[id] : (dt_ / eps0_);
                        Ex[id] += be * (corr - d);
                    }
                    if (in_pml_z(k)) {
                        real d = diff_zm(Hyp, id, sK, inv_dz);
                        real corr = d / kz_[k];
                        real& p = psi_Ex_z_[id]; p = bEz_[k] * p + cEz_[k] * d; corr += p;
                        const real be = bEx_mg_ ? (*bEx_mg_)[id] : (dt_ / eps0_);
                        Ex[id] += be * (-(corr - d));
                    }

                    // Ey: dHx/dz, dHz/dx
                    if (in_pml_z(k)) {
                        real d = diff_zm(Hxp, id, sK, inv_dz);
                        real corr = d / kz_[k];
                        real& p = psi_Ey_z_[id]; p = bEz_[k] * p + cEz_[k] * d; corr += p;
                        const real be = bEy_mg_ ? (*bEy_mg_)[id] : (dt_ / eps0_);
                        Ey[id] += be * (corr - d);
                    }
                    if (in_pml_x(i)) {
                        real d = diff_xm(Hzp, id, sI, inv_dx);
                        real corr = d / kx_[i];
                        real& p = psi_Ey_x_[id]; p = bEx_[i] * p + cEx_[i] * d; corr += p;
                        const real be = bEy_mg_ ? (*bEy_mg_)[id] : (dt_ / eps0_);
                        Ey[id] += be * (-(corr - d));
                    }

                    // Ez: dHy/dx, dHx/dy
                    if (in_pml_x(i)) {
                        real d = diff_xm(Hyp, id, sI, inv_dx);
                        real corr = d / kx_[i];
                        real& p = psi_Ez_x_[id]; p = bEx_[i] * p + cEx_[i] * d; corr += p;
                        const real be = bEz_mg_ ? (*bEz_mg_)[id] : (dt_ / eps0_);
                        Ez[id] += be * (corr - d);
                    }
                    if (in_pml_y(j)) {
                        real d = diff_ym(Hxp, id, sJ, inv_dy);
                        real corr = d / ky_[j];
                        real& p = psi_Ez_y_[id]; p = bEy_[j] * p + cEy_[j] * d; corr += p;
                        const real be = bEz_mg_ ? (*bEz_mg_)[id] : (dt_ / eps0_);
                        Ez[id] += be * (-(corr - d));
                    }
                }
    }

private:
    // MODIFY: build_profiles_and_coeffs uses local spacing
    void build_profiles_and_coeffs(const CPMLConfig& cfg) {
        // σ_max calculation now needs to account for varying spacing
        // Option 1: Use average spacing in PML region
        auto avg_pml_spacing = [&](const std::vector<real>& spacing, int N) {
            real sum = 0;
            int count = 0;
            for (int i = 0; i < npml_; ++i) {
                sum += spacing[i] + spacing[N - 1 - i];
                count += 2;
            }
            return sum / count;
        };
        
        real dx_avg = avg_pml_spacing(dx_array_, NxT_);
        real dy_avg = avg_pml_spacing(dy_array_, NyT_);
        real dz_avg = avg_pml_spacing(dz_array_, NzT_);
        
        auto sigma_max = [&](real ds)->real {
            real d = npml_ * ds;
            return ((cfg.m + 1.0) * std::log(1.0 / cfg.Rerr) * eps0_ * c0_) / (2.0 * d);
        };
        
        real sxmax = sigma_max(dx_avg);
        real symax = sigma_max(dy_avg);
        real szmax = sigma_max(dz_avg);


        auto fill_1d = [&](int N, int npml, real smax, const std::vector<real>& spacing_arr,
            std::vector<real>& s, std::vector<real>& k, std::vector<real>& a) {
                for (int n = 0; n < N; ++n) {
                    int d = std::min(n, N - 1 - n);
                    if (d >= npml) { s[n] = 0.0; k[n] = 1.0; a[n] = 0.0; continue; }

                    real r = real(npml - (d + 0.5)) / real(npml);  // (0,1]
                    real rm = std::pow(r, cfg.m);
                    real kap = 1.0 + (cfg.kappa_max - 1.0) * rm;     // κ(n)
                    real sig = smax * kap * rm;

                    // CHANGED: Use local spacing for alpha
                    real ds_local = spacing_arr[n];
                    real alf = cfg.alpha_linear ? (cfg.alpha0 * (c0_ / ds_local) * (1.0 - r))
                                                : (cfg.alpha0 * (c0_ / ds_local) * std::pow(1.0 - r, cfg.m));
                    s[n] = sig; k[n] = kap; a[n] = alf;
                }
            };

        fill_1d((int)NxT_, npml_, sxmax, dx_array_, sx_, kx_, ax_);
        fill_1d((int)NyT_, npml_, symax, dy_array_, sy_, ky_, ay_);
        fill_1d((int)NzT_, npml_, szmax, dz_array_, sz_, kz_, az_);

        auto build_bc = [&](    const std::vector<real>& s,
                                const std::vector<real>& k,
                                const std::vector<real>& a,
                                real norm,                      // ★ Electric side pass ε0_, magnetic side pass μ0_
                                std::vector<real>& b,
                                std::vector<real>& c    ) {
                int N = (int)s.size();
                for (int n = 0; n < N; ++n) {

                    // Normalized σ (unit 1/s)
                    real sig_n = s[n] / norm;
                    real kap = k[n];
                    real alf = a[n];

                    real be = std::exp(-(sig_n / kap + alf) * dt_);    // Recursive coefficient
                    be = std::max<real>(be, 1e-30);

                    real ce = 0.0;
                    if (std::abs(sig_n) > 1e-30 || std::abs(alf) > 1e-30) {
                        ce = (sig_n * (be - 1.0)) / (kap * (sig_n + kap * alf));
                    }
                    b[n] = be; c[n] = ce;
                }
            };

        build_bc(sx_, kx_, ax_, eps0_, bEx_, cEx_);
        build_bc(sy_, ky_, ay_, eps0_, bEy_, cEy_);
        build_bc(sz_, kz_, az_, eps0_, bEz_, cEz_);

        build_bc(sx_, kx_, ax_, mu0_, bHx_, cHx_);
        build_bc(sy_, ky_, ay_, mu0_, bHy_, cHy_);
        build_bc(sz_, kz_, az_, mu0_, bHz_, cHz_);

    }
};

// === Factory ===
inline std::unique_ptr<IBoundary> make_boundary(
    const BoundaryParams& params,
    size_t NxCore, size_t NyCore, size_t NzCore,
    const GridSpacing& grid_spacing, 
    real dt, real eps0, real mu0, real c0)
{
    switch (params.type) {
    case BcType::PEC:
        return std::make_unique<PECBoundary>(NxCore, NyCore, NzCore);
    case BcType::CPML_RC:
    default:
        return std::make_unique<CPMLBoundary>(NxCore, NyCore, NzCore, grid_spacing, dt, eps0, mu0, c0, params.cpml);
    }
}
