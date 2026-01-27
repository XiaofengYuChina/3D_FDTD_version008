// fdtd_stepper.hpp - Core FDTD time-stepping functions

#pragma once
#include <cstddef>
#include "omp_config.hpp"
#include "global_function.hpp"

template<typename Real>
inline void fdtd_update_H(
    size_t NxT, size_t NyT, size_t NzT,
    const Real* __restrict dx_array,
    const Real* __restrict dy_array,
    const Real* __restrict dz_array,
    const Real* __restrict inv_dx_array,
    const Real* __restrict inv_dy_array,
    const Real* __restrict inv_dz_array,
    const Real* __restrict aHx, const Real* __restrict bHx,
    const Real* __restrict aHy, const Real* __restrict bHy,
    const Real* __restrict aHz, const Real* __restrict bHz,
    Real* __restrict Ex, Real* __restrict Ey, Real* __restrict Ez,
    Real* __restrict Hx, Real* __restrict Hy, Real* __restrict Hz)
{
    using namespace fdtd_math;
    const size_t sI = NyT * NzT;
    const size_t sJ = NzT;
    const size_t sK = 1;

#if FDTD_OMP_ENABLED
#pragma omp parallel for
#endif
    for (int i = 0; i < (int)(NxT - 1); ++i)
        for (int j = 0; j < (int)(NyT - 1); ++j)
            for (int k = 0; k < (int)(NzT - 1); ++k) {
                size_t id = idx3((size_t)i, (size_t)j, (size_t)k, NyT, NzT);
                const Real inv_dx = inv_dx_array[i];
                const Real inv_dy = inv_dy_array[j];
                const Real inv_dz = inv_dz_array[k];

                Real curlEx = diff_y(Ez, id, sJ, inv_dy) - diff_z(Ey, id, sK, inv_dz);
                Real curlEy = diff_z(Ex, id, sK, inv_dz) - diff_x(Ez, id, sI, inv_dx);
                Real curlEz = diff_x(Ey, id, sI, inv_dx) - diff_y(Ex, id, sJ, inv_dy);

                Hx[id] = aHx[id] * Hx[id] - bHx[id] * curlEx;
                Hy[id] = aHy[id] * Hy[id] - bHy[id] * curlEy;
                Hz[id] = aHz[id] * Hz[id] - bHz[id] * curlEz;
            }
}

template<typename Real>
inline void fdtd_update_E(
    size_t NxT, size_t NyT, size_t NzT,
    const Real* __restrict dx_array,
    const Real* __restrict dy_array,
    const Real* __restrict dz_array,
    const Real* __restrict inv_dx_array,
    const Real* __restrict inv_dy_array,
    const Real* __restrict inv_dz_array,
    const Real* __restrict aEx, const Real* __restrict bEx,
    const Real* __restrict aEy, const Real* __restrict bEy,
    const Real* __restrict aEz, const Real* __restrict bEz,
    Real* __restrict Ex, Real* __restrict Ey, Real* __restrict Ez,
    Real* __restrict Hx, Real* __restrict Hy, Real* __restrict Hz,
    const Real* __restrict Jx, const Real* __restrict Jy, const Real* __restrict Jz)
{
    using namespace fdtd_math;
    const size_t sI = NyT * NzT;
    const size_t sJ = NzT;
    const size_t sK = 1;

#if FDTD_OMP_ENABLED
#pragma omp parallel for
#endif
    for (int i = 1; i < (int)NxT; ++i)
        for (int j = 1; j < (int)NyT; ++j)
            for (int k = 1; k < (int)NzT; ++k) {
                size_t id = idx3((size_t)i, (size_t)j, (size_t)k, NyT, NzT);
                const Real inv_dx = inv_dx_array[i];
                const Real inv_dy = inv_dy_array[j];
                const Real inv_dz = inv_dz_array[k];

                Real curlHx = diff_ym(Hz, id, sJ, inv_dy) - diff_zm(Hy, id, sK, inv_dz);
                Real curlHy = diff_zm(Hx, id, sK, inv_dz) - diff_xm(Hz, id, sI, inv_dx);
                Real curlHz = diff_xm(Hy, id, sI, inv_dx) - diff_ym(Hx, id, sJ, inv_dy);

                Ex[id] = aEx[id] * Ex[id] + bEx[id] * (curlHx - Jx[id]);
                Ey[id] = aEy[id] * Ey[id] + bEy[id] * (curlHy - Jy[id]);
                Ez[id] = aEz[id] * Ez[id] + bEz[id] * (curlHz - Jz[id]);
            }
}
