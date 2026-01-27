// omp_config.hpp - OpenMP configuration and fallback definitions

#pragma once

#if defined(_OPENMP)
#include <omp.h>
#define FDTD_OMP_ENABLED 1
#else
#define FDTD_OMP_ENABLED 0
inline int omp_get_max_threads() { return 1; }
#endif

// MSVC doesn't support OpenMP collapse clause
#if defined(_MSC_VER)
    #define OMP_PARALLEL_FOR _Pragma("omp parallel for")
#else
    #define OMP_PARALLEL_FOR _Pragma("omp parallel for collapse(3)")
#endif
