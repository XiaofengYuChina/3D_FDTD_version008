// omp_config.hpp —— Parallel Switch

#pragma once

#if defined(_OPENMP)
#include <omp.h>
#define FDTD_OMP_ENABLED 1
#else
#define FDTD_OMP_ENABLED 0

//Provide a minimal 'stand-in' so the code can compile even without using #ifdef
inline int omp_get_max_threads() { return 1; }
#endif

// MSVC doesn't support OpenMP collapse clause - only enable on GCC/Clang
#if defined(_MSC_VER)
    #define OMP_PARALLEL_FOR _Pragma("omp parallel for")
#else
    #define OMP_PARALLEL_FOR _Pragma("omp parallel for collapse(3)")
#endif
