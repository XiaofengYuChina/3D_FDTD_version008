cmake_minimum_required(VERSION 3.20)
project(3D_FDTD_version008 LANGUAGES CXX)

# C++20
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# === Optional switch: whether to enable OpenMP parallelization ===
option(ENABLE_OPENMP "Enable OpenMP parallelization" ON)

# Executable
add_executable(3D_FDTD_v008 src/3D_FDTD_version008.cpp)

# General compilation optimization options
if(MSVC)
    target_compile_options(3D_FDTD_v008 PRIVATE /DNOMINMAX /permissive-)
    target_compile_options(3D_FDTD_v008 PRIVATE
        "$<$<CONFIG:Release>:/O2>"
        "$<$<CONFIG:Debug>:/RTC1>"
    )
else()
    target_compile_options(3D_FDTD_v008 PRIVATE -march=native)
    target_compile_options(3D_FDTD_v008 PRIVATE
        "$<$<CONFIG:Release>:-O3>"
    )
endif()

# Enable OpenMP as needed (after enabling, compiler will automatically define _OPENMP)
if (ENABLE_OPENMP)
  find_package(OpenMP)
  if (OpenMP_CXX_FOUND)
    target_link_libraries(3D_FDTD_v008 PRIVATE OpenMP::OpenMP_CXX)
  else()
    message(WARNING "ENABLE_OPENMP=ON but OpenMP not found; building WITHOUT OpenMP")
  endif()
endif()

# Test executable
add_executable(test_auto_mesh tests/test_auto_mesh.cpp)
target_include_directories(test_auto_mesh PRIVATE ${CMAKE_SOURCE_DIR}/src)

if(MSVC)
    target_compile_options(test_auto_mesh PRIVATE /DNOMINMAX)
else()
    target_compile_options(test_auto_mesh PRIVATE -march=native)
endif()

# Mesh override verification test
add_executable(verify_mesh_override tests/verify_mesh_override.cpp)
target_include_directories(verify_mesh_override PRIVATE ${CMAKE_SOURCE_DIR}/src)

if(MSVC)
    target_compile_options(verify_mesh_override PRIVATE /DNOMINMAX)
else()
    target_compile_options(verify_mesh_override PRIVATE -march=native)
endif()

