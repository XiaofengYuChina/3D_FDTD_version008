// user_config.hpp — All user-configurable simulation parameters
// This is the ONLY file you need to modify to configure your simulation.

#pragma once

#include <vector>
#include <string>
#include <limits>
#include "tls_materials.hpp"

namespace UserConfig {

// Physical domain boundaries (meters) - excludes PML boundaries which are added automatically
constexpr double DOMAIN_X_MIN = 0.0;
constexpr double DOMAIN_X_MAX = 4000e-9;
constexpr double DOMAIN_Y_MIN = 0.0;
constexpr double DOMAIN_Y_MAX = 4000e-9;
constexpr double DOMAIN_Z_MIN = 0.0;
constexpr double DOMAIN_Z_MAX = 1000e-9;

constexpr double DOMAIN_SIZE_X = DOMAIN_X_MAX - DOMAIN_X_MIN;
constexpr double DOMAIN_SIZE_Y = DOMAIN_Y_MAX - DOMAIN_Y_MIN;
constexpr double DOMAIN_SIZE_Z = DOMAIN_Z_MAX - DOMAIN_Z_MIN;

constexpr double LAMBDA0 = 1500e-9;  // Central wavelength (m)

// Mesh mode: 0 = UNIFORM, 1 = AUTO_NONUNIFORM (ppw-based with interface snapping)
constexpr int MESH_MODE = 1;

// Uniform mesh settings (MESH_MODE = 0)
constexpr double DX_UNIFORM = 10e-9;
constexpr double DY_UNIFORM = 10e-9;
constexpr double DZ_UNIFORM = 10e-9;

// Auto mesh settings (MESH_MODE = 1)
// Accuracy level: 1->6ppw, 2->10ppw, 3->14ppw, >=3: 14+(acc-3)*4 ppw
constexpr int    AUTO_MESH_ACCURACY = 2;
constexpr double AUTO_MESH_PPW_OVERRIDE = 20.0;  // Direct PPW (if > 0, overrides accuracy)
constexpr double AUTO_MESH_DX_MIN = 2e-9;
constexpr double AUTO_MESH_DX_MAX = 100e-9;
constexpr double AUTO_MESH_MAX_GRADING_RATIO = 1.41421356237;  // sqrt(2)

// Grid cell counts (UNIFORM mode only)
constexpr size_t NX = 100;
constexpr size_t NY = 100;
constexpr size_t NZ = 100;

// Time stepping
constexpr double CFL_FACTOR = 0.99;  // Must be < 1 for stability
constexpr size_t N_STEPS = 10000;
constexpr size_t SAVE_EVERY = 10;

// Boundary conditions: 0 = PEC, 1 = CPML (absorbing)
constexpr int BOUNDARY_TYPE = 1;

// CPML parameters
constexpr int    CPML_NPML = 8;
constexpr double CPML_M = 3.0;
constexpr double CPML_RERR = 1e-8;
constexpr double CPML_ALPHA0 = 0.2;
constexpr double CPML_KAPPA_MAX = 10.0;
constexpr bool   CPML_ALPHA_LINEAR = true;

// Source types: 0=Dipole, 1=PlaneWave
// Waveform: 0=Ricker, 1=GaussianModulatedSine, 2=RickerLikeGaussian2nd, 3=ContinuousWave
// Polarization: 0=Ex, 1=Ey, 2=Ez, 3=Hx, 4=Hy, 5=Hz

struct DipoleSourceDef {
    bool enabled;
    double x, y, z;         // Position (m)
    double amplitude;       // Peak current (A)
    double frequency;       // Hz (0 = use wavelength)
    double wavelength;      // m (used if frequency <= 0)
    double tau;             // Time constant (s), 0 = auto
    double df_fwhm;         // Bandwidth FWHM (Hz), used if tau <= 0
    double t0_factor;       // t0 = t0_factor * tau_eff
    int waveform;
    int polarization;
};

inline const std::vector<DipoleSourceDef> DIPOLE_SOURCES = {
    {true, 1100e-9, 2000e-9, 500e-9, 1e-5, 0.0, 1500e-9, 5e-15, 2.4e13, 3.0, 0, 2},
};

struct PlaneWaveSourceDef {
    bool enabled;
    double injection_position;  // Position along propagation axis (m)
    double amplitude;           // Peak E-field (V/m)
    double frequency;           // Hz (0 = use wavelength)
    double wavelength;          // m
    double tau;                 // Time constant (s), 0 = auto
    double df_fwhm;             // Bandwidth FWHM (Hz)
    double t0_factor;
    int waveform;
    int direction;              // 0=+X, 1=-X, 2=+Y, 3=-Y, 4=+Z, 5=-Z
    int polarization;           // Perpendicular to propagation
    double x_min, x_max, y_min, y_max, z_min, z_max;  // Injection region (0s = full domain)
};

inline const std::vector<PlaneWaveSourceDef> PLANE_WAVE_SOURCES = {};

// Detector field components: 0=Ex, 1=Ey, 2=Ez, 3=Hx, 4=Hy, 5=Hz,
//   6=E_magnitude, 7=H_magnitude, 8=Sx, 9=Sy, 10=Sz, 11=S_magnitude
// Slice planes: 0=XY (z=const), 1=XZ (y=const), 2=YZ (x=const)

inline const std::string RUN_TAG = "3D_FDTD_v008_output";

struct MeshDetectorDef {
    bool enabled;
    std::string name;
    double slice_position;  // Position along normal axis (m)
    int slice_plane;
    bool export_mesh_lines;
    bool export_spacing;
};

inline const std::vector<MeshDetectorDef> MESH_DETECTORS = {
    {true, "mesh_info", 500e-9, 0, true, true},
};

struct FieldMovie2DDef {
    bool enabled;
    std::string name;
    int field_component;
    int slice_plane;
    double slice_position;  // m
    std::string frame_pattern;
};

inline const std::vector<FieldMovie2DDef> FIELD_MOVIE_2D_DETECTORS = {
    {true, "Ez_movie", 2, 0, 500e-9, "ez_%04d.raw"},
};

struct PointFieldDetectorDef {
    bool enabled;
    std::string name;
    double x, y, z;  // m
    std::vector<int> components;
};

inline const std::vector<PointFieldDetectorDef> POINT_FIELD_DETECTORS = {
    {true, "Ez_probe", 2000e-9, 2000e-9, 500e-9, {2}},
};

// Legacy compatibility
constexpr double Z_SLICE_Z = 500e-9;
inline const std::string FRAME_PATTERN = "ez_%04d.raw";
constexpr bool WRITE_FLOAT64 = true;
constexpr double PROBE_X = 2000e-9;
constexpr double PROBE_Y = 2000e-9;
constexpr double PROBE_Z = 500e-9;

// Two-Level System (TLS) - materials defined in tls_materials.hpp
// Reference by name in STRUCTURES using tls_material field ("" for no TLS)
constexpr bool TLS_ENABLED = true;
constexpr bool TLS_ENABLE_CLAMP = false;

// Stability monitor
constexpr size_t STABILITY_MONITOR_INTERVAL = 10;
constexpr double STABILITY_MAX_E = 1e10;              // Max E-field (V/m)
constexpr double STABILITY_MAX_P = 1e5;               // Max polarization (C/m^2)
constexpr double STABILITY_CONSERVATION_ERR = 0.01;   // Max fractional error
constexpr double STABILITY_GROWTH_RATE = 10.0;        // Max growth per 100 steps

// Structure types: "box", "sphere", "cylinder"
struct StructureDef {
    std::string type;
    double n;  // Refractive index

    // box: x0, x1, y0, y1, z0, z1
    // sphere: cx, cy, cz, radius, -, -
    // cylinder: cx, cy, radius, z0, z1, -
    double params[6];

    std::string tls_material = "";  // TLS material name or "" for none
};

inline const std::vector<StructureDef> STRUCTURES = {
    {"cylinder", 3.3, {2000e-9, 2000e-9, 1000e-9, 400e-9, 600e-9, 0}, "gain_1500nm"},
};

constexpr double BACKGROUND_N = 1.0;

// Mesh override regions - highest priority, override auto mesh settings
constexpr bool MESH_OVERRIDE_ENABLED = true;

struct MeshOverrideRegion {
    bool enabled;
    double x0, x1, y0, y1, z0, z1;  // Bounds (m)
    double dx, dy, dz;              // Target mesh size (m)
};

// When regions overlap, the smallest dx/dy/dz wins
inline const std::vector<MeshOverrideRegion> MESH_OVERRIDE_REGIONS = {
    {true, 1000e-9, 3000e-9, 1000e-9, 3000e-9, 400e-9, 600e-9, 20e-9, 20e-9, 20e-9},
};

namespace PhysicalConstants {
    constexpr double PI   = 3.14159265358979323846;
    constexpr double C0   = 299792458.0;
    constexpr double EPS0 = 8.854187817e-12;
    constexpr double MU0  = 1.2566370614359173e-6;
    constexpr double HBAR = 1.054571817e-34;
}

// Parallel domain decomposition V2
// Splits domain into subdomains along chosen axis with halo exchange
// Use for large grids; serial mode (OpenMP loops) is better for small/medium grids
constexpr bool PARALLEL_ENABLED = false;
constexpr int PARALLEL_NUM_DOMAINS = 8;     // Recommended: match CPU core count
constexpr int PARALLEL_DECOMP_AXIS = 0;     // 0=X, 1=Y, 2=Z (choose largest axis)
constexpr int PARALLEL_HALO_WIDTH = 1;      // 1 is sufficient for standard FDTD

} // namespace UserConfig
