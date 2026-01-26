// user_config.hpp — All user-configurable simulation parameters
//
// This is the ONLY file you need to modify to configure your simulation.
// All parameters are centralized here for easy access.
//
// Usage: Modify the values below, then rebuild and run the simulation.

#pragma once

#include <vector>
#include <string>
#include <limits>  // for std::numeric_limits (NaN for default physical coords)
#include "tls_materials.hpp"  // TLS material library

namespace UserConfig {

// ============================================================================
//                         1. DOMAIN AND GRID SETTINGS
// ============================================================================

// Physical domain boundaries (meters)
// This defines the simulation region using absolute physical coordinates.
// All positions (sources, detectors, structures) use the same coordinate system.
// The domain excludes PML boundaries which are added automatically.
constexpr double DOMAIN_X_MIN = 0.0;        // Physical domain x minimum (m)
constexpr double DOMAIN_X_MAX = 4000e-9;    // Physical domain x maximum (m) = 2 um
constexpr double DOMAIN_Y_MIN = 0.0;        // Physical domain y minimum (m)
constexpr double DOMAIN_Y_MAX = 4000e-9;    // Physical domain y maximum (m) = 2 um
constexpr double DOMAIN_Z_MIN = 0.0;        // Physical domain z minimum (m)
constexpr double DOMAIN_Z_MAX = 1000e-9;    // Physical domain z maximum (m) = 1 um

// Derived domain sizes (for backward compatibility in some calculations)
constexpr double DOMAIN_SIZE_X = DOMAIN_X_MAX - DOMAIN_X_MIN;
constexpr double DOMAIN_SIZE_Y = DOMAIN_Y_MAX - DOMAIN_Y_MIN;
constexpr double DOMAIN_SIZE_Z = DOMAIN_Z_MAX - DOMAIN_Z_MIN;

// Central wavelength (meters) - used for mesh generation and source
constexpr double LAMBDA0 = 1500e-9;          // 300 nm

// Mesh mode selection:
// 0 = UNIFORM           - constant grid spacing
// 1 = AUTO_NONUNIFORM   - ppw-based auto mesh with interface snapping (recommended)
constexpr int MESH_MODE = 1;

// --- Uniform mesh settings (for MESH_MODE = 0) ---
constexpr double DX_UNIFORM = 10e-9;        // 10 nm uniform spacing
constexpr double DY_UNIFORM = 10e-9;
constexpr double DZ_UNIFORM = 10e-9;

// --- Auto mesh settings (for MESH_MODE = 1) ---
// Mesh accuracy level determines points-per-wavelength (PPW):
//   acc=1 -> 6 ppw  (fast, less accurate)
//   acc=2 -> 10 ppw (balanced)
//   acc=3 -> 14 ppw (high accuracy)
//   acc>=3: 14 + (acc-3)*4 ppw (e.g., acc=8 -> 34 ppw)
constexpr int    AUTO_MESH_ACCURACY = 2;               // Mesh accuracy level (1-8)
constexpr double AUTO_MESH_PPW_OVERRIDE = 20.0;        // Direct PPW override (if > 0, overrides AUTO_MESH_ACCURACY)
constexpr double AUTO_MESH_DX_MIN = 2e-9;              // Minimum spacing (m)
constexpr double AUTO_MESH_DX_MAX = 100e-9;            // Maximum spacing (m)
constexpr double AUTO_MESH_MAX_GRADING_RATIO = 1.41421356237;  // sqrt(2) by default

// Grid cell counts (used for UNIFORM mode only)
// For AUTO_NONUNIFORM, these are computed from domain size
constexpr size_t NX = 100;
constexpr size_t NY = 100;
constexpr size_t NZ = 100;


// ============================================================================
//                         2. TIME STEPPING SETTINGS
// ============================================================================

// CFL safety factor (must be < 1 for stability)
constexpr double CFL_FACTOR = 0.99;

// Total number of time steps
constexpr size_t N_STEPS = 10000;

// Save data every N steps
constexpr size_t SAVE_EVERY = 10;


// ============================================================================
//                         3. BOUNDARY CONDITIONS (CPML)
// ============================================================================

// Boundary type: 0 = PEC (perfect electric conductor), 1 = CPML (absorbing)
constexpr int BOUNDARY_TYPE = 1;            // 1 = CPML (recommended)

// CPML parameters
constexpr int    CPML_NPML = 8;            // PML thickness in cells
constexpr double CPML_M = 3.0;              // Polynomial grading order
constexpr double CPML_RERR = 1e-8;          // Target reflection coefficient
constexpr double CPML_ALPHA0 = 0.2;         // Alpha parameter
constexpr double CPML_KAPPA_MAX = 10.0;     // Maximum kappa
constexpr bool   CPML_ALPHA_LINEAR = true;  // Linear alpha distribution


// ============================================================================
//                         4. SOURCE CONFIGURATION
// ============================================================================
//
// Source Types Available:
//   0 = Dipole (point current source)
//   1 = PlaneWave (plane wave source)
//
// Waveform Types:
//   0 = Ricker (Mexican hat wavelet)
//   1 = GaussianModulatedSine
//   2 = RickerLikeGaussian2nd
//   3 = ContinuousWave (with smooth turn-on)
//
// Polarization:
//   0 = Ex, 1 = Ey, 2 = Ez, 3 = Hx, 4 = Hy, 5 = Hz

// -------------------- Dipole Source Definition --------------------
struct DipoleSourceDef {
    bool enabled;               // Whether this source is active
    double x, y, z;             // Position (meters)
    double amplitude;           // Peak current (A)
    double frequency;           // Frequency (Hz), 0 = use wavelength
    double wavelength;          // Wavelength (m), used if frequency <= 0
    double tau;                 // Time constant (s), 0 = auto
    double df_fwhm;             // Bandwidth FWHM (Hz), used if tau <= 0
    double t0_factor;           // t0 = t0_factor * tau_eff
    int waveform;               // 0=Ricker, 1=GaussianModulatedSine, 2=RickerLikeGaussian2nd, 3=CW
    int polarization;           // 0=Ex, 1=Ey, 2=Ez
};

// Define dipole sources here
inline const std::vector<DipoleSourceDef> DIPOLE_SOURCES = {
    // Default dipole source
    {true, 1100e-9, 2000e-9, 500e-9,   // enabled, x, y, z
     1e-5,                              // amplitude (A)
     0.0, 1500e-9,                      // frequency (Hz), wavelength (m)
     5e-15, 2.4e13,                     // tau (s), df_fwhm (Hz)
     3.0,                               // t0_factor
     0,                                 // waveform (Ricker)
     2},                                // polarization (Ez)
};

// -------------------- Plane Wave Source Definition --------------------
struct PlaneWaveSourceDef {
    bool enabled;               // Whether this source is active
    double injection_position;  // Position along propagation axis (meters)
    double amplitude;           // Peak E-field (V/m)
    double frequency;           // Frequency (Hz), 0 = use wavelength
    double wavelength;          // Wavelength (m)
    double tau;                 // Time constant (s), 0 = auto
    double df_fwhm;             // Bandwidth FWHM (Hz)
    double t0_factor;           // t0 = t0_factor * tau_eff
    int waveform;               // 0=Ricker, 1=GaussianModulatedSine, 2=RickerLikeGaussian2nd, 3=CW
    int direction;              // 0=+X, 1=-X, 2=+Y, 3=-Y, 4=+Z, 5=-Z
    int polarization;           // 0=Ex, 1=Ey, 2=Ez (perpendicular to propagation)
    // Injection region bounds (0,0,0,0,0,0 = full domain)
    double x_min, x_max, y_min, y_max, z_min, z_max;
};

// Define plane wave sources here
inline const std::vector<PlaneWaveSourceDef> PLANE_WAVE_SOURCES = {
    // Example: disabled plane wave
    // {true, 200e-9,                    // enabled, injection_position
    //  1.0,                             // amplitude (V/m)
    //  0.0, 500e-9,                     // frequency (Hz), wavelength (m)
    //  0.0, 0.0,                        // tau (s), df_fwhm (Hz) - auto
    //  3.0,                             // t0_factor
    //  1,                               // waveform (GaussianModulatedSine)
    //  4,                               // direction (+Z)
    //  0,                               // polarization (Ex)
    //  0, 0, 0, 0, 0, 0},               // injection region (full domain)
};

// ============================================================================
//                         5. DETECTOR CONFIGURATION
// ============================================================================
//
// Detector Types Available:
//   - MeshDetector: Outputs refractive index distribution and mesh grid info
//   - FieldMovie2D: Records 2D slices of field components over time
//   - PointFieldDetector: Records field time series at a single point
//
// Field Components:
//   0=Ex, 1=Ey, 2=Ez, 3=Hx, 4=Hy, 5=Hz,
//   6=E_magnitude, 7=H_magnitude, 8=Sx, 9=Sy, 10=Sz, 11=S_magnitude
//
// Slice Planes:
//   0=XY (z=const), 1=XZ (y=const), 2=YZ (x=const)

// Output directory tag (results saved to frames/<RUN_TAG>/)
inline const std::string RUN_TAG = "3D_FDTD_v008_output";

// -------------------- Mesh Detector Definition --------------------
struct MeshDetectorDef {
    bool enabled;               // Whether this detector is active
    std::string name;           // Output directory name
    double slice_position;      // Position along normal axis (meters)
    int slice_plane;            // 0=XY, 1=XZ, 2=YZ
    bool export_mesh_lines;     // Export cell boundary coordinates
    bool export_spacing;        // Export dx, dy, dz arrays
};

// Define mesh detectors here
inline const std::vector<MeshDetectorDef> MESH_DETECTORS = {
    {true, "mesh_info", 500e-9, 0, true, true},  // XY plane at z=500nm
};

// -------------------- Field Movie 2D Definition --------------------
struct FieldMovie2DDef {
    bool enabled;               // Whether this detector is active
    std::string name;           // Output directory name
    int field_component;        // 0=Ex, 1=Ey, 2=Ez, 3=Hx, 4=Hy, 5=Hz, etc.
    int slice_plane;            // 0=XY, 1=XZ, 2=YZ
    double slice_position;      // Position in physical coords (meters)
    std::string frame_pattern;  // Output filename pattern
};

// Define field movie detectors here
inline const std::vector<FieldMovie2DDef> FIELD_MOVIE_2D_DETECTORS = {
    // Ez field on XY plane at z=500nm
    {true, "Ez_movie", 2, 0, 500e-9, "ez_%04d.raw"},
    // Uncomment to add more:
    // {true, "Ex_movie", 0, 0, 500e-9, "ex_%04d.raw"},
    // {true, "Hy_movie", 4, 0, 500e-9, "hy_%04d.raw"},
};

// -------------------- Point Field Detector Definition --------------------
struct PointFieldDetectorDef {
    bool enabled;               // Whether this detector is active
    std::string name;           // Output directory name
    double x, y, z;             // Position (meters)
    std::vector<int> components; // Field components to record (0=Ex, 1=Ey, 2=Ez, etc.)
};

// Define point field detectors here
inline const std::vector<PointFieldDetectorDef> POINT_FIELD_DETECTORS = {
    // Ez probe at center
    {true, "Ez_probe", 2000e-9, 2000e-9, 500e-9, {2}},  // Ez only
    // Uncomment to add more:
    // {true, "E_probe", 1500e-9, 1500e-9, 500e-9, {0, 1, 2}},  // Ex, Ey, Ez
};

// -------------------- Legacy compatibility settings --------------------
// 2D slice detector: z-slice position in physical coordinates (meters)
constexpr double Z_SLICE_Z = 500e-9;    // Z position of the XY slice (m)

// Frame output pattern and precision
inline const std::string FRAME_PATTERN = "ez_%04d.raw";
constexpr bool WRITE_FLOAT64 = true;        // true = double, false = float

// 1D probe position in physical coordinates (meters)
// Must specify explicit positions within the domain
constexpr double PROBE_X = 2000e-9;     // Probe x position (m) - domain center
constexpr double PROBE_Y = 2000e-9;     // Probe y position (m)
constexpr double PROBE_Z = 500e-9;      // Probe z position (m)


// ============================================================================
//                         6. TWO-LEVEL SYSTEM (GAIN MEDIUM) CONFIGURATION
// ============================================================================
//
// TLS materials are defined in: tls_materials.hpp
// Reference them by name in STRUCTURES using the tls_material field.
//
// Example:
//   In tls_materials.hpp: {"gain_1500nm", 1500e-9, 7e12, 1e-12, 1e25, 1.0}
//   In STRUCTURES below:  {"box", 3.0, {...}, "gain_1500nm"}  // with TLS
//                         {"box", 4.0, {...}, ""}             // no TLS

// Global enable for two-level system (master switch)
// If false, ALL TLS will be disabled regardless of per-structure settings
constexpr bool TLS_ENABLED = true;

// Enable population clamping (optional, for debugging)
constexpr bool TLS_ENABLE_CLAMP = false;


// ============================================================================
//                         7. STABILITY MONITOR THRESHOLDS
// ============================================================================

// Stability monitor interval (steps)
// How often to check stability and output physics monitoring to console
constexpr size_t STABILITY_MONITOR_INTERVAL = 10;

// Maximum E-field magnitude before declaring instability (V/m)
constexpr double STABILITY_MAX_E = 1e10;

// Maximum polarization magnitude before declaring instability (C/m^2)
constexpr double STABILITY_MAX_P = 1e5;

// Maximum conservation error (fractional)
constexpr double STABILITY_CONSERVATION_ERR = 0.01;  // 1%

// Maximum allowed growth rate per 100 steps
constexpr double STABILITY_GROWTH_RATE = 10.0;


// ============================================================================
//                         8. STRUCTURE DEFINITIONS
// ============================================================================
//
// Define your simulation structures here.
// Each structure is specified by its type, position, dimensions, and material.
//
// Structure types:
//   "box"      - rectangular box
//   "sphere"   - sphere
//   "cylinder" - cylinder along z-axis
//
// Material is specified by refractive index n (assuming non-magnetic, mu_r = 1)

struct StructureDef {
    std::string type;           // "box", "sphere", or "cylinder"
    double n;                   // Refractive index

    // For box: x0, x1, y0, y1, z0, z1 (bounds in meters, relative to core origin)
    // For sphere: cx, cy, cz, radius (center and radius in meters)
    // For cylinder: cx, cy, radius, z0, z1 (center, radius, and z-bounds)
    double params[6];

    // ===== Two-Level System (TLS) binding =====
    // Set to a TLS material name defined in TLS_MATERIALS to add gain medium
    // Leave empty ("") for pure dielectric (no TLS)
    std::string tls_material = "";
};

// Define structures here
// Coordinates are relative to the physical domain origin (after PML)
// Use negative values for "relative to center" positioning
//
// Example configurations:

// Structure definition format:
//   {type, n, {params...}, tls_material}
//
// tls_material: Name of TLS material from TLS_MATERIALS, or "" for no TLS
//
// Examples:
//   {"box", 2.0, {...}, ""}              // Pure dielectric, no TLS
//   {"box", 3.0, {...}, "gain_1500nm"}   // With TLS using "gain_1500nm" material
//
inline const std::vector<StructureDef> STRUCTURES = {
    // Cylinder with TLS (gain medium) - references "gain_1500nm" from TLS_MATERIALS
    // params = {cx, cy, radius, z0, z1, 0}
    {"cylinder", 3.3, {2000e-9, 2000e-9, 1000e-9, 400e-9, 600e-9, 0}, "gain_1500nm"},

    // Example: Box WITHOUT TLS (pure dielectric) - empty string means no TLS
    // {"box", 4.0, {900e-9, 1100e-9, 900e-9, 1100e-9, 400e-9, 600e-9}, ""},

    // Example: Multiple structures with different TLS materials
    // {"box", 3.0, {100e-9, 300e-9, 100e-9, 300e-9, 400e-9, 600e-9}, "gain_1300nm"},
};

// Background material refractive index (default: air, n=1)
constexpr double BACKGROUND_N = 1.0;


// ============================================================================
//                         10. MESH OVERRIDE REGIONS (HIGHEST PRIORITY)
// ============================================================================
//
// Define custom mesh refinement regions here.
// These have the HIGHEST priority and will override auto mesh settings.
// Useful for:
// - Accurately resolving curved boundaries (like cylinders, spheres)
// - Fine-tuning mesh in specific regions of interest
// - Creating uniform mesh in certain areas
//
// Set MESH_OVERRIDE_ENABLED to true to activate this feature.

constexpr bool MESH_OVERRIDE_ENABLED = true;

// Mesh Override Region definition
struct MeshOverrideRegion {
    bool enabled;               // Whether this region is active
    double x0, x1;              // X bounds (meters, relative to physical domain origin)
    double y0, y1;              // Y bounds (meters, relative to physical domain origin)
    double z0, z1;              // Z bounds (meters, relative to physical domain origin)
    double dx;                  // Target mesh size in X direction (meters)
    double dy;                  // Target mesh size in Y direction (meters)
    double dz;                  // Target mesh size in Z direction (meters)
};

// Define mesh override regions here
// Coordinates are relative to the physical domain origin (after PML)
//
// Example: Refine mesh around the cylinder boundary for accurate curved surface
inline const std::vector<MeshOverrideRegion> MESH_OVERRIDE_REGIONS = {
    // Region 1: Fine mesh around cylinder (r=100nm at x=700nm, y=700nm)
    // Covers the cylinder plus a margin, with 5nm uniform mesh
    // {true,   600e-9, 800e-9,   600e-9, 800e-9,   400e-9, 600e-9,   5e-9, 5e-9, 5e-9},
    // Region 2: Fine mesh around first box (n=4)
    {true,   1000e-9, 3000e-9,   1000e-9, 3000e-9,   400e-9, 600e-9,   20e-9, 20e-9, 20e-9},
    // Region 3: Fine mesh around third box (n=3)
    // {true,   1200e-9, 1400e-9,   1200e-9, 1400e-9,   400e-9, 600e-9,   5e-9, 5e-9, 5e-9},
};

// NOTE: Override regions can overlap. When they do, the SMALLEST dx/dy/dz wins.
// This allows you to define a moderately fine region, then an even finer sub-region.


// ============================================================================
//                         9. PHYSICAL CONSTANTS
// ============================================================================
// These are fundamental constants - typically you don't need to modify them

namespace PhysicalConstants {
    constexpr double PI   = 3.14159265358979323846;
    constexpr double C0   = 299792458.0;              // Speed of light (m/s)
    constexpr double EPS0 = 8.854187817e-12;          // Vacuum permittivity (F/m)
    constexpr double MU0  = 1.2566370614359173e-6;    // Vacuum permeability (H/m)
    constexpr double HBAR = 1.054571817e-34;          // Reduced Planck constant (J·s)
}


// ============================================================================
//                     11. PARALLEL DOMAIN DECOMPOSITION V2
// ============================================================================
//
// TRUE domain decomposition for improved performance on large grids.
// The simulation domain is split into multiple subdomains along a chosen axis.
//
// V2 Architecture (efficient):
// ┌──────────────────┐  ┌──────────────────┐  ┌──────────────────┐
// │ PML │ 物理 │Halo│  │Halo│ 物理 │Halo│  │Halo│ 物理 │ PML │
// │     │ 区域 │    │  │    │ 区域 │    │  │    │ 区域 │     │
// │     域0         │  │       域1       │  │         域2     │
// └──────────────────┘  └──────────────────┘  └──────────────────┘
//                   ↕ O(N²) halo exchange ↕
//
// Key features of V2:
// - Each subdomain has its OWN independent PML data (psi arrays)
// - Only halo regions exchanged per step: O(N²) instead of O(N³)
// - NO scatter/gather operations during computation
// - Subdomains compute completely independently
//
// When to use parallel mode:
// - Large grids where domain decomposition improves cache efficiency
// - Preparing for future MPI distributed computing
// - Testing domain decomposition correctness
//
// When to use serial mode:
// - Small to medium grids (OpenMP loop parallelization is sufficient)
// - When TLS (two-level system) is heavily used (currently requires global ops)
//
// Requirements:
// - OpenMP must be enabled during compilation (-fopenmp flag)
// - Grid must be large enough: at least (4 + 2*halo_width) cells per domain

// Enable parallel domain decomposition V2
constexpr bool PARALLEL_ENABLED = false;

// Number of parallel domains
// Recommended: 2-8 domains, ideally matching CPU core count
constexpr int PARALLEL_NUM_DOMAINS = 8;

// Decomposition axis:
// 0 = X axis, 1 = Y axis, 2 = Z axis
// Choose the axis with the largest number of cells for best load balancing
constexpr int PARALLEL_DECOMP_AXIS = 0;

// Halo (ghost cell) width for inter-domain communication
// 1 is sufficient for standard FDTD (nearest-neighbor stencil)
constexpr int PARALLEL_HALO_WIDTH = 1;

} // namespace UserConfig
