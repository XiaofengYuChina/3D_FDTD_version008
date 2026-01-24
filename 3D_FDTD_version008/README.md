# 3D FDTD Electromagnetic Simulator

A high-performance 3D Finite-Difference Time-Domain (FDTD) electromagnetic simulator written in C++20 with Python visualization tools. This project simulates electromagnetic wave propagation in 3D space with support for custom geometries, materials, sources, and detectors.

## Features

- **3D Yee-grid electromagnetic simulation** with CFL-stable time stepping
- **Configurable structure, materials, and sources** via header-based configuration
- **CPU-optimized field update kernels** with optional OpenMP parallelization
- **Multiple detector types**: 2D field slices, 1D point probes, refractive index profiles, Poynting vector flux, and total EM energy
- **Python visualization scripts** for generating movies, plots, and analysis
- **CPML (Convolutional Perfectly Matched Layer)** boundary conditions for absorbing boundaries

## Project Structure

```
3D_FDTD_version008/
├── src/                          # C++ source code
│   ├── 3D_FDTD_version008.cpp   # Main simulation loop
│   ├── params.hpp                # Simulation parameters and configuration
│   ├── boundary.hpp              # Boundary condition implementations
│   ├── fdtd_stepper.hpp          # FDTD field update equations
│   ├── detectors.hpp             # Detector implementations
│   ├── source.hpp                # Source implementations
│   ├── structure_material.hpp    # Material and geometry definitions
│   ├── global_function.hpp       # Utility functions
│   └── omp_config.hpp            # OpenMP configuration
├── Py_Visualization/             # Python visualization scripts
│   ├── make_Ez_movie.py          # Generate Ez field evolution movie
│   ├── make_Ez_probe.py          # Plot Ez point probe time series
│   ├── make_em_total_energy.py   # Plot total EM energy vs time
│   ├── make_box_poynting.py      # Visualize Poynting vector flux
│   └── make_index_profile.py     # Plot refractive index profile
├── CMakeLists.txt                # CMake build configuration
├── CMakePresets.json             # CMake presets for different build configurations
└── requirements.txt               # Python dependencies

Output files are written to: frames/3D_FDTD_v008_output/
```

## Requirements

### C++ Build Requirements
- **CMake 3.20+**
- **C++20 compatible compiler**:
  - Windows: Visual Studio 2022 Build Tools (or full VS) with MSVC
  - Linux/macOS: GCC 10+ or Clang 12+
- **(Optional) OpenMP** for parallelization support

### Python Requirements
- **Python 3.7+**
- Dependencies listed in `requirements.txt`:
  - numpy
  - matplotlib
  - imageio
  - ffmpeg (for MP4 generation)

## Building the Project

### Windows (Visual Studio)

1. **Configure the project:**
   ```bash
   cmake --preset x64-release
   ```

2. **Build the executable:**
   ```bash
   cmake --build out/build/x64-release --config Release
   ```

3. **The executable will be created at:**
   ```
   out/build/x64-release/Release/3D_FDTD_v008.exe
   ```

### Linux/macOS

1. **Configure the project:**
   ```bash
   cmake -S . -B out/build -DCMAKE_BUILD_TYPE=Release
   ```

2. **Build the executable:**
   ```bash
   cmake --build out/build --config Release
   ```

3. **The executable will be created at:**
   ```
   out/build/3D_FDTD_v008
   ```

### Optional: Enable OpenMP

To enable OpenMP parallelization for faster simulations:

**Windows:**
```bash
cmake --preset x64-release -DENABLE_OPENMP=ON
cmake --build out/build/x64-release --config Release
```

**Linux/macOS:**
```bash
cmake -S . -B out/build -DCMAKE_BUILD_TYPE=Release -DENABLE_OPENMP=ON
cmake --build out/build --config Release
```

## Running the Simulation

1. **Run the C++ simulation:**
   ```bash
   # Windows
   out/build/x64-release/Release/3D_FDTD_v008.exe
   
   # Linux/macOS
   ./out/build/3D_FDTD_v008
   ```

2. **The simulation will:**
   - Print progress messages showing `Ez(center)` values at regular intervals
   - Generate output files in `frames/3D_FDTD_v008_output/`:
     - `Ez_center/` - 2D slices of Ez field at z-center
     - `Ez_probe/` - 1D time series of Ez at probe point
     - `n_center/` - Refractive index profile
     - `box_flux80/` - Poynting vector flux data
     - `EM_total_energy/` - Total EM energy in box

3. **Example output:**
   ```
   OpenMP NOT enabled
   step 0 Ez(center)=...
   step 50 Ez(center)=...
   ...
   Done. dt = ... fs
   Done. t_total = ... fs
   ```

## Python Visualization

**Note:** The Python visualization scripts are **not automatically run** after the C++ simulation. You need to run them manually after the simulation completes.

### Setup Python Environment

1. **Install Python dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

2. **FFmpeg Installation (for MP4 generation):**
   
   The `imageio-ffmpeg` package (included in requirements.txt) should automatically provide ffmpeg binaries. However, if you encounter issues, you can install ffmpeg manually:
   
   **Option 1: Using imageio-ffmpeg (Recommended - Already in requirements.txt)**
   ```bash
   pip install imageio-ffmpeg
   ```
   This package bundles ffmpeg binaries and should work automatically.
   
   **Option 2: System Installation**
   
   **Windows:**
   - Download from: https://www.gyan.dev/ffmpeg/builds/ (or https://ffmpeg.org/download.html)
   - Extract and add the `bin` folder to your system PATH
   - Or use package managers:
     - **Chocolatey:** `choco install ffmpeg`
     - **Scoop:** `scoop install ffmpeg`
   
   **Linux:**
   ```bash
   sudo apt-get update
   sudo apt-get install ffmpeg
   # or for other distributions:
   sudo yum install ffmpeg  # CentOS/RHEL
   sudo pacman -S ffmpeg    # Arch Linux
   ```
   
   **macOS:**
   ```bash
   brew install ffmpeg
   ```
   
   **Note:** If ffmpeg is not available, the `make_Ez_movie.py` script will skip MP4 generation but still create a GIF file.

### Running Visualization Scripts

Navigate to the `Py_Visualization/` directory and run the desired scripts:

```bash
cd Py_Visualization

# Generate Ez field evolution movie (MP4 and GIF)
# Requires: numpy, matplotlib, imageio, ffmpeg (for MP4)
python make_Ez_movie.py

# Plot Ez point probe time series
# Requires: numpy, matplotlib
python make_Ez_probe.py

# Plot total EM energy vs time
# Requires: numpy, matplotlib
python make_em_total_energy.py

# Visualize Poynting vector flux
# Requires: numpy, matplotlib
python make_box_poynting.py

# Plot refractive index profile
# Requires: numpy, matplotlib
python make_index_profile.py
```

**Dependencies Summary:**
- **All scripts require:** `numpy`, `matplotlib` (installed via `requirements.txt`)
- **make_Ez_movie.py additionally requires:** `imageio` (for GIF) and `ffmpeg` (for MP4)
- All other imports (`json`, `glob`, `struct`, `pathlib`) are Python standard library

**Important:** The Python scripts read from hardcoded paths relative to the `Py_Visualization/` directory:
- They expect output at: `../out/build/x64-release/frames/3D_FDTD_v008_output/`
- If you build with a different configuration, you may need to update the paths in the Python scripts

### Generated Visualizations

- `ez_evolution.mp4` / `ez_evolution.gif` - Animated field evolution
- `Ez_probe.jpg` - Time series plot of Ez at probe point
- `EM_total_energy.jpg` - Total EM energy vs time
- `refractive_index.jpg` - Refractive index profile
- Other analysis plots as configured

## Configuration

Simulation parameters can be edited in `src/params.hpp`:

- **Grid size**: `Nx`, `Ny`, `Nz` (default: 100×100×100)
- **Grid spacing**: `dx`, `dy`, `dz` (default: 10 nm)
- **Time steps**: `nSteps` (default: 1000)
- **Boundary conditions**: PML thickness, CPML parameters
- **Source parameters**: Frequency, pulse width, position
- **Detector locations**: Probe points, slice planes
- **Output settings**: Save frequency, output directory

After modifying parameters, rebuild the project and run again.

## Workflow Summary

1. **Configure simulation** in `src/params.hpp`
2. **Build the project** using CMake
3. **Run the C++ simulation** to generate output data
4. **Run Python visualization scripts** to create plots and movies
5. **Analyze results** from generated images and data files

## Output Data Format

- **2D slices**: Binary `.raw` files (float64) with metadata in `metadata.json`
- **1D probes**: Binary `.bin` files (float64) with metadata in `metadata_ts.json`
- **Time series**: CSV files with time and value columns
- **Metadata**: JSON files containing grid dimensions, parameters, and file patterns

## Notes

- The simulation uses a Yee grid with staggered electric and magnetic fields
- CFL condition is automatically enforced for stability
- PML boundaries are included in the total grid size
- All output files are written to `frames/` directory relative to the executable location
