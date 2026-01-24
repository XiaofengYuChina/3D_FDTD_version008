# Tests Directory

This directory contains unit tests and verification tools for the 3D FDTD simulation.

## C++ Tests

- **test_auto_mesh_v2.cpp** - Unit tests for the Auto Non-uniform Mesh Generator V2
  - Tests interface snapping, override regions, grading ratios, PPW mapping
  - Compile: Add to CMakeLists.txt or compile separately with: `g++ -std=c++17 -I../src test_auto_mesh_v2.cpp -o test_auto_mesh`

- **verify_mesh_override.cpp** - Integration test for mesh override alignment
  - Validates user's real case with cylinder structure and override box
  - Ensures fine mesh is confined to the override region

## Python Tests

- **make_Ez_test.py** - Enhanced Ez visualization with logarithmic color scaling
  - Alternative visualization with custom colormap and dark theme
  - Outputs: ez_evolution_enhanced.mp4

## Running Tests

```bash
# Python tests (from Py_Visualization directory)
cd ../Py_Visualization
python ../tests/make_Ez_test.py
```
