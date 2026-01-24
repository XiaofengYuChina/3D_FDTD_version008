"""
load_grid.py - Load and parse grid spacing data from FDTD simulation

Supports multiple mesh modes:
- UNIFORM: Constant spacing in all directions
- MANUAL_GRADED: Center-fine, edge-coarse grading
- AUTO_NONUNIFORM: Structure-aware automatic mesh
"""

import numpy as np
from pathlib import Path
import json


# Build mode selection: "auto", "debug", or "release"
# Change this to force a specific build mode
BUILD_MODE = "auto"


def find_output_path(run_tag="3D_FDTD_v008_output", subfolder=None, build_mode=None):
    """
    Find the output directory (handles different build configurations).

    Parameters:
    -----------
    run_tag : str
        The output folder name (e.g., "3D_FDTD_v008_output")
    subfolder : str, optional
        Subfolder within run_tag (e.g., "Ez_center", "Ez_probe", "n_center")
    build_mode : str, optional
        Build mode: "auto", "debug", or "release". If None, uses global BUILD_MODE.

    Returns:
    --------
    Path : The found path, or first option if none exist
    """
    mode = build_mode if build_mode is not None else BUILD_MODE

    # Define paths based on build mode
    if mode == "debug":
        possible_paths = [
            Path(f"../out/build/x64-debug/Debug/frames/{run_tag}"),
            Path(f"../build/Debug/frames/{run_tag}"),
        ]
    elif mode == "release":
        possible_paths = [
            Path(f"../out/build/x64-release/Release/frames/{run_tag}"),
            Path(f"../build/Release/frames/{run_tag}"),
        ]
    else:  # auto - try all paths
        possible_paths = [
            Path(f"../out/build/x64-debug/Debug/frames/{run_tag}"),
            Path(f"../out/build/x64-release/Release/frames/{run_tag}"),
            Path(f"../build/Release/frames/{run_tag}"),
            Path(f"../build/Debug/frames/{run_tag}"),
            Path(f"../build/frames/{run_tag}"),
            Path(f"../frames/{run_tag}"),
            Path(f"frames/{run_tag}"),
        ]

    for p in possible_paths:
        if p.exists():
            if subfolder:
                return p / subfolder
            return p

    # Return first option even if it doesn't exist
    first = possible_paths[0]
    if subfolder:
        return first / subfolder
    return first


def load_grid_spacing(run_tag="3D_FDTD_v008_output"):
    """
    Load grid spacing and compute physical coordinates.

    Returns a dictionary with:
    - dx, dy, dz: Cell spacing arrays (meters)
    - x, y, z: Cell boundary positions (meters)
    - x_center, y_center, z_center: Cell center positions (meters)
    """
    base_path = find_output_path(run_tag)
    grid_file = base_path / "grid_spacing.txt"

    if not grid_file.exists():
        raise FileNotFoundError(f"Grid spacing file not found: {grid_file}")

    # Parse the file with section markers
    with open(grid_file, 'r') as f:
        lines = f.readlines()

    # Parse sections
    sections = {
        'dx': [], 'dy': [], 'dz': [],
        'x_bounds': [], 'y_bounds': [], 'z_bounds': []
    }
    current_section = 'dx'

    for line in lines:
        line = line.strip()
        if not line:
            continue

        # Check for section markers
        if line.startswith('#'):
            line_lower = line.lower()
            if 'dy' in line_lower and 'x_bounds' not in line_lower and 'y_bounds' not in line_lower and 'z_bounds' not in line_lower:
                current_section = 'dy'
            elif 'dz' in line_lower and 'x_bounds' not in line_lower and 'y_bounds' not in line_lower and 'z_bounds' not in line_lower:
                current_section = 'dz'
            elif 'x_bounds' in line_lower:
                current_section = 'x_bounds'
            elif 'y_bounds' in line_lower:
                current_section = 'y_bounds'
            elif 'z_bounds' in line_lower:
                current_section = 'z_bounds'
            elif 'dx' in line_lower:
                current_section = 'dx'
            continue

        # Parse numeric value
        try:
            val = float(line)
            sections[current_section].append(val)
        except ValueError:
            continue

    # Convert to numpy arrays
    dx = np.array(sections['dx'])
    dy = np.array(sections['dy'])
    dz = np.array(sections['dz'])

    # Use stored bounds if available, otherwise compute from spacing
    if len(sections['x_bounds']) == len(dx) + 1:
        x = np.array(sections['x_bounds'])
    else:
        x = np.zeros(len(dx) + 1)
        for i in range(len(dx)):
            x[i+1] = x[i] + dx[i]

    if len(sections['y_bounds']) == len(dy) + 1:
        y = np.array(sections['y_bounds'])
    else:
        y = np.zeros(len(dy) + 1)
        for j in range(len(dy)):
            y[j+1] = y[j] + dy[j]

    if len(sections['z_bounds']) == len(dz) + 1:
        z = np.array(sections['z_bounds'])
    else:
        z = np.zeros(len(dz) + 1)
        for k in range(len(dz)):
            z[k+1] = z[k] + dz[k]

    return {
        'dx': dx, 'dy': dy, 'dz': dz,
        'x': x, 'y': y, 'z': z,
        'x_center': (x[:-1] + x[1:]) / 2,
        'y_center': (y[:-1] + y[1:]) / 2,
        'z_center': (z[:-1] + z[1:]) / 2
    }


def load_mesh_info(run_tag="3D_FDTD_v008_output"):
    """Load mesh metadata from JSON file"""
    base_path = find_output_path(run_tag)
    json_file = base_path / "mesh_info.json"

    if not json_file.exists():
        return None

    with open(json_file, 'r') as f:
        return json.load(f)


def compute_grading_ratios(spacing):
    """Compute grading ratio between adjacent cells"""
    if len(spacing) < 2:
        return np.array([1.0])

    ratios = np.zeros(len(spacing) - 1)
    for i in range(len(spacing) - 1):
        ratios[i] = max(spacing[i+1] / spacing[i], spacing[i] / spacing[i+1])
    return ratios


def print_grid_summary(grid, mesh_info=None):
    """Print a summary of the grid properties"""
    dx, dy, dz = grid['dx'], grid['dy'], grid['dz']
    x, y, z = grid['x'], grid['y'], grid['z']

    print("\n" + "="*60)
    print("MESH SUMMARY")
    print("="*60)

    if mesh_info:
        print(f"Mode: {mesh_info.get('mode', 'Unknown')}")
        print(f"Lambda: {mesh_info.get('lambda0_nm', 'N/A')} nm")
        print(f"Structures: {mesh_info.get('num_structures', 'N/A')}")

    print(f"\nDimensions: {len(dx)} x {len(dy)} x {len(dz)} cells")
    print(f"Domain: {x[-1]*1e9:.1f} x {y[-1]*1e9:.1f} x {z[-1]*1e9:.1f} nm")

    print("\nSpacing ranges (nm):")
    print(f"  dx: {dx.min()*1e9:.2f} - {dx.max()*1e9:.2f} (ratio: {dx.max()/dx.min():.2f}x)")
    print(f"  dy: {dy.min()*1e9:.2f} - {dy.max()*1e9:.2f} (ratio: {dy.max()/dy.min():.2f}x)")
    print(f"  dz: {dz.min()*1e9:.2f} - {dz.max()*1e9:.2f} (ratio: {dz.max()/dz.min():.2f}x)")

    # Grading ratios
    rx = compute_grading_ratios(dx)
    ry = compute_grading_ratios(dy)
    rz = compute_grading_ratios(dz)

    print("\nMax grading ratios:")
    print(f"  dx: {rx.max():.4f}", "(OK)" if rx.max() < 1.2 else "(WARNING: >1.2)")
    print(f"  dy: {ry.max():.4f}", "(OK)" if ry.max() < 1.2 else "(WARNING: >1.2)")
    print(f"  dz: {rz.max():.4f}", "(OK)" if rz.max() < 1.2 else "(WARNING: >1.2)")

    print("="*60 + "\n")


# Example usage
if __name__ == "__main__":
    try:
        grid = load_grid_spacing()
        mesh_info = load_mesh_info()
        print_grid_summary(grid, mesh_info)
    except FileNotFoundError as e:
        print(f"Error: {e}")
        print("Make sure you've run the simulation first.")
