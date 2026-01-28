#!/usr/bin/env python3
"""
plot_mesh_detector.py - Visualize refractive index distribution and mesh grid

Usage:
    python plot_mesh_detector.py --name mesh_z_500nm
    python plot_mesh_detector.py --name mesh_y_2000nm
    python plot_mesh_detector.py --name mesh_z_500nm --show-grid
    python plot_mesh_detector.py --name mesh_y_2000nm --show-grid

Features:
- Reads metadata.json automatically to determine slice plane and dimensions
- Shows refractive index distribution with colorbar
- Optional mesh grid overlay
- Physical coordinate system (nm)
"""

import argparse
import json
import numpy as np
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

from load_grid import find_output_path


def load_metadata(detector_path):
    """Load metadata.json from detector directory."""
    meta_path = detector_path / "metadata.json"
    if not meta_path.exists():
        raise FileNotFoundError(f"Metadata not found: {meta_path}")
    with open(meta_path, 'r') as f:
        return json.load(f)


def load_refractive_index(detector_path, metadata):
    """Load refractive index data from binary file."""
    data_path = detector_path / "refractive_index.raw"
    if not data_path.exists():
        raise FileNotFoundError(f"Data not found: {data_path}")

    dtype = np.float64 if metadata.get('dtype', 'float64') == 'float64' else np.float32
    dim1 = metadata['slice_dim1']
    dim2 = metadata['slice_dim2']

    data = np.fromfile(data_path, dtype=dtype)
    return data.reshape((dim1, dim2))


def load_mesh_bounds(detector_path):
    """Load mesh boundary coordinates."""
    bounds = {}
    for axis in ['x', 'y', 'z']:
        path = detector_path / f"{axis}_bounds.txt"
        if path.exists():
            with open(path, 'r') as f:
                lines = [l.strip() for l in f if l.strip() and not l.startswith('#')]
                bounds[axis] = np.array([float(x) for x in lines]) * 1e9  # Convert to nm
    return bounds


def main():
    parser = argparse.ArgumentParser(description='Plot refractive index and mesh grid')
    parser.add_argument('--name', type=str, default='mesh_info',
                        help='Detector name (folder name, default: mesh_info)')
    parser.add_argument('--run-tag', type=str, default='3D_FDTD_v008_output',
                        help='Run tag (default: 3D_FDTD_v008_output)')
    parser.add_argument('--show-grid', action='store_true', default=True,
                        help='Show mesh grid lines (default: on)')
    parser.add_argument('--no-grid', action='store_true',
                        help='Hide mesh grid lines')
    parser.add_argument('--grid-color', type=str, default='black',
                        help='Grid line color (default: black)')
    parser.add_argument('--grid-alpha', type=float, default=0.5,
                        help='Grid line transparency (default: 0.5)')
    parser.add_argument('--grid-width', type=float, default=0.8,
                        help='Grid line width (default: 0.8)')
    parser.add_argument('--cmap', type=str, default='viridis',
                        help='Colormap (default: viridis)')
    parser.add_argument('--save', type=str, default=None,
                        help='Save to file (e.g., mesh_plot.png)')
    parser.add_argument('--dpi', type=int, default=150,
                        help='DPI for saved figure (default: 150)')
    parser.add_argument('--show', action='store_true',
                        help='Show interactive plot')
    args = parser.parse_args()

    # Find detector directory
    detector_path = find_output_path(run_tag=args.run_tag, subfolder=args.name)
    if not detector_path.exists():
        print(f"Error: Detector directory not found: {detector_path}")
        return 1

    print(f"Loading data from: {detector_path}")

    # Load metadata
    try:
        metadata = load_metadata(detector_path)
    except FileNotFoundError as e:
        print(f"Error: {e}")
        return 1

    # Load refractive index data
    try:
        n_data = load_refractive_index(detector_path, metadata)
    except FileNotFoundError as e:
        print(f"Error: {e}")
        return 1

    # Load mesh bounds
    bounds = load_mesh_bounds(detector_path)

    # Determine axis labels based on slice plane
    slice_plane = metadata.get('slice_plane', 'XY')
    dim1_label = metadata.get('dim1_label', 'x')
    dim2_label = metadata.get('dim2_label', 'y')
    slice_coord = metadata.get('slice_physical_nm', 0)

    print(f"Slice plane: {slice_plane} at {slice_coord:.1f} nm")
    print(f"Data shape: {n_data.shape}")
    print(f"n range: [{n_data.min():.3f}, {n_data.max():.3f}]")

    # Get coordinate arrays for pcolormesh
    if dim1_label in bounds and dim2_label in bounds:
        x_bounds = bounds[dim1_label]
        y_bounds = bounds[dim2_label]
    else:
        # Fallback to uniform grid
        Lx = metadata.get('Lx_phys_m', 1e-6) * 1e9
        Ly = metadata.get('Ly_phys_m', 1e-6) * 1e9
        x_bounds = np.linspace(0, Lx, n_data.shape[0] + 1)
        y_bounds = np.linspace(0, Ly, n_data.shape[1] + 1)

    # Create figure
    fig, ax = plt.subplots(figsize=(10, 8))

    # Plot refractive index using pcolormesh (handles non-uniform grid)
    X, Y = np.meshgrid(x_bounds, y_bounds, indexing='ij')
    im = ax.pcolormesh(X, Y, n_data, cmap=args.cmap, shading='flat')

    # Add colorbar
    cbar = plt.colorbar(im, ax=ax, label='Refractive Index n')

    # Draw mesh grid (default: on, use --no-grid to disable)
    show_grid = args.show_grid and not args.no_grid
    if show_grid:
        # Vertical lines (constant x)
        vlines = [[(x, y_bounds[0]), (x, y_bounds[-1])] for x in x_bounds]
        # Horizontal lines (constant y)
        hlines = [[(x_bounds[0], y), (x_bounds[-1], y)] for y in y_bounds]

        all_lines = vlines + hlines
        lc = LineCollection(all_lines, colors=args.grid_color,
                           alpha=args.grid_alpha, linewidths=args.grid_width)
        ax.add_collection(lc)

    # Labels and title
    ax.set_xlabel(f'{dim1_label.upper()} (nm)')
    ax.set_ylabel(f'{dim2_label.upper()} (nm)')
    ax.set_title(f'Refractive Index - {slice_plane} plane at {dim1_label if slice_plane=="YZ" else (dim2_label if slice_plane=="XZ" else "z")}={slice_coord:.1f} nm')
    ax.set_aspect('equal')

    plt.tight_layout()

    # Save or show
    script_dir = Path(__file__).resolve().parent
    if args.save:
        plt.savefig(args.save, dpi=args.dpi, bbox_inches='tight')
        print(f"Saved to: {args.save}")

    if args.show:
        plt.switch_backend('TkAgg')
        plt.show()
    elif not args.save:
        # Default: save to script directory with detector name
        default_save = script_dir / f"{args.name}_refractive_index_plot.png"
        plt.savefig(default_save, dpi=args.dpi, bbox_inches='tight')
        print(f"Saved to: {default_save}")

    return 0


if __name__ == '__main__':
    exit(main())
