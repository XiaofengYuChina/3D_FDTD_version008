#!/usr/bin/env python3
"""
plot_field_movie_2d.py - Visualize 2D field evolution as movie or snapshots

Usage:
    python plot_field_movie_2d.py --name Ez_movie
    python plot_field_movie_2d.py --name Ez_movie --movie --fps 30
    python plot_field_movie_2d.py --name Ex_movie --frame 100 --save ex_frame100.png

Features:
- Reads metadata.json automatically to determine field component and slice
- Generates movie or individual frame snapshots
- Physical coordinate system (nm)
- Automatic color scale or user-defined range
"""

import argparse
import json
import glob
import re
import numpy as np
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter

try:
    import imageio_ffmpeg
    matplotlib.rcParams['animation.ffmpeg_path'] = imageio_ffmpeg.get_ffmpeg_exe()
except ImportError:
    pass

from load_grid import find_output_path, load_grid_spacing


def load_metadata(detector_path):
    """Load metadata.json from detector directory."""
    meta_path = detector_path / "metadata.json"
    if not meta_path.exists():
        raise FileNotFoundError(f"Metadata not found: {meta_path}")
    with open(meta_path, 'r') as f:
        return json.load(f)


def load_frame(filepath, dim1, dim2, dtype):
    """Load a single frame from binary file."""
    data = np.fromfile(filepath, dtype=dtype)
    return data.reshape((dim1, dim2))


def extract_frame_number(filepath):
    """Extract frame number from filename."""
    basename = Path(filepath).name
    match = re.search(r'(\d+)', basename)
    if match:
        return int(match.group(1))
    return 0


def get_frame_files(detector_path, pattern):
    """Get sorted list of frame files."""
    # Convert pattern like "frame_%04d.raw" to glob pattern "frame_*.raw"
    glob_pattern = re.sub(r'%\d*d', '*', pattern)
    files = sorted(glob.glob(str(detector_path / glob_pattern)))
    return files


def main():
    parser = argparse.ArgumentParser(description='Plot 2D field movie or snapshots')
    parser.add_argument('--name', type=str, required=True,
                        help='Detector name (folder name, e.g., Ez_movie)')
    parser.add_argument('--run-tag', type=str, default='3D_FDTD_v008_output',
                        help='Run tag (default: 3D_FDTD_v008_output)')
    parser.add_argument('--movie', action='store_true',
                        help='Generate movie instead of single frame')
    parser.add_argument('--frame', type=int, default=-1,
                        help='Frame number to plot (-1 = last frame)')
    parser.add_argument('--fps', type=int, default=20,
                        help='Frames per second for movie (default: 20)')
    parser.add_argument('--vmin', type=float, default=None,
                        help='Color scale minimum (auto if not set)')
    parser.add_argument('--vmax', type=float, default=None,
                        help='Color scale maximum (auto if not set)')
    parser.add_argument('--cmap', type=str, default='RdBu_r',
                        help='Colormap (default: RdBu_r)')
    parser.add_argument('--save', type=str, default=None,
                        help='Save to file')
    parser.add_argument('--dpi', type=int, default=150,
                        help='DPI for saved figure (default: 150)')
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

    # Extract info from metadata
    field_component = metadata.get('field_component', 'Unknown')
    slice_plane = metadata.get('slice_plane', 'XY')
    slice_nm = metadata.get('slice_physical_nm', 0)
    dim1 = metadata['slice_dim1']
    dim2 = metadata['slice_dim2']
    dim1_label = metadata.get('dim1_label', 'x')
    dim2_label = metadata.get('dim2_label', 'y')
    dtype = np.float64 if metadata.get('dtype', 'float64') == 'float64' else np.float32
    frame_pattern = metadata.get('frame_pattern', 'frame_%04d.raw')
    save_every = metadata.get('save_every', 1)

    print(f"Field: {field_component}, Plane: {slice_plane} at {slice_nm:.1f} nm")
    print(f"Dimensions: {dim1} x {dim2}")

    # Get frame files
    frame_files = get_frame_files(detector_path, frame_pattern)
    if not frame_files:
        print(f"Error: No frame files found matching pattern: {frame_pattern}")
        return 1

    print(f"Found {len(frame_files)} frames")

    # Load grid for physical coordinates
    try:
        grid = load_grid_spacing(args.run_tag)
        # Get bounds based on slice plane
        if dim1_label == 'x':
            x_bounds = grid['x'] * 1e9  # to nm
        elif dim1_label == 'y':
            x_bounds = grid['y'] * 1e9
        else:
            x_bounds = grid['z'] * 1e9

        if dim2_label == 'y':
            y_bounds = grid['y'] * 1e9
        elif dim2_label == 'z':
            y_bounds = grid['z'] * 1e9
        else:
            y_bounds = grid['x'] * 1e9
    except Exception as e:
        print(f"Warning: Could not load grid, using uniform: {e}")
        Lx = metadata.get('Lx_phys_m', 1e-6) * 1e9
        Ly = metadata.get('Ly_phys_m', 1e-6) * 1e9
        x_bounds = np.linspace(0, Lx, dim1 + 1)
        y_bounds = np.linspace(0, Ly, dim2 + 1)

    # Ensure bounds match dimensions
    if len(x_bounds) != dim1 + 1:
        x_bounds = np.linspace(x_bounds[0], x_bounds[-1], dim1 + 1)
    if len(y_bounds) != dim2 + 1:
        y_bounds = np.linspace(y_bounds[0], y_bounds[-1], dim2 + 1)

    X, Y = np.meshgrid(x_bounds, y_bounds, indexing='ij')

    if args.movie:
        # Generate movie
        output_file = args.save or str(detector_path / f"{args.name}_movie.mp4")

        fig, ax = plt.subplots(figsize=(10, 8))

        # Load first frame to set up plot
        data = load_frame(frame_files[0], dim1, dim2, dtype)

        # Determine color scale
        if args.vmin is None or args.vmax is None:
            # Sample frames to determine scale
            sample_indices = np.linspace(0, len(frame_files)-1, min(10, len(frame_files)), dtype=int)
            all_max = 0
            for idx in sample_indices:
                d = load_frame(frame_files[idx], dim1, dim2, dtype)
                all_max = max(all_max, np.abs(d).max())
            vmin = args.vmin if args.vmin is not None else -all_max
            vmax = args.vmax if args.vmax is not None else all_max
        else:
            vmin, vmax = args.vmin, args.vmax

        im = ax.pcolormesh(X, Y, data, cmap=args.cmap, vmin=vmin, vmax=vmax, shading='flat')
        cbar = plt.colorbar(im, ax=ax, label=field_component)
        ax.set_xlabel(f'{dim1_label.upper()} (nm)')
        ax.set_ylabel(f'{dim2_label.upper()} (nm)')
        title = ax.set_title(f'{field_component} - {slice_plane} plane - Frame 0')
        ax.set_aspect('equal')
        plt.tight_layout()

        # Create movie
        print(f"Generating movie: {output_file}")
        writer = FFMpegWriter(fps=args.fps, metadata={'title': f'{field_component} Evolution'})

        with writer.saving(fig, output_file, dpi=args.dpi):
            for i, fpath in enumerate(frame_files):
                data = load_frame(fpath, dim1, dim2, dtype)
                im.set_array(data.ravel())
                frame_num = extract_frame_number(fpath)
                title.set_text(f'{field_component} - {slice_plane} plane - Step {frame_num * save_every}')
                writer.grab_frame()

                if (i + 1) % 50 == 0:
                    print(f"  Processed {i+1}/{len(frame_files)} frames...")

        print(f"Movie saved to: {output_file}")
        plt.close()

    else:
        # Plot single frame
        if args.frame < 0:
            frame_idx = len(frame_files) - 1
        else:
            frame_idx = min(args.frame, len(frame_files) - 1)

        frame_path = frame_files[frame_idx]
        frame_num = extract_frame_number(frame_path)
        data = load_frame(frame_path, dim1, dim2, dtype)

        print(f"Plotting frame {frame_idx} (step {frame_num * save_every})")
        print(f"Field range: [{data.min():.3e}, {data.max():.3e}]")

        fig, ax = plt.subplots(figsize=(10, 8))

        vmin = args.vmin if args.vmin is not None else -np.abs(data).max()
        vmax = args.vmax if args.vmax is not None else np.abs(data).max()

        im = ax.pcolormesh(X, Y, data, cmap=args.cmap, vmin=vmin, vmax=vmax, shading='flat')
        cbar = plt.colorbar(im, ax=ax, label=field_component)
        ax.set_xlabel(f'{dim1_label.upper()} (nm)')
        ax.set_ylabel(f'{dim2_label.upper()} (nm)')
        ax.set_title(f'{field_component} - {slice_plane} plane at z={slice_nm:.1f} nm - Step {frame_num * save_every}')
        ax.set_aspect('equal')
        plt.tight_layout()

        output_file = args.save or str(detector_path / f"{args.name}_frame{frame_idx:04d}.png")
        plt.savefig(output_file, dpi=args.dpi, bbox_inches='tight')
        print(f"Saved to: {output_file}")
        plt.close()

    return 0


if __name__ == '__main__':
    exit(main())
