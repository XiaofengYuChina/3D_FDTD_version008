#!/usr/bin/env python3
"""
make_Ez_movie_physical.py - Generate Ez field evolution movie with physical coordinates

Features:
- Physical coordinate system (nm)
- Structure boundaries drawn as thick black lines
- Source position marker (uses actual coordinates from mesh_info.json)
- Physical domain boundaries (uses phys_domain_x/y_min/max_nm from metadata)
- Step/time annotation on each frame (step=N  t=X.XXXX fs)
- Non-uniform grid support via pcolormesh
"""

import argparse
import json
import glob
import sys
import re
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter
from matplotlib.patches import Rectangle, Circle
from matplotlib.lines import Line2D
import imageio_ffmpeg

from load_grid import load_grid_spacing, find_output_path

matplotlib.rcParams['animation.ffmpeg_path'] = imageio_ffmpeg.get_ffmpeg_exe()


def extract_step_from_filename(filepath):
    """Extract step number from filename like ez_0120.raw -> 120"""
    import os
    basename = os.path.basename(filepath)
    match = re.search(r'ez_(\d+)\.raw', basename)
    if match:
        return int(match.group(1))
    return 0


def draw_structure_boundaries(ax, structures, z_slice_nm, linewidth=2.0, color='black'):
    """Draw structure boundaries at the current z-slice."""
    if structures is None or z_slice_nm is None:
        return []

    artists = []
    for s in structures:
        stype = s['type']
        params = s['params_nm']

        if stype == 'box':
            x0, x1, y0, y1, z0, z1 = params
            if z0 <= z_slice_nm <= z1:
                rect = Rectangle((x0, y0), x1 - x0, y1 - y0,
                                fill=False, edgecolor=color,
                                linewidth=linewidth, linestyle='-', zorder=10)
                ax.add_patch(rect)
                artists.append(rect)

        elif stype == 'cylinder':
            cx, cy, radius, z0, z1, _ = params
            if z0 <= z_slice_nm <= z1:
                circle = Circle((cx, cy), radius,
                               fill=False, edgecolor=color,
                               linewidth=linewidth, linestyle='-', zorder=10)
                ax.add_patch(circle)
                artists.append(circle)

        elif stype == 'sphere':
            cx, cy, cz, radius, _, _ = params
            dz = abs(z_slice_nm - cz)
            if dz < radius:
                r_slice = np.sqrt(radius**2 - dz**2)
                circle = Circle((cx, cy), r_slice,
                               fill=False, edgecolor=color,
                               linewidth=linewidth, linestyle='-', zorder=10)
                ax.add_patch(circle)
                artists.append(circle)

    return artists


def main():
    parser = argparse.ArgumentParser(description='Generate Ez field evolution movie')
    parser.add_argument('--fps', type=int, default=20, help='Frames per second (default: 20)')
    parser.add_argument('--vmin', type=float, default=-100, help='Color scale minimum (default: -100)')
    parser.add_argument('--vmax', type=float, default=100, help='Color scale maximum (default: 100)')
    parser.add_argument('--output', type=str, default='ez_evolution_physical.mp4', help='Output filename')
    args = parser.parse_args()

    # Find output directory
    data_path = find_output_path(subfolder="Ez_center")
    if not data_path.exists():
        print(f"Error: Directory not found: {data_path}")
        sys.exit(1)

    # Load grid
    try:
        grid = load_grid_spacing()
    except FileNotFoundError as e:
        print(f"Error: {e}")
        sys.exit(1)

    # Read Ez metadata
    meta_file = data_path / "metadata.json"
    if not meta_file.exists():
        print(f"Error: metadata.json not found in {data_path}")
        sys.exit(1)

    with open(meta_file, "r", encoding="utf-8") as f:
        meta = json.load(f)

    NxT, NyT = int(meta["NxT"]), int(meta["NyT"])
    npml = int(meta["npml"])

    # Load mesh_info.json (required)
    base_path = find_output_path()
    mesh_info_file = base_path / "mesh_info.json"
    if not mesh_info_file.exists():
        print(f"Error: mesh_info.json not found in {base_path}")
        sys.exit(1)

    with open(mesh_info_file, "r", encoding="utf-8") as f:
        mesh_info = json.load(f)

    structures = mesh_info.get("structures", [])
    grid_type = mesh_info.get("grid_type", "bounds")

    # Validate grid dimensions
    if grid_type == "bounds":
        if len(grid["x"]) != NxT + 1 or len(grid["y"]) != NyT + 1:
            print(f"Error: Grid dimension mismatch")
            sys.exit(1)

    # Source position (actual coordinates)
    source_x_actual_nm = mesh_info.get("source_x_actual_nm")
    source_y_actual_nm = mesh_info.get("source_y_actual_nm")

    # Physical domain boundaries
    phys_domain_x_min_nm = mesh_info.get("phys_domain_x_min_nm", 0.0)
    phys_domain_x_max_nm = mesh_info.get("phys_domain_x_max_nm")
    phys_domain_y_min_nm = mesh_info.get("phys_domain_y_min_nm", 0.0)
    phys_domain_y_max_nm = mesh_info.get("phys_domain_y_max_nm")

    if phys_domain_x_max_nm is None:
        print("Error: phys_domain_x_max_nm not found in mesh_info.json")
        sys.exit(1)

    # Build coordinate arrays
    if grid_type == "bounds":
        pml_offset_x = grid["x"][npml]
        pml_offset_y = grid["y"][npml]
        x_nm = (grid["x"] - pml_offset_x) * 1e9
        y_nm = (grid["y"] - pml_offset_y) * 1e9
    else:
        x_center = grid.get("x_center", grid["x"])
        y_center = grid.get("y_center", grid["y"])
        pml_offset_x = x_center[npml] - grid["dx"][npml] / 2
        pml_offset_y = y_center[npml] - grid["dy"][npml] / 2
        x_center_nm = (x_center - pml_offset_x) * 1e9
        y_center_nm = (y_center - pml_offset_y) * 1e9
        half_dx = np.concatenate([[grid["dx"][0]/2], grid["dx"]/2])
        half_dy = np.concatenate([[grid["dy"][0]/2], grid["dy"]/2])
        x_nm = np.concatenate([[x_center_nm[0] - half_dx[0]*1e9],
                               x_center_nm + half_dx[1:]*1e9])
        y_nm = np.concatenate([[y_center_nm[0] - half_dy[0]*1e9],
                               y_center_nm + half_dy[1:]*1e9])

    # Create meshgrid for pcolormesh
    X_mesh, Y_mesh = np.meshgrid(x_nm, y_nm)

    # Get z-slice info
    z_slice_physical_nm = meta.get("z_slice_physical_nm")

    # Load frame files
    frame_files = sorted(glob.glob(str(data_path / "ez_*.raw")))
    if len(frame_files) == 0:
        print("Error: No frame files found!")
        sys.exit(1)

    print(f"Found {len(frame_files)} frames")

    # Validate first frame shape
    first_frame = np.fromfile(frame_files[0], dtype=np.float64)
    if first_frame.size != NxT * NyT:
        print(f"Error: Frame size mismatch ({first_frame.size} vs {NxT * NyT})")
        sys.exit(1)

    # Load dt for time annotation
    dt_fs = None
    if "dt_fs" in mesh_info:
        dt_fs = mesh_info["dt_fs"]
    elif "dt" in mesh_info:
        dt_fs = mesh_info["dt"] * 1e15
    elif "dt_fs" in meta:
        dt_fs = meta["dt_fs"]
    elif "dt" in meta:
        dt_fs = meta["dt"] * 1e15
    if dt_fs is None:
        dt_fs = 0.01  # Fallback

    # Create figure
    fig, ax = plt.subplots(figsize=(8, 7))

    z_data = first_frame.reshape(NxT, NyT).T
    im = ax.pcolormesh(X_mesh, Y_mesh, z_data, cmap="seismic",
                       vmin=args.vmin, vmax=args.vmax, shading='flat')

    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label("Ez (arb. units)")
    ax.set_xlabel("x (nm)")
    ax.set_ylabel("y (nm)")
    ax.set_aspect('equal', adjustable='box')

    if z_slice_physical_nm is not None:
        ax.set_title(f"Ez @ z = {z_slice_physical_nm:.1f} nm")
    else:
        ax.set_title("Ez @ z = center")

    # Physical domain boundary lines
    ax.axvline(phys_domain_x_min_nm, color="k", ls="--", lw=1.5, alpha=0.7)
    ax.axvline(phys_domain_x_max_nm, color="k", ls="--", lw=1.5, alpha=0.7)
    ax.axhline(phys_domain_y_min_nm, color="k", ls="--", lw=1.5, alpha=0.7)
    ax.axhline(phys_domain_y_max_nm, color="k", ls="--", lw=1.5, alpha=0.7)

    # Structure boundaries
    draw_structure_boundaries(ax, structures, z_slice_physical_nm, linewidth=2.5, color='black')

    # Source marker
    if source_x_actual_nm is not None and source_y_actual_nm is not None:
        ax.plot(source_x_actual_nm, source_y_actual_nm, 'x', color='lime', markersize=14,
                markeredgewidth=3, zorder=20)
        ax.axvline(source_x_actual_nm, color='lime', ls=':', lw=1, alpha=0.6, zorder=5)
        ax.axhline(source_y_actual_nm, color='lime', ls=':', lw=1, alpha=0.6, zorder=5)

    # Time annotation text
    first_step = extract_step_from_filename(frame_files[0])
    first_t_fs = first_step * dt_fs
    time_text = ax.text(
        0.02, 0.98, f"step={first_step}   t={first_t_fs:.4f} fs",
        transform=ax.transAxes, ha="left", va="top", fontsize=10, color="white",
        bbox=dict(facecolor="black", alpha=0.6, edgecolor="none"), zorder=30
    )

    # Legend
    legend_elements = [
        Line2D([0], [0], color='k', ls='--', lw=1.5, label='Physical domain'),
    ]
    if source_x_actual_nm is not None:
        legend_elements.append(
            Line2D([0], [0], color='lime', ls=':', lw=1, marker='x', markersize=8,
                   markeredgewidth=2, label=f'Source ({source_x_actual_nm:.0f}, {source_y_actual_nm:.0f}) nm')
        )
    ax.legend(handles=legend_elements, loc='upper right', fontsize=8)

    # Generate MP4
    print(f"Generating {args.output}...")
    try:
        writer = FFMpegWriter(fps=args.fps, codec="libx264", bitrate=4000)
        with writer.saving(fig, args.output, dpi=140):
            for i, fp in enumerate(frame_files):
                z = np.fromfile(fp, dtype=np.float64).reshape(NxT, NyT).T
                im.set_array(z.ravel())

                step = extract_step_from_filename(fp)
                t_fs = step * dt_fs
                time_text.set_text(f"step={step}   t={t_fs:.4f} fs")

                writer.grab_frame()
                if (i + 1) % 50 == 0:
                    print(f"  Progress: {i + 1}/{len(frame_files)}")

        print(f"Saved: {args.output}")
    except Exception as e:
        print(f"Error generating movie: {e}")
    finally:
        plt.close(fig)


if __name__ == "__main__":
    main()
