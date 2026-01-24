"""
make_index_profile.py - Visualize refractive index profile with mesh grid overlay

Features:
- Refractive index distribution at center z-slice
- Mesh grid overlay (dashed lines)
- PML regions highlighted in semi-transparent purple
- Physical/PML boundary marked with thick solid lines
- Physical coordinates (nm) on axes
"""

import matplotlib
matplotlib.use('Agg')  # Non-interactive backend - save only, no display

import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.collections import LineCollection
from pathlib import Path
from load_grid import find_output_path  # Use shared path detection with BUILD_MODE support


def load_grid_spacing(base_path):
    """Load grid spacing data from grid_spacing.txt"""
    grid_file = base_path / "grid_spacing.txt"

    if not grid_file.exists():
        return None

    with open(grid_file, 'r') as f:
        lines = f.readlines()

    sections = {
        'dx': [], 'dy': [], 'dz': [],
        'x_bounds': [], 'y_bounds': [], 'z_bounds': []
    }
    current_section = 'dx'

    for line in lines:
        line = line.strip()
        if not line:
            continue

        if line.startswith('#'):
            line_lower = line.lower()
            if 'dy' in line_lower and 'bounds' not in line_lower:
                current_section = 'dy'
            elif 'dz' in line_lower and 'bounds' not in line_lower:
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

        try:
            val = float(line)
            sections[current_section].append(val)
        except ValueError:
            continue

    dx = np.array(sections['dx'])
    dy = np.array(sections['dy'])
    dz = np.array(sections['dz'])

    # Compute bounds from spacing if not provided
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

    return {'dx': dx, 'dy': dy, 'dz': dz, 'x': x, 'y': y, 'z': z}


def plot_index_with_mesh(n_xy, grid, npml, kslice,
                         mesh_step=1,
                         show_mesh=True,
                         show_pml=True,
                         figsize=(10, 8),
                         dpi=150,
                         z_slice_physical_nm=None,
                         Lx_phys_nm=None, Ly_phys_nm=None):
    """
    Plot refractive index profile with mesh grid overlay.

    Parameters:
    -----------
    n_xy : 2D array
        Refractive index data (NxT x NyT)
    grid : dict
        Grid spacing data with 'x', 'y' arrays (in meters)
    npml : int
        Number of PML cells
    kslice : int
        Z-slice index
    mesh_step : int
        Draw every nth mesh line (1 = all lines, 5 = every 5th line)
    show_mesh : bool
        Whether to show mesh grid lines
    show_pml : bool
        Whether to highlight PML regions
    z_slice_physical_nm : float or None
        Physical z-coordinate of slice in nm (if available from metadata)
    Lx_phys_nm : float or None
        User-configured physical domain size in X (nm)
    Ly_phys_nm : float or None
        User-configured physical domain size in Y (nm)
    """

    NxT, NyT = n_xy.shape

    # Convert to nm for display, with PHYSICAL coordinates (0 at physical domain boundary)
    # Physical coordinate = Total coordinate - PML offset
    pml_offset_x = grid['x'][npml]  # PML offset in meters
    pml_offset_y = grid['y'][npml]

    # Convert to physical coordinates (nm), 0 at physical domain boundary
    x_nm = (grid['x'] - pml_offset_x) * 1e9
    y_nm = (grid['y'] - pml_offset_y) * 1e9

    # Physical domain boundaries
    # Use user-configured size if available, otherwise use grid-derived size
    x_phys_start = 0.0  # x_nm[npml] = 0 after offset
    y_phys_start = 0.0  # y_nm[npml] = 0 after offset

    if Lx_phys_nm is not None and Lx_phys_nm > 0:
        x_phys_end = Lx_phys_nm  # Use user-configured size
    else:
        x_phys_end = x_nm[NxT - npml]  # Fallback to grid-derived size

    if Ly_phys_nm is not None and Ly_phys_nm > 0:
        y_phys_end = Ly_phys_nm  # Use user-configured size
    else:
        y_phys_end = y_nm[NyT - npml]  # Fallback to grid-derived size

    # Full domain boundaries (PML regions are now negative)
    x_min, x_max = x_nm[0], x_nm[-1]
    y_min, y_max = y_nm[0], y_nm[-1]

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Plot refractive index using pcolormesh for non-uniform grid
    X, Y = np.meshgrid(x_nm, y_nm)
    im = ax.pcolormesh(X, Y, n_xy.T, shading='flat', cmap='viridis')

    # Add colorbar
    cbar = plt.colorbar(im, ax=ax, label='Refractive Index n', shrink=0.9)

    # === Highlight PML regions (semi-transparent purple) ===
    if show_pml:
        pml_color = 'purple'
        pml_alpha = 0.25

        # Left PML (negative x region)
        ax.fill_betweenx([y_min, y_max], x_min, x_phys_start,
                         color=pml_color, alpha=pml_alpha, zorder=2)
        # Right PML
        ax.fill_betweenx([y_min, y_max], x_phys_end, x_max,
                         color=pml_color, alpha=pml_alpha, zorder=2)
        # Bottom PML (negative y region, excluding corners already covered)
        ax.fill_between([x_phys_start, x_phys_end], y_min, y_phys_start,
                        color=pml_color, alpha=pml_alpha, zorder=2)
        # Top PML (excluding corners already covered)
        ax.fill_between([x_phys_start, x_phys_end], y_phys_end, y_max,
                        color=pml_color, alpha=pml_alpha, zorder=2)

    # === Draw mesh grid lines (dashed) ===
    if show_mesh:
        mesh_color = 'white'
        mesh_alpha = 0.4
        mesh_linewidth = 0.1

        # Vertical lines (x = const)
        for i in range(0, len(x_nm), mesh_step):
            ax.axvline(x_nm[i], color=mesh_color, linestyle='-',
                      linewidth=mesh_linewidth, alpha=mesh_alpha, zorder=3)

        # Horizontal lines (y = const)
        for j in range(0, len(y_nm), mesh_step):
            ax.axhline(y_nm[j], color=mesh_color, linestyle='-',
                      linewidth=mesh_linewidth, alpha=mesh_alpha, zorder=3)

    # === Draw physical/PML boundary (thick solid line) ===
    boundary_color = 'red'
    boundary_linewidth = 2.0

    # Physical region boundary (rectangle) - now starts at (0, 0)
    rect = plt.Rectangle((x_phys_start, y_phys_start),
                         x_phys_end - x_phys_start,
                         y_phys_end - y_phys_start,
                         fill=False, edgecolor=boundary_color,
                         linewidth=boundary_linewidth, linestyle='-',
                         zorder=4, label='Physical Domain')
    ax.add_patch(rect)

    # Full domain boundary (outer rectangle)
    rect_outer = plt.Rectangle((x_min, y_min),
                               x_max - x_min,
                               y_max - y_min,
                               fill=False, edgecolor='black',
                               linewidth=1.5, linestyle='-',
                               zorder=4, label='Full Domain (incl. PML)')
    ax.add_patch(rect_outer)

    # === Labels and title ===
    ax.set_xlabel('x (nm) [Physical coordinates]', fontsize=12)
    ax.set_ylabel('y (nm) [Physical coordinates]', fontsize=12)

    # Title with physical z-coordinate if available
    if z_slice_physical_nm is not None:
        title = f'Refractive Index Profile @ z = {z_slice_physical_nm:.1f} nm\n'
    else:
        title = f'Refractive Index Profile @ z-slice k={kslice}\n'
    title += f'Grid: {NxT} × {NyT} cells, PML: {npml} cells'
    ax.set_title(title, fontsize=12, fontweight='bold')

    # Set axis limits
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.set_aspect('equal')

    # Add legend
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D

    legend_elements = [
        Patch(facecolor='purple', alpha=0.25, label='CPML Region'),
        Line2D([0], [0], color='red', linewidth=2, label='Physical Domain'),
        Line2D([0], [0], color='black', linewidth=1.5, label='Full Domain'),
    ]
    if show_mesh:
        legend_elements.append(
            Line2D([0], [0], color='white', linewidth=0.5, linestyle='--',
                   label=f'Mesh Grid (every {mesh_step})'))

    ax.legend(handles=legend_elements, loc='upper right', fontsize=9,
              framealpha=0.9)

    plt.tight_layout()
    return fig, ax


def main():
    # === Configuration ===
    run_tag = "3D_FDTD_v008_output"
    mesh_step = 1      # Draw every nth mesh line (increase for dense grids)
    show_mesh = True   # Toggle mesh grid display
    show_pml = True    # Toggle PML region highlight

    # Find output path
    base_path = find_output_path(run_tag)
    root = base_path / "n_center"

    print(f"Looking for data in: {root}")

    if not root.exists():
        print(f"Error: Directory not found: {root}")
        print("Please run the FDTD simulation first.")
        return

    # Read metadata
    meta_file = root / "metadata.json"
    if not meta_file.exists():
        print(f"Error: metadata.json not found in {root}")
        return

    with open(meta_file, "r", encoding="utf-8") as f:
        meta = json.load(f)

    NxT = int(meta["NxT"])
    NyT = int(meta["NyT"])
    NzT = int(meta["NzT"])
    kslice = int(meta["kslice"])
    npml = int(meta["npml"])
    dtype = np.float64 if meta["dtype"] == "float64" else np.float32
    pattern = meta["framePattern"]

    # Get physical z-coordinate if available (new metadata field)
    z_slice_physical_nm = meta.get("z_slice_physical_nm", None)

    # Get user-configured physical domain size (for correct boundary display)
    Lx_phys_nm = meta.get("Lx_phys_nm", None)
    Ly_phys_nm = meta.get("Ly_phys_nm", None)
    Lz_phys_nm = meta.get("Lz_phys_nm", None)

    print(f"Grid: {NxT} x {NyT} x {NzT}")
    print(f"K-slice index: {kslice}")
    if z_slice_physical_nm is not None:
        print(f"Z-slice physical coordinate: {z_slice_physical_nm:.2f} nm")
    print(f"PML thickness: {npml} cells")
    if Lx_phys_nm is not None:
        print(f"Physical domain size: {Lx_phys_nm:.1f} x {Ly_phys_nm:.1f} x {Lz_phys_nm:.1f} nm")

    # Find and load raw data file
    raw_path = root / (pattern.replace("%04d", "0000"))
    if not raw_path.exists():
        candidates = sorted(root.glob("n_*.raw"))
        if not candidates:
            print("Error: Cannot find n_*.raw files")
            return
        raw_path = candidates[0]

    print(f"Loading: {raw_path}")

    data = np.fromfile(raw_path, dtype=dtype)
    if data.size != NxT * NyT:
        print(f"Error: Size mismatch - loaded {data.size}, expected {NxT * NyT}")
        return

    n_xy = data.reshape((NxT, NyT))

    # Load grid spacing
    grid = load_grid_spacing(base_path)

    if grid is None:
        print("Warning: grid_spacing.txt not found, using uniform grid assumption")
        # Create uniform grid as fallback
        dx_uniform = 10e-9  # Assume 10 nm
        grid = {
            'x': np.arange(NxT + 1) * dx_uniform,
            'y': np.arange(NyT + 1) * dx_uniform,
            'dx': np.ones(NxT) * dx_uniform,
            'dy': np.ones(NyT) * dx_uniform,
        }

    # Auto-adjust mesh_step for very dense grids
   # if NxT > 200 and mesh_step == 1:
   #     mesh_step = max(1, NxT // 50)
   #     print(f"Auto-adjusted mesh_step to {mesh_step} for readability")

    # Create plot
    print("Generating plot...")
    fig, ax = plot_index_with_mesh(
        n_xy, grid, npml, kslice,
        mesh_step=mesh_step,
        show_mesh=show_mesh,
        show_pml=show_pml,
        figsize=(10, 8),
        dpi=150,
        z_slice_physical_nm=z_slice_physical_nm,
        Lx_phys_nm=Lx_phys_nm,
        Ly_phys_nm=Ly_phys_nm
    )

    # Save figure
    output_file = "refractive_index_with_mesh.png"
    fig.savefig(output_file, dpi=1080, bbox_inches='tight')
    print(f"Saved: {output_file}")

    # Also save a version without mesh for comparison
    fig2, ax2 = plot_index_with_mesh(
        n_xy, grid, npml, kslice,
        mesh_step=mesh_step,
        show_mesh=False,
        show_pml=show_pml,
        figsize=(10, 8),
        dpi=150,
        z_slice_physical_nm=z_slice_physical_nm,
        Lx_phys_nm=Lx_phys_nm,
        Ly_phys_nm=Ly_phys_nm
    )
    fig2.savefig("refractive_index_no_mesh.png", dpi=1080, bbox_inches='tight')
    print("Saved: refractive_index_no_mesh.png")
    plt.close('all')


if __name__ == "__main__":
    main()
