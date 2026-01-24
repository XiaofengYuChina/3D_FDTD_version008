#!/usr/bin/env python3
"""
plot_point_field.py - Visualize point field detector time series

Usage:
    python plot_point_field.py --name Ez_probe
    python plot_point_field.py --name E_probe --save e_probe.png
    python plot_point_field.py --name Ez_probe --fft --save spectrum.png

Features:
- Reads metadata.json automatically to determine field components
- Plots time series of all recorded components
- Optional FFT spectrum analysis
- Physical time axis (fs or ps)
"""

import argparse
import json
import numpy as np
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from load_grid import find_output_path


def load_metadata(detector_path):
    """Load metadata.json from detector directory."""
    meta_path = detector_path / "metadata.json"
    if not meta_path.exists():
        raise FileNotFoundError(f"Metadata not found: {meta_path}")
    with open(meta_path, 'r') as f:
        return json.load(f)


def load_time_series(detector_path, filename, dtype):
    """Load time series data from binary file."""
    filepath = detector_path / filename
    if not filepath.exists():
        raise FileNotFoundError(f"Data file not found: {filepath}")
    return np.fromfile(filepath, dtype=dtype)


def main():
    parser = argparse.ArgumentParser(description='Plot point field detector time series')
    parser.add_argument('--name', type=str, required=True,
                        help='Detector name (folder name, e.g., Ez_probe)')
    parser.add_argument('--run-tag', type=str, default='3D_FDTD_v008_output',
                        help='Run tag (default: 3D_FDTD_v008_output)')
    parser.add_argument('--fft', action='store_true',
                        help='Show FFT spectrum instead of time series')
    parser.add_argument('--log-fft', action='store_true',
                        help='Use log scale for FFT magnitude')
    parser.add_argument('--xlim', type=float, nargs=2, default=None,
                        help='X-axis limits (e.g., --xlim 0 100)')
    parser.add_argument('--ylim', type=float, nargs=2, default=None,
                        help='Y-axis limits')
    parser.add_argument('--freq-max', type=float, default=None,
                        help='Maximum frequency for FFT plot (THz)')
    parser.add_argument('--save', type=str, default=None,
                        help='Save to file')
    parser.add_argument('--dpi', type=int, default=150,
                        help='DPI for saved figure (default: 150)')
    parser.add_argument('--figsize', type=float, nargs=2, default=[12, 6],
                        help='Figure size in inches (default: 12 6)')
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
    components = metadata.get('components', ['Ez'])
    files = metadata.get('files', [f"{c}_ts.bin" for c in components])
    dt = metadata.get('dt', 1e-18)  # Time step in seconds
    save_every = metadata.get('save_every', 1)
    n_steps = metadata.get('n_steps', 0)
    dtype = np.float64 if metadata.get('dtype', 'float64') == 'float64' else np.float32

    # Probe position
    probe_x = metadata.get('probe_x_physical_nm', 0)
    probe_y = metadata.get('probe_y_physical_nm', 0)
    probe_z = metadata.get('probe_z_physical_nm', 0)

    print(f"Probe position: ({probe_x:.1f}, {probe_y:.1f}, {probe_z:.1f}) nm")
    print(f"Components: {components}")
    print(f"dt = {dt*1e15:.4f} fs, save_every = {save_every}")

    # Load time series data for each component
    data = {}
    for comp, fname in zip(components, files):
        try:
            data[comp] = load_time_series(detector_path, fname, dtype)
            print(f"  {comp}: {len(data[comp])} samples, range [{data[comp].min():.3e}, {data[comp].max():.3e}]")
        except FileNotFoundError as e:
            print(f"Warning: {e}")

    if not data:
        print("Error: No data loaded")
        return 1

    # Create time array
    n_samples = len(list(data.values())[0])
    dt_eff = dt * save_every
    t = np.arange(n_samples) * dt_eff

    # Choose time unit
    t_max = t[-1]
    if t_max < 1e-12:
        t_unit = 'fs'
        t_scale = 1e15
    else:
        t_unit = 'ps'
        t_scale = 1e12

    t_scaled = t * t_scale

    # Create figure
    fig, ax = plt.subplots(figsize=tuple(args.figsize))

    if args.fft:
        # FFT analysis
        freq = np.fft.rfftfreq(n_samples, dt_eff)
        freq_THz = freq * 1e-12  # Convert to THz

        for comp, values in data.items():
            spectrum = np.abs(np.fft.rfft(values))
            if args.log_fft:
                spectrum = np.log10(spectrum + 1e-30)
            ax.plot(freq_THz, spectrum, label=comp, linewidth=1)

        ax.set_xlabel('Frequency (THz)')
        if args.log_fft:
            ax.set_ylabel('log10(|FFT|)')
        else:
            ax.set_ylabel('|FFT|')
        ax.set_title(f'Spectrum at ({probe_x:.0f}, {probe_y:.0f}, {probe_z:.0f}) nm')

        if args.freq_max:
            ax.set_xlim(0, args.freq_max)
        else:
            # Auto limit to meaningful range
            ax.set_xlim(0, min(freq_THz[-1], 1000))

        ax.grid(True, alpha=0.3)
        ax.legend()

    else:
        # Time series plot
        for comp, values in data.items():
            ax.plot(t_scaled, values, label=comp, linewidth=0.8)

        ax.set_xlabel(f'Time ({t_unit})')
        ax.set_ylabel('Field Value')
        ax.set_title(f'Time Series at ({probe_x:.0f}, {probe_y:.0f}, {probe_z:.0f}) nm')
        ax.grid(True, alpha=0.3)
        ax.legend()

        if args.xlim:
            ax.set_xlim(args.xlim)

    if args.ylim:
        ax.set_ylim(args.ylim)

    plt.tight_layout()

    # Save or show
    if args.save:
        plt.savefig(args.save, dpi=args.dpi, bbox_inches='tight')
        print(f"Saved to: {args.save}")
    else:
        suffix = "_spectrum" if args.fft else "_timeseries"
        default_save = detector_path / f"{args.name}{suffix}.png"
        plt.savefig(default_save, dpi=args.dpi, bbox_inches='tight')
        print(f"Saved to: {default_save}")

    plt.close()
    return 0


if __name__ == '__main__':
    exit(main())
