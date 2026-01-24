#!/usr/bin/env python3
"""
visual_gain_dynamics.py - Visualize population inversion dynamics

VOLUME DENSITY FORMULATION:
  The CSV contains integrated values (atoms), computed from densities [m^-3]
  by integrating over cell volumes: total_Nu = ∫ Nu(r) dV

Definitions:
  inv_sym = (Nu - Ng) / (Nu + Ng)   -- symmetric normalized inversion, range [-1, 1]
  inv_01  = Nu / (Nu + Ng)          -- 0-to-1 population fraction, range [0, 1]

Note: This script recalculates inv_sym from Nu and Ng for accuracy.
"""

import numpy as np
import pandas as pd
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from load_grid import find_output_path


def main():
    # Find data directory
    data_dir = find_output_path(subfolder="populations")
    csv_file = data_dir / 'inversion_time.csv'

    if not csv_file.exists():
        print(f"Error: inversion_time.csv not found in {data_dir}")
        exit(1)

    # Read CSV
    df = pd.read_csv(csv_file)
    df.columns = df.columns.str.strip()

    times_fs = df['t'].values * 1e15
    inversion = df['inversion'].values
    Nu = df['total_Nu'].values
    Ng = df['total_Ng'].values

    # Recalculate normalized inversion (bypass potential C++ threshold issues)
    Ntotal = Nu + Ng
    inv_sym = np.where(Ntotal > 1e-30, (Nu - Ng) / Ntotal, 0.0)
    inv_01 = np.where(Ntotal > 1e-30, Nu / Ntotal, 0.5)

    print(f"Loaded {len(df)} data points from {csv_file}")
    print(f"  Inversion range: {inversion.min():.3e} to {inversion.max():.3e}")
    print(f"  inv_sym range: [{inv_sym.min():.4f}, {inv_sym.max():.4f}]")

    # Create plots
    fig, axes = plt.subplots(3, 1, figsize=(12, 9))

    # Plot 1: Inversion (Nu - Ng)
    axes[0].plot(times_fs, inversion, 'b-', linewidth=2)
    axes[0].axhline(0, color='r', linestyle='--', alpha=0.5)
    axes[0].set_ylabel('Integrated Inversion (atoms)')
    axes[0].set_title('Population Inversion Dynamics', fontweight='bold')
    axes[0].grid(True, alpha=0.3)

    if inversion[0] > 0:
        axes[0].text(0.5, 0.95, 'POSITIVE INVERSION (GAIN)',
                     transform=axes[0].transAxes, fontsize=14, ha='center',
                     bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.8))

    # Plot 2: Populations (integrated: ∫N dV)
    axes[1].plot(times_fs, Nu, 'r-', label='Nu (upper)', linewidth=2)
    axes[1].plot(times_fs, Ng, 'b-', label='Ng (ground)', linewidth=2)
    axes[1].set_ylabel('Integrated Population (atoms)')
    axes[1].set_title('Upper and Ground State Populations', fontweight='bold')
    axes[1].legend(fontsize=11)
    axes[1].grid(True, alpha=0.3)

    # Plot 3: Normalized Inversion
    axes[2].plot(times_fs, inv_sym, 'g-', linewidth=2)
    axes[2].axhline(0, color='r', linestyle='--', alpha=0.5)
    axes[2].set_ylabel('Normalized Inversion')
    axes[2].set_xlabel('Time (fs)')
    axes[2].set_title('Symmetric Normalized Inversion: (Nu-Ng)/(Nu+Ng)', fontweight='bold')
    axes[2].grid(True, alpha=0.3)

    margin = max(0.05, 0.1 * (inv_sym.max() - inv_sym.min()))
    axes[2].set_ylim(inv_sym.min() - margin, inv_sym.max() + margin)

    plt.tight_layout()
    plt.savefig('population_dynamics.png', dpi=150, bbox_inches='tight')
    print("Saved: population_dynamics.png")
    plt.close()

    # Save dedicated inv_sym figure
    fig2, ax2 = plt.subplots(figsize=(10, 6))
    ax2.plot(times_fs, inv_sym, 'g-', linewidth=2, label='inv_sym = (Nu-Ng)/(Nu+Ng)')
    ax2.axhline(0, color='r', linestyle='--', alpha=0.5, label='zero')

    ax2.annotate(f'Initial: {inv_sym[0]:.4f}', xy=(times_fs[0], inv_sym[0]),
                 xytext=(times_fs[0] + 2, inv_sym[0] + 0.02), fontsize=10,
                 bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    ax2.annotate(f'Final: {inv_sym[-1]:.4f}', xy=(times_fs[-1], inv_sym[-1]),
                 xytext=(times_fs[-1] - 5, inv_sym[-1] + 0.02), fontsize=10,
                 bbox=dict(boxstyle='round', facecolor='white', alpha=0.8), ha='right')

    ax2.set_xlabel('Time (fs)')
    ax2.set_ylabel('Symmetric Normalized Inversion')
    ax2.set_title('Symmetric Normalized Inversion: (Nu-Ng)/(Nu+Ng)', fontweight='bold')
    ax2.legend(fontsize=11, loc='best')
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim(inv_sym.min() - margin, inv_sym.max() + margin)

    fig2.savefig('tls_normalized_inversion_sym.png', dpi=150, bbox_inches='tight')
    print("Saved: tls_normalized_inversion_sym.png")
    plt.close(fig2)

    # Summary
    print(f"\nResults:")
    print(f"  Inversion: {inversion[0]:.3e} -> {inversion[-1]:.3e}")
    print(f"  inv_sym: {inv_sym[0]:.4f} -> {inv_sym[-1]:.4f}")
    print(f"  Status: {'Nu > Ng (gain)' if inversion[0] > 0 else 'Ng > Nu (absorption)'}")


if __name__ == "__main__":
    main()
