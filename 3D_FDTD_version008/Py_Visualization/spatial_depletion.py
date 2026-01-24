#!/usr/bin/env python3
"""
spatial_depletion.py - Memory-efficient spatial gain depletion analysis
Loads only first & last frames for comparison
"""

import matplotlib
matplotlib.use('Agg')  # Non-interactive backend - save only, no display

import numpy as np
import matplotlib.pyplot as plt
import json
from pathlib import Path
from load_grid import find_output_path

print("Spatial gain depletion analysis")

# Find data directory (auto-detect debug/release)
data_dir = find_output_path(subfolder="populations")
print(f"Looking for data in: {data_dir}")

if not data_dir.exists():
    print(f"Error: Directory not found: {data_dir}")
    print("Please run the FDTD simulation first.")
    exit(1)

# Load metadata
meta_file = data_dir / 'metadata_populations.json'
if not meta_file.exists():
    print(f"Error: metadata_populations.json not found in {data_dir}")
    exit(1)

with open(meta_file) as f:
    meta = json.load(f)

NxT, NyT, NzT = meta['NxT'], meta['NyT'], meta['NzT']
grid_size = NxT * NyT * NzT
n_timesteps = meta['nSteps'] // meta['saveEvery']

print(f"Grid: {NxT} x {NyT} x {NzT}")
print(f"Timesteps: {n_timesteps}")
print("Loading first and last timestep")

# Calculate byte offsets
bytes_per_frame = grid_size * 8  # float64
offset_last = (n_timesteps - 1) * bytes_per_frame

# Load first timestep
Nu_bin = data_dir / 'populations_Nu.bin'
Ng_bin = data_dir / 'populations_Ng.bin'

if not Nu_bin.exists() or not Ng_bin.exists():
    print(f"Error: Binary files not found in {data_dir}")
    exit(1)

with open(Nu_bin, 'rb') as f:
    Nu_first = np.fromfile(f, dtype=np.float64, count=grid_size)
with open(Ng_bin, 'rb') as f:
    Ng_first = np.fromfile(f, dtype=np.float64, count=grid_size)

# Load last timestep
with open(Nu_bin, 'rb') as f:
    f.seek(offset_last)
    Nu_last = np.fromfile(f, dtype=np.float64, count=grid_size)
with open(Ng_bin, 'rb') as f:
    f.seek(offset_last)
    Ng_last = np.fromfile(f, dtype=np.float64, count=grid_size)

print("Loaded")

# Reshape and calculate inversion
Nu_first = Nu_first.reshape(NxT, NyT, NzT)
Ng_first = Ng_first.reshape(NxT, NyT, NzT)
Nu_last = Nu_last.reshape(NxT, NyT, NzT)
Ng_last = Ng_last.reshape(NxT, NyT, NzT)

inv_first = Nu_first - Ng_first
inv_last = Nu_last - Ng_last

# Take center slice (z = middle)
z_mid = NzT // 2
slice_first = inv_first[:, :, z_mid]
slice_last = inv_last[:, :, z_mid]
change = slice_last - slice_first

# Stats
print(f"\nTotal inversion change: {np.sum(inv_last) - np.sum(inv_first):.3e}")
print(f"Percent change: {100 * (np.sum(inv_last) - np.sum(inv_first)) / np.sum(inv_first):.1f}%")

# Plot
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

vmax = np.max(np.abs([slice_first, slice_last]))

# Initial
im1 = axes[0].imshow(slice_first.T, cmap='RdBu_r', vmin=-vmax, vmax=vmax, origin='lower')
axes[0].set_title('Initial Inversion (t=0)', fontsize=13, fontweight='bold')
axes[0].set_xlabel('x')
axes[0].set_ylabel('y')
plt.colorbar(im1, ax=axes[0], label='Nu - Ng [m⁻³]')

# Final
im2 = axes[1].imshow(slice_last.T, cmap='RdBu_r', vmin=-vmax, vmax=vmax, origin='lower')
axes[1].set_title(f'Final Inversion (t={n_timesteps - 1})', fontsize=13, fontweight='bold')
axes[1].set_xlabel('x')
axes[1].set_ylabel('y')
plt.colorbar(im2, ax=axes[1], label='Nu - Ng [m⁻³]')

# Change
im3 = axes[2].imshow(change.T, cmap='RdBu_r', origin='lower')
axes[2].set_title('Inversion Change (Depletion)', fontsize=13, fontweight='bold')
axes[2].set_xlabel('x')
axes[2].set_ylabel('y')
plt.colorbar(im3, ax=axes[2], label='Δ(Nu - Ng) [m⁻³]')

plt.tight_layout()
plt.savefig('spatial_gain_depletion.png', dpi=150, bbox_inches='tight')
print("\nSaved: spatial_gain_depletion.png")
plt.close()
