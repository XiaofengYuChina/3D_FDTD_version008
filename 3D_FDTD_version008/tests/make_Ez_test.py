"""
make_Ez_test.py - Enhanced Ez field visualization with logarithmic color scaling
"""

import matplotlib
matplotlib.use('Agg')  # Non-interactive backend - save only, no display

import json
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter
from matplotlib.colors import SymLogNorm, LinearSegmentedColormap
import imageio_ffmpeg
from load_grid import find_output_path

matplotlib.rcParams['animation.ffmpeg_path'] = imageio_ffmpeg.get_ffmpeg_exe()

# Find output directory (auto-detect debug/release)
data_path = find_output_path(subfolder="Ez_center")
print(f"Looking for data in: {data_path}")

if not data_path.exists():
    print(f"Error: Directory not found: {data_path}")
    print("Please run the FDTD simulation first.")
    exit(1)

# Read metadata
meta_file = data_path / "metadata.json"
if not meta_file.exists():
    print(f"Error: metadata.json not found in {data_path}")
    exit(1)

with open(meta_file, "r") as f:
    meta = json.load(f)

NxT, NyT = int(meta["NxT"]), int(meta["NyT"])
Nx, Ny, npml = meta["Nx"], meta["Ny"], meta["npml"]

frame_files = sorted(glob.glob(str(data_path / "ez_*.raw")))
print(f"Found {len(frame_files)} frames")

if len(frame_files) == 0:
    print("Error: No frame files found!")
    print(f"Expected files like: {data_path}/ez_0000.raw")
    exit(1)

# Find global max
maxabs = 0.0
for fp in frame_files:
    z = np.fromfile(fp, dtype=np.float64).reshape(NxT, NyT).T
    m = float(np.max(np.abs(z)))
    if m > maxabs:
        maxabs = m
if maxabs == 0.0:
    maxabs = 1.0

print(f"Global max|Ez| = {maxabs:.2e}")

# Colormap selection
colors = [
    '#000080',  # Navy (strong negative)
    '#0000CD',  # Medium blue
    '#1E90FF',  # Dodger blue
    '#00008B',  # Dark blue
    '#4B0082',  # Indigo
    '#2F4F4F',  # Dark slate gray (near zero)
    '#8B4513',  # Saddle brown
    '#B8860B',  # Dark goldenrod
    '#FF8C00',  # Dark orange
    '#FF4500',  # Orange red
    '#8B0000',  # Dark red (strong positive)
]
cmap = LinearSegmentedColormap.from_list('visible_diverging', colors, N=512)

# Logarithmic color scaling
norm = SymLogNorm(linthresh=maxabs / 1000, vmin=-maxabs, vmax=maxabs, base=10)

# Create figure
fig, ax = plt.subplots(figsize=(12, 10), facecolor='#2b2b2b')
ax.set_facecolor('#1a1a1a')

im = ax.imshow(np.zeros((NyT, NxT)), cmap=cmap, norm=norm,
               origin="upper", aspect="auto")
cbar = plt.colorbar(im, ax=ax, extend='both')
cbar.set_label("Ez (arb. units, log scale)", fontsize=12, color='white')
cbar.ax.tick_params(colors='white')  # White colorbar ticks

ax.set_xlabel("i (cell index)", fontsize=12, color='white')
ax.set_ylabel("j (cell index)", fontsize=12, color='white')
ax.tick_params(colors='white')  # White axis ticks

# PML boundaries
for x in [npml - 0.5, npml + Nx - 0.5]:
    ax.axvline(x, color='lime', linestyle='--', linewidth=1.5, alpha=0.9)
for y in [npml - 0.5, npml + Ny - 0.5]:
    ax.axhline(y, color='lime', linestyle='--', linewidth=1.5, alpha=0.9)

# Text box for per-frame max
text_box = ax.text(0.02, 0.98, '', transform=ax.transAxes,
                   fontsize=11, verticalalignment='top', color='white',
                   bbox=dict(boxstyle='round', facecolor='#404040', alpha=0.9, edgecolor='lime'))

# Generate movie
print("Generating MP4...")
try:
    writer = FFMpegWriter(fps=20, codec="libx264", bitrate=6000)
    with writer.saving(fig, "ez_evolution_enhanced.mp4", dpi=150):
        for i, fp in enumerate(frame_files):
            z = np.fromfile(fp, dtype=np.float64).reshape(NxT, NyT).T
            frame_max = np.max(np.abs(z))

            im.set_data(z)
            ax.set_title(f"Ez Field @ z=Nz/2 (Frame {i}/{len(frame_files)})",
                         fontsize=14, color='white', pad=15)
            text_box.set_text(f'Frame max: {frame_max:.2e}')

            writer.grab_frame()
            if (i + 1) % 50 == 0:
                print(f"MP4 progress: {i + 1}/{len(frame_files)}")
    print("Saved: ez_evolution_enhanced.mp4")
except Exception as e:
    print(f"Error: {e}")
finally:
    plt.close(fig)

print("Done!")
