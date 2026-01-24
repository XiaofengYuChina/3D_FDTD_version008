"""
make_Ez_movie.py - Generate Ez field evolution movie (MP4 and GIF)
"""

import matplotlib
matplotlib.use('Agg')  # Non-interactive backend - save only, no display

import json
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter
import imageio.v3 as iio
import imageio_ffmpeg
from load_grid import find_output_path

# Tell matplotlib where to find ffmpeg
matplotlib.rcParams['animation.ffmpeg_path'] = imageio_ffmpeg.get_ffmpeg_exe()
print(f"ffmpeg: {imageio_ffmpeg.get_ffmpeg_exe()}")

# Find output directory (auto-detect debug/release)
data_path = find_output_path(subfolder="Ez_center")
print(f"Looking for data in: {data_path}")

if not data_path.exists():
    print(f"Error: Directory not found: {data_path}")
    print("Please run the FDTD simulation first.")
    exit(1)

# 1) Read metadata
meta_file = data_path / "metadata.json"
if not meta_file.exists():
    print(f"Error: metadata.json not found in {data_path}")
    exit(1)

with open(meta_file, "r", encoding="utf-8") as f:
    meta = json.load(f)

NxT, NyT = int(meta["NxT"]), int(meta["NyT"])
Nx, Ny, npml = meta["Nx"], meta["Ny"], meta["npml"]

frame_files = sorted(glob.glob(str(data_path / "ez_*.raw")))
print(f"Found {len(frame_files)} frame files")

# Check if we have any frames
if len(frame_files) == 0:
    print("Error: No frame files found!")
    print(f"Expected files like: {data_path}/ez_0000.raw")
    exit(1)

# 2) Scan all frames to determine global color scale
maxabs = 0.0
for fp in frame_files:
    z = np.fromfile(fp, dtype=np.float64).reshape(NxT, NyT).T
    m = float(np.max(np.abs(z)))
    if m > maxabs:
        maxabs = m
if maxabs == 0.0:
    maxabs = 1.0

print(f"Global max |Ez| = {maxabs:.2e}")

# 3) Generate MP4
print("Generating MP4...")
fig, ax = plt.subplots()
im = ax.imshow(np.zeros((NyT, NxT)), cmap="seismic", vmin=-100, vmax=+100, origin="upper")
cbar = plt.colorbar(im, ax=ax)
cbar.set_label("Ez (arb. units)")
ax.set_title("Ez @ z = Nz/2")
ax.set_xlabel("i (cell index)")
ax.set_ylabel("j (cell index)")

# Physical region boundaries
xlines = [npml - 0.5, npml + Nx - 0.5]
ylines = [npml - 0.5, npml + Ny - 0.5]
for x in xlines:
    ax.axvline(x, color="k", linestyle="--", linewidth=1)
for y in ylines:
    ax.axhline(y, color="k", linestyle="--", linewidth=1)

# Try to generate MP4 (requires ffmpeg)
try:
    writer = FFMpegWriter(fps=20, codec="libx264", bitrate=4000)
    with writer.saving(fig, "ez_evolution.mp4", dpi=140):
        for i, fp in enumerate(frame_files):
            z = np.fromfile(fp, dtype=np.float64).reshape(NxT, NyT).T
            im.set_data(z)
            writer.grab_frame()
            if (i + 1) % 50 == 0:
                print(f"MP4 progress: {i + 1}/{len(frame_files)}")
    print("Saved: ez_evolution.mp4")
except FileNotFoundError:
    print("Warning: ffmpeg not found. MP4 generation skipped.")
    print("Install: pip install imageio-ffmpeg")
except Exception as e:
    print(f"Warning: MP4 generation failed: {e}")
finally:
    plt.close(fig)

# 4) Optional: Save GIF
if len(frame_files) > 0:
    print("Generating GIF...")
    imgs = []
    for fp in frame_files:
        z = np.fromfile(fp, dtype=np.float64).reshape(NxT, NyT).T
        g = (z / maxabs) * 0.5 + 0.5  # normalize to 0..1
        g = np.clip(g, 0.0, 1.0)
        imgs.append((g * 255).astype(np.uint8))
    iio.imwrite("ez_evolution.gif", imgs, duration=0.05)
    print("Saved: ez_evolution.gif")
else:
    print("Skipping GIF: no frames available")

print("Done!")
