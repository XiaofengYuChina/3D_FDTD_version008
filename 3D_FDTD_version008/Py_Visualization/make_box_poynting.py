"""
make_box_poynting.py - Plot refractive index profile (simple version)
"""

import matplotlib
matplotlib.use('Agg')  # Non-interactive backend - save only, no display

import json
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from load_grid import find_output_path

# Find output directory (auto-detect debug/release)
root = find_output_path(subfolder="n_center")
print(f"Looking for data in: {root}")

if not root.exists():
    print(f"Error: Directory not found: {root}")
    print("Please run the FDTD simulation first.")
    exit(1)

# Read metadata
meta_file = root / "metadata.json"
if not meta_file.exists():
    print(f"Error: metadata.json not found in {root}")
    exit(1)

with open(meta_file, "r", encoding="utf-8") as f:
    meta = json.load(f)

NxT = int(meta["NxT"])
NyT = int(meta["NyT"])
NzT = int(meta["NzT"])
kslice = int(meta["kslice"])
npml = int(meta["npml"])
dtype = np.float64 if meta["dtype"] == "float64" else np.float32
pattern = meta["framePattern"]

# Locate the first frame
raw_path = root / (pattern.replace("%04d", "0000"))
if not raw_path.exists():
    # Fallback: in case multiple frames exist, pick the first n_*.raw
    candidates = sorted(root.glob("n_*.raw"))
    if not candidates:
        print("Error: Cannot find n_*.raw files")
        exit(1)
    raw_path = candidates[0]

print(f"Loading: {raw_path}")

# Read binary (row-major, i outer loop, j inner loop)
data = np.fromfile(raw_path, dtype=dtype)
if data.size != NxT * NyT:
    print(f"Error: Size mismatch - loaded {data.size} elements, expected NxT*NyT={NxT * NyT}")
    exit(1)

n_xy = data.reshape((NxT, NyT))  # restore (i, j) layout

# Optional: crop out PML region to view only the physical domain
i0, i1 = npml, NxT - npml
j0, j1 = npml, NyT - npml
n_phys = n_xy[i0:i1, j0:j1]

plt.figure(figsize=(6, 5))
im = plt.imshow(n_phys.T, origin="lower", aspect="equal")
plt.title(f"Refractive index n @ k={kslice}")
plt.xlabel("i (x)")
plt.ylabel("j (y)")
plt.colorbar(im, label="n")
plt.tight_layout()
plt.savefig("refractive_index.jpg", dpi=300)
print("Saved: refractive_index.jpg")
plt.close()
