"""
make_Ez_probe.py - Plot Ez time series at center point
"""

import matplotlib
matplotlib.use('Agg')  # Non-interactive backend - save only, no display

import json
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from load_grid import find_output_path

# Find output directory (auto-detect debug/release)
root = find_output_path(subfolder="Ez_probe")
print(f"Looking for data in: {root}")

if not root.exists():
    print(f"Error: Directory not found: {root}")
    print("Please run the FDTD simulation first.")
    exit(1)

# Read metadata
meta_file = root / "metadata_ts.json"
if not meta_file.exists():
    print(f"Error: metadata_ts.json not found in {root}")
    exit(1)

meta = json.loads(meta_file.read_text())
dt = float(meta["dt"])
saveEvery = int(meta["saveEvery"])

bin_path = root / "ez_center_ts.bin"
if not bin_path.exists():
    print(f"Error: ez_center_ts.bin not found in {root}")
    exit(1)

# Read float64 sequence
data = np.fromfile(bin_path, dtype=np.float64)
print(f"Loaded {len(data)} data points")

# Corresponding time axis (unit: seconds)
t = np.arange(len(data)) * (saveEvery * dt)

plt.figure()
plt.plot(t * 1e15, data)
plt.xlabel("Time (fs)")
plt.ylabel("Ez(center)")
plt.title("Center Ez vs Time")
plt.grid(True)
plt.tight_layout()
plt.savefig("Ez_probe.jpg", dpi=300)
print("Saved: Ez_probe.jpg")
plt.close()
