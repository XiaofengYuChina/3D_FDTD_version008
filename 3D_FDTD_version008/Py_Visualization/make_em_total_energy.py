"""
make_em_total_energy.py - Plot total EM energy in simulation box vs time
"""

import matplotlib
matplotlib.use('Agg')  # Non-interactive backend - save only, no display

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from load_grid import find_output_path

# Find output directory (auto-detect debug/release)
base_path = find_output_path(subfolder="EM_total_energy")
print(f"Looking for data in: {base_path}")

if not base_path.exists():
    print(f"Error: Directory not found: {base_path}")
    print("Please run the FDTD simulation first.")
    exit(1)

# Find CSV file
csv_file = base_path / "energy_time.csv"
if not csv_file.exists():
    print(f"Error: energy_time.csv not found in {base_path}")
    exit(1)

print(f"Loading: {csv_file}")

# Read CSV directly, skip the header row
data = np.loadtxt(csv_file, delimiter=",", skiprows=1)

# Split into two columns
t = data[:, 0] * 1e15   # Time array (float) -> fs
U = data[:, 1] * 1e20   # Energy array (float) -> 10^-20 J

print(f"Loaded {len(t)} data points")

# Plot
plt.figure(figsize=(6, 4))
plt.plot(t, U, label="EM energy in box")
plt.xlabel("Time (fs)")
plt.ylabel("Energy (10^-20 J)")
plt.title("Total EM energy in box vs time")
plt.legend()
plt.grid(True)
plt.savefig("EM_total_energy.jpg", dpi=300)
print("Saved: EM_total_energy.jpg")
plt.close()
