
#!/usr/bin/env python3

# Script to modify dust_density.inp for CPD grid cells and visualize results

import numpy as np
import matplotlib.pyplot as plt

# Grid parameters (must match your original model)
nr = 150        # Number of radial grid cells
ntheta = 80     # Number of theta grid cells (after refinement)
nphi = 116      # Number of phi grid cells (updated from your script)

# Define radial grid (assuming same as original script)
au = 1.496e13   # AU in cm
r_in = 0.2 * au
r_out = 90 * au
r_i = np.linspace(r_in, r_out, nr + 1)  # Cell walls
r_c = 0.5 * (r_i[:-1] + r_i[1:])        # Cell centers

# Read the existing dust_density.inp file
with open('dust_density.inp', 'r') as f:
    lines = f.readlines()

# Verify file structure
n_cells = int(lines[1].strip())  # Total number of cells
n_dust_species = int(lines[2].strip())  # Number of dust species
expected_cells = nr * ntheta * nphi
if n_cells != expected_cells:
    raise ValueError(f"Expected {expected_cells} cells, but found {n_cells} in dust_density.inp")

# Extract density arrays from the file
dust1 = np.array([float(line.strip()) for line in lines[3:3 + n_cells]])           # Small grains
dust2 = np.array([float(line.strip()) for line in lines[3 + n_cells:3 + 2 * n_cells]])  # Large grains

# CPD grid cells: r = 37.2 au, phi ≈ 0, and a few z values near midplane
cpd_rmin = 36.2 * au  # 36.2 AU in cm
cpd_rmax = 38.2 * au  # 38.2 AU in cm

# Radial indices for CPD (ir_cpd = [59, 60, 61, 62, 63] from your script)
ir_cpd = [59, 60, 61, 62, 63]  # Corresponds to r_c ≈ 36.2 to 38.2 AU
# Verify: r_c[59] ≈ 36.1 au, r_c[63] ≈ 38.5 au

# CPD azimuthal range (near phi = 0 and phi = 2π)
iphi_cpd = [0, 1, 114, 115]  # Adjusted for nphi = 116

# CPD vertical range (near midplane)
itheta_cpd = range(77, 80)  # Near midplane (z ≈ 0)

# Collect indices of all CPD cells
cpd_indices = []
for iphi in iphi_cpd:
    for itheta in itheta_cpd:
        for ir in ir_cpd:
            idx = ir + itheta * nr + iphi * nr * ntheta  # Fortran-style index
            if idx < n_cells:  # Ensure index is within bounds
                cpd_indices.append(idx)
            else:
                print(f"Warning: Index {idx} exceeds grid size {n_cells-1}")

# Print number of CPD cells
print(f"Number of CPD cells: {len(cpd_indices)}")

# Print original density values for the selected cells
print("Original density values for CPD cells:")
print("Index | Small grains (dust1) | Large grains (dust2)")
print("---------------------------------------------------")
for idx in cpd_indices:
    print(f"{idx:5d} | {dust1[idx]:18.6e} | {dust2[idx]:18.6e}")

# Modify density values (increase by factor of 100)
mag_factor = 1
for idx in cpd_indices:
    dust1[idx] *= mag_factor  # Small grains
    dust2[idx] *= mag_factor  # Large grains

print(f"Density values modified by a factor of {mag_factor}")
# Update lines with modified densities
for idx in cpd_indices:
    lines[3 + idx] = f'{dust1[idx]:13.6e}\n'              # Small grains
    lines[3 + n_cells + idx] = f'{dust2[idx]:13.6e}\n'    # Large grains

# Write the modified data back to dust_density.inp
with open('dust_density.inp', 'w') as f:
    f.writelines(lines)

# Report the number of modified cells
print(f"Modified {len(cpd_indices)} cells for CPD at r = 37.2 au, phi ≈ 0.")

############################################################################################################
