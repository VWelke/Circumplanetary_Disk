import numpy as np
from astropy.constants import M_sun, L_sun, R_sun, au, pc, G, M_jup, R_jup
from astropy.units import K, g, cm
import matplotlib.pyplot as plt

# Convert constants to cgs units
au = au.cgs.value
pc = pc.cgs.value
M_sun = M_sun.cgs.value
T_sun = 3000
L_sun = L_sun.cgs.value
R_sun = R_sun.cgs.value
M_jup = M_jup.cgs.value
R_jup = R_jup.cgs.value

# Grid refinement function
def grid_refine_mid_plane(theta_orig, nlev, nspan):
    theta = theta_orig.copy()
    rev = theta[0] > theta[1]
    for ilev in range(nlev):
        theta_new = 0.5 * (theta[1:nspan+1] + theta[:nspan])
        theta_ref = np.hstack((theta, theta_new))
        theta_ref.sort()
        theta = theta_ref
        if rev:
            theta = theta[::-1]
    return theta

# Define the parameters of the model
nphot = 1e8
nphot_scat = 1e7
nr = 150
ntheta = 80
nphi = 116
r_in = 0.2 * au
r_out = 90 * au

# Corrected theta grid (full range [0, π])
theta_i = np.linspace(0.0, np.pi, ntheta + 1)
theta_i = grid_refine_mid_plane(theta_i, nlev_thetain=4, nspan_thetain=3)

# Corrected phi grid (strictly closed to 2π)
phi_i = np.linspace(0.0, 2 * np.pi, nphi + 1, endpoint=True)
phi_i[-1] = 2 * np.pi  # Ensure exact closure

# Cell center positions
r_c = 0.5 * (r_i[:-1] + r_i[1:])
theta_c = 0.5 * (theta_i[:-1] + theta_i[1:])
phi_c = 0.5 * (phi_i[:-1] + phi_i[1:])

# Create 3D grid
qq = np.meshgrid(r_c, theta_c, phi_c, indexing='ij')
rr, tt, pp = qq
zr = np.pi / 2.0 - tt  # z = π/2 - θ

# Density model
ndustspec = 2
sigmag0 = 1e3
sigmad0 = sigmag0 * 0.01
fracbig = 0.99
plsig = -1.0
hr0 = 0.05
plh = 0.1
hrbigsett = 0.02

sigmad = sigmad0 * (rr / au) ** plsig
sigmadsm = sigmad * (1.0 - fracbig)
sigmadbg = sigmad * fracbig
hhrsm = hr0 * (rr / au) ** plh
hhrbg = hrbigsett * (rr / au) ** plh
hhsm = hhrsm * rr
hhbg = hhrbg * rr
rhodsm = (sigmadsm / (np.sqrt(2.0 * np.pi) * hhsm)) * np.exp(-(zr**2 / hhrsm**2) / 2.0)
rhodbg = (sigmadbg / (np.sqrt(2.0 * np.pi) * hhbg)) * np.exp(-(zr**2 / hhrbg**2) / 2.0)

# CPD parameters
sigmad02 = 2127
plsig2 = -1.2
plh2 = 1.15
hr02 = 0.1
hrbigsett2 = 0.05

# CPD region mask
r_min = 36.2 * au
r_max = 38.2 * au
phi_width = np.deg2rad(5.0)  # CPD azimuthal width (adjust as needed)
phi_min = 2 * np.pi - phi_width
phi_max = phi_width
mask_r = (rr >= r_min) & (rr <= r_max)
mask_phi = (pp >= phi_min) | (pp <= phi_max)
mask = mask_r & mask_phi

# Shifted radial coordinates for CPD
r_masked = rr[mask]
pp_masked = pp[mask]
r_mid = 37.2 * au
r_shifted = np.sqrt(r_masked**2 + r_mid**2 - 2 * r_masked * r_mid * np.cos(pp_masked))

# Update density model with CPD
sigmad[mask] += sigmad02 * (r_shifted / au) ** plsig2
sigmadsm[mask] = sigmad[mask] * (1.0 - fracbig)
sigmadbg[mask] = sigmad[mask] * fracbig
hhrsm[mask] += hr02 * (r_shifted / au) ** plh2
hhrbg[mask] += hrbigsett2 * (r_shifted / au) ** plh2
hhsm[mask] = hhrsm[mask] * r_shifted
hhbg[mask] = hhrbg[mask] * r_shifted
rhodsm[mask] = (sigmadsm[mask] / (np.sqrt(2.0 * np.pi) * hhsm[mask])) * np.exp(-(zr[mask]**2 / hhrsm[mask]**2) / 2.0)
rhodbg[mask] = (sigmadbg[mask] / (np.sqrt(2.0 * np.pi) * hhbg[mask])) * np.exp(-(zr[mask]**2 / hhrbg[mask]**2) / 2.0)

# Write input files
with open('radmc3d.inp', 'w+') as f:
    f.write('nphot = %d\n' % nphot)
    f.write('nphot_scat = %d\n' % nphot_scat)
    f.write('scattering_mode_max = 1\n')
    f.write('iranfreqmode = 1\n')
    f.write('modified_random_walk = 1\n')
    f.write('phibound_rotational = 1\n')  # Enable periodic boundary in phi

with open('amr_grid.inp', 'w+') as f:
    f.write('1\n')
    f.write('0\n')
    f.write('111\n')
    f.write('0\n')
    f.write('1 1 1\n')
    f.write('%d %d %d\n' % (nr, ntheta, nphi))
    for value in r_i:
        f.write('%13.6e\n' % value)
    for value in theta_i:
        f.write('%17.10e\n' % value)
    for value in phi_i:
        f.write('%13.6e\n' % value)

with open('dust_density.inp', 'w+') as f:
    f.write('1\n')
    f.write('%d\n' % (nr * ntheta * nphi))
    f.write('2\n')
    data = rhodsm.ravel(order='F')
    np.savetxt(f, data.T, fmt=['%13.6e'])
    data = rhodbg.ravel(order='F')
    np.savetxt(f, data.T, fmt=['%13.6e'])

# Star and planet parameters
mstar = 1.65 * M_sun
rstar = 1.6 * R_sun
tstar = 4000
pstar = np.array([0, 0, 0])

mplanet = 3 * M_jup
rplanet = 1.17 * R_jup
tplanet = 1200
pplanet = np.array([37.2 * au, 0, 0])

# Wavelengths
lam = np.logspace(np.log10(0.1), np.log10(1e4), 150)
nlam = lam.size

# Write stars.inp
with open('stars.inp', 'w+') as f:
    f.write('2\n')
    f.write('2 %d\n\n' % nlam)
    f.write('%13.6e %13.6e %13.6e %13.6e %13.6e\n\n' % (rstar, mstar, pstar[0], pstar[1], pstar[2]))
    f.write('%13.6e %13.6e %13.6e %13.6e %13.6e\n\n' % (rplanet, mplanet, pplanet[0], pplanet[1], pplanet[2]))
    for value in lam:
        f.write('%13.6e\n' % value)
    f.write('\n%13.6e\n' % (-tstar))
    f.write('\n%13.6e\n' % (-tplanet))

# Write dustopac.inp
with open('dustopac.inp', 'w+') as f:
    f.write('2\n')
    f.write('2\n')
    f.write('============================================================================\n')
    f.write('1\n')
    f.write('0\n')
    f.write('0.1_micron\n')
    f.write('----------------------------------------------------------------------------\n')
    f.write('1\n')
    f.write('0\n')
    f.write('100_micron\n')
    f.write('----------------------------------------------------------------------------\n')

# Write wavelength_micron.inp
with open('wavelength_micron.inp', 'w+') as f:
    f.write('%d\n' % nlam)
    for value in lam:
        f.write('%13.6e\n' % value)

# Plot density profile
plt.figure(figsize=(10, 6))
plt.plot(rr.flatten() / au, np.log10(sigmad.flatten()), label='Dust Density')
plt.xlabel('Radius (au)')
plt.ylabel('log(Dust Density)')
plt.title('Dust Density Profile')
plt.legend()
plt.grid(True)
plt.show()