
import numpy as np


# Import some astronomical constants
from astropy.constants import M_sun, L_sun, R_sun, au, pc, G, M_jup, R_jup
from astropy.units import K, g, cm


# Convert constants to cgs units
au = au.cgs.value     # Astronomical Unit       [cm]
pc = pc.cgs.value     # Parsec                  [cm]
M_sun = M_sun.cgs.value  # Solar mass              [g]
T_sun = 3000      # Solar temperature       [K]
L_sun = L_sun.cgs.value  # Solar luminosity        [erg/s]
R_sun = R_sun.cgs.value  # Solar radius            [cm]
M_jup = M_jup.cgs.value  # Jupiter mass            [g]
R_jup = R_jup.cgs.value  # Jupiter radius          [cm]   


# Add grid refinement before making the coordinates
# add additional points between existing points
# nspan: The number of points from the start of the array to consider for refinement.
# midplane has mirror symmetry, so only need to refine one side
def grid_refine_mid_plane(theta_orig, nlev, nspan):
    theta = theta_orig.copy()
    rev = theta[0] > theta[1]
    for ilev in range(nlev):
        theta_new = 0.5 * (theta[1:nspan+1] + theta[:nspan])   # #refine disk surface
        theta_ref = np.hstack((theta, theta_new))
        theta_ref.sort()
        theta = theta_ref
        if rev:
            theta = theta[::-1]
    return theta

# if want to find gird_refinement in phi, go back to model 8


nlev_thetain = 4 
nspan_thetain= 3  


# Define the parameters of the model

    # radmc3d.inp parameter : main settings for RADMC-3D
nphot    = 1e8  #for the thermal monte carto simulation
nphot_scat = 1e8  #for the scattering monte carto simulation
#multiple CPU cores, may need cluster

    # Grid : defines layout of space

# number go back to CPD_simple_1_test_no_planet_grid_refinement
nr       = 150 
ntheta   = 80
nphi     = 200


# Radius for PPD not CPD
# r:  inner CPD (0.2 au to 2.2 au) , gap (2.2 au to 26 au) , outer dust ring (26 au to 90 au)
# (Poblette et al 2022), but doesnt include the inner CPD rim, prob need calculagte truncation radius from star mass
# CPD locates 37.2 au away from star (r = 37.2 )
r_in      = 0.2*au   # 3 times Jupyter radius
r_out     = 90*au 

theta_up  = np.pi*0.5 - 0.7e0  
        # Coordinate array
r_i       = np.linspace(r_in,r_out,nr+1)
#r_i       = np.logspace(np.log10(r_in),np.log10(r_out),nr+1)  #+1 because it is not counting cell centers, but the walls
#r_i       = grid_refine_inner_edge(r_i,nlev_rin,nspan_rin)

theta_i   = np.linspace(theta_up,0.5e0*np.pi,ntheta+1, endpoint =False)  # theta goes to pi/2 lets z starts from zero 
theta_i   = np.append(theta_i, 0.5e0*np.pi)  # Now exactly periodic
#theta_i   = grid_refine_mid_plane(theta_i, nlev_thetain, nspan_thetain)


print(np.pi/2.e0-theta_i)


phi_i = np.linspace(0., 2.*np.pi, nphi+1, endpoint=False)
phi_i = np.append(phi_i, 2.*np.pi)  # Now exactly periodic

        # Cell center position array

r_c       = 0.5 * ( r_i[:-1] + r_i[1:] )
theta_c   = 0.5 * ( theta_i[:-1] + theta_i[1:] )
phi_c     = 0.5 * ( phi_i[:-1] + phi_i[1:] )

ntheta       = len(theta_c)  
nphi         = len(phi_c)

print(phi_i)
print(phi_c)
print( f'the number of phi grid cells are {nphi}')
print(f'The number of theta grid cells are {ntheta}')
        # Make the grid
            # takes in the center positions of the cells and returns a 3D matrix of the grid
            # indexing='ij' means that the first two indices of the 3D matrix are the r and theta coordinates
qq       = np.meshgrid(r_c,theta_c,phi_c,indexing='ij')
            # Extract the coordinates (r,theta, z) from the 3D matrix
rr       = qq[0]  # final r coor defined by cell center
tt       = qq[1] # final theta coor defined by cell center, just for defining zr
pp       = qq[2]
zr       = np.pi/2.e0 - qq[1]    # z = pi/2 - theta, essentially frame rotated by 90 degrees, and z is from 0 to 0.7 radians



    # Density: dust values (for each species) for each cell in the grid
                                                                                                                                                             
        # number of dust species
ndustspec = 2


sigmag0  = 1e3               # Sigma gas at 1 AU
sigmad0  = sigmag0 * 0.01    # Sigma dust at 1 AU
fracbig  = 0.99              # Fraction of dust that is the big grain dust
plsig    = -1.0e0            # Powerlaw of the surface density
hr0      = 0.05              # H_p/r at 1 AU
plh      = 0.1               # Powerlaw of flaring
hrbigsett= 0.02              # The big grains are settled a bit more than the small grains
        # dust density function

sigmad   = sigmad0 * (rr/au)**plsig
sigmadsm = sigmad*(1.-fracbig)
sigmadbg = sigmad*fracbig
hhrsm    = hr0 * (rr/au)**plh
hhrbg    = hrbigsett * (rr/au)**plh
hhsm     = hhrsm * rr
hhbg     = hhrbg * rr
rhodsm   = ( sigmadsm / (np.sqrt(2.e0*np.pi)*hhsm) ) * np.exp(-(zr**2/hhrsm**2)/2.e0)
rhodbg   = ( sigmadbg / (np.sqrt(2.e0*np.pi)*hhbg) ) * np.exp(-(zr**2/hhrbg**2)/2.e0)



#CPD parameters
#sigma_g0 =  10**3 #(g/cm**2)   # gas surface density at 1 au
sigmad02 = 0.016 #(g/cm**2) #g/cm^2 # dust surface density at 1 au
plsig2 = -1.2# power law index for the dust surface density
plh2  = 1.15 # power law index for the dust scale height
hr02 = 0.1
hrbigsett2 =0.05

# Make the dust density model for CPD


# Define the radial range for the CPD
r_min = 34.7 * au  # Inner radius of the CPD region
r_max = 39.7 * au  # Outer radius of the CPD region
phi_max = np.arctan(2.5/37.2) # Azimuthal range of the CPD region


# Create a mask for the specified range of rr
mask_r = (rr >= r_min) & (rr <= r_max)
mask_phi = (pp >= np.pi/2 - phi_max) & (pp <= np.pi/2 + phi_max)
print(f'The number of phi grid cells in the mask are {np.sum(mask_phi)}')
print(f'The number of r grid cells in the mask are {np.sum(mask_r)}')
mask = mask_r & mask_phi
print(f'The number of cells in the CPD region is {np.sum(mask)}')
# Extract the masked radial values
r_masked = rr[mask]
pp_masked = pp[mask]

# Define the midpoint of the CPD region
r_mid = 37.2*au # Midpoint of the radial range

# Shift the radial array by the midpoint and take the positive values

r_shifted = np.sqrt(r_masked**2 + r_mid**2 - 2*r_masked*r_mid*np.cos(pp_masked-np.pi/2))
print(np.min(r_shifted/au), np.max(r_shifted/au))

# Update the density model within the specified range
# Apply the shifted and positive radial values
# ongoing discussion about shape (Sun etal 2024)
sigmad[mask] += sigmad02 * (r_shifted / au) ** plsig2
sigmadsm[mask] = sigmad[mask] * (1. - fracbig)
sigmadbg[mask] = sigmad[mask] * fracbig
hhrsm[mask] += hr02 * (r_shifted / au) ** plh2
hhrbg[mask] += hrbigsett2 * (r_shifted / au) ** plh2
hhsm[mask] = hhrsm[mask] * r_shifted
hhbg[mask] = hhrbg[mask] * r_shifted
rhodsm[mask] = (sigmadsm[mask] / (np.sqrt(2.e0 * np.pi) * hhsm[mask])) * np.exp(-(zr[mask] ** 2 / hhrsm[mask] ** 2) / 2.e0)
rhodbg[mask] = (sigmadbg[mask] / (np.sqrt(2.e0 * np.pi) * hhbg[mask])) * np.exp(-(zr[mask] ** 2 / hhrbg[mask] ** 2) / 2.e0)

# add azimuthal averaging
azav  = False                  # Switch on the azimuthal averaging

if azav:
    rho_dust_3d_av = []
    for rho in rhodsm,rhodbg:
        rhoav = rho.mean(axis=2)
        rho_dust_3d_av.append(rhoav)
    rhodsm,rhodbg = rho_dust_3d_av
    phi_i = [phi_i[0],phi_i[-1]]
    nphi = 1








import matplotlib.pyplot as plt
def plot_small_dust_density_CPD(r_shifted, rhodsm):
    plt.figure(figsize=(10, 6))
    plt.plot(r_shifted.flatten(), np.log10(sigmad[mask]), label='Small Dust Density (rhodsm)')
    plt.xlabel('Radius (au)')
    plt.ylabel('Dust Density')
    plt.title('Small Dust Density Profile')
    plt.legend()
    plt.grid(True)
    plt.show()
#plot_small_dust_density_CPD(r_shifted, rhodsm)

def plot_small_dust_density(rr, rhodsm):
    plt.figure(figsize=(10, 6))
    plt.plot(rr.flatten() / au, np.log10(rhodsm.flatten()), label='Small Dust Density (rhodsm)')
    plt.xlabel('Radius (au)')
    plt.ylabel('Dust Density')
    plt.title('Small Dust Density Profile')
    plt.legend()
    plt.grid(True)
    plt.show()

#plot_small_dust_density(rr, rhodsm)
    # Star and planet parameters  
        # Star parameters

mstar    = 1.65*M_sun  #1.65 Msun
rstar    = 1.6*R_sun  #Rsun
tstar    = 4000 #7650 #K
pstar    = np.array([0,0.,0.])  # 37.2 au at R later


        # Planet parameters
mplanet  = 3*M_jup
rplanet  = 1.17*R_jup
tplanet  = 1200



# Place the planet at the cell center
ir = np.argmin(np.abs(r_c - 37.2*au))
iphi = np.argmin(np.abs(phi_c - np.pi/2))
print(f'ir = {ir}, iphi = {iphi}')
pplanet = [r_c[ir], 0.0, phi_c[iphi]]  # Cell center coordinates
print(f'pplanet = {pplanet}')

    # Wavelengths
    # ALMA wavelength range for different bands in microns
    # https://www.eso.org/public/teles-instr/alma/receiver-bands/
lam1     = 0.1e0
lam2     = 25e0
lam3     = 25e1
lam4     = 1.0e4
n12      = 20
n23      = 30
n34      = 100
lam12    = np.logspace(np.log10(lam1),np.log10(lam2),n12,endpoint=False)
lam23    = np.logspace(np.log10(lam2),np.log10(lam3),n23,endpoint=False)
lam34    = np.logspace(np.log10(lam3),np.log10(lam4),n34,endpoint=True)
lam      = np.concatenate([lam12,lam23,lam34])
nlam     = lam.size

# Check grid periodicity
print("Phi grid exact 2pi?:", phi_i[-1] == 2*np.pi)  # Should be True
print("Theta exact pi/2?:", theta_i[-1] == np.pi/2)  # Should be True

# Write to .inp files

# Main setting
with open('radmc3d.inp', 'w+') as f:  # Open the file in write-plus (can both read and write the file) mode
    f.write('nphot = %d\n' % (nphot))  # means nphot = <value of nphot>, %d is placeholder later replaced by value of nphot, \n is new line 
    f.write('nphot_scat = %d\n' % (nphot_scat))  
    f.write('scattering_mode_max = 1\n')  # if scattering opacity is included anywhere
    f.write('iranfreqmode = 1\n')  # frequency grid is based on input files
    f.write('modified_random_walk =1\n') #for optically thick

# Grid
with open('amr_grid.inp','w+') as f:
    f.write('1\n')                       # iformat
    f.write('0\n')                       # AMR grid style  (0=regular grid, no AMR)
    f.write('111\n')                     # Coordinate system: spherical, let all axes be active
    f.write('0\n')                       # gridinfo
    if nphi>1:
        f.write('1 1 1\n')                  # Include r,theta,phi coordinates
    else:
        f.write('1 1 0\n')                  # Include r,theta coordinates
    f.write('%d %d %d\n'%(nr,ntheta,nphi))  # Size of grid
    for value in r_i:
        f.write('%13.6e\n'%(value))      # X coordinates (cell walls)
    for value in theta_i:
        f.write('%17.10e\n'%(value))     # Y coordinates (cell walls) (use higher precision here)
    for value in phi_i:
        f.write('%13.6e\n'%(value))      # Z coordinates (cell walls)
#

with open('dust_density.inp','w+') as f:
    f.write('1\n')                       # Format number
    f.write('%d\n'%(nr*ntheta*nphi))     # Nr of cells
    f.write('2\n')                       # Nr of dust species
    data = rhodsm.ravel(order='F')         # Create a 1-D view, fortran-style indexing
    np.savetxt(f,data.T,fmt=['%13.6e'])  # The data
    data = rhodbg.ravel(order='F')         # Create a 1-D view, fortran-style indexing
    np.savetxt(f,data.T,fmt=['%13.6e'])  # The data

# Stars and Planets
with open('stars.inp', 'w+') as f:
    f.write('2\n')  # 2: lambda in microns, 1: freq in hz
    f.write('2 %d\n\n' % (nlam))  # number of stars, number of wavelengths
    f.write('%13.6e %13.6e %13.6e %13.6e %13.6e\n\n'%(rstar, mstar, pstar[0], pstar[1], pstar[2]))  # star's r, M, pos(x, y, z)
    f.write('%13.6e %13.6e %13.6e %13.6e %13.6e\n\n' % (rplanet, mplanet, pplanet[0], pplanet[1], pplanet[2]))  # planet's r, M, pos(x, y, z)
    for value in lam:
        f.write('%13.6e\n'%(value))
    f.write('\n%13.6e\n'%(-tstar))  # write negative tstar
    f.write('\n%13.6e\n' % (-tplanet))  # write negative tplanet



# Dust opacity
with open('dustopac.inp','w+') as f:
    f.write('2               Format number of this file\n')
    f.write('2               Nr of dust species\n')
    f.write('============================================================================\n')
    f.write('1               Way in which this dust species is read\n')
    f.write('0               0=Thermal grain\n')
    f.write('0.1_micron      Extension of name of dustkappa_***.inp file\n')
    f.write('----------------------------------------------------------------------------\n')
    f.write('1               Way in which this dust species is read\n')
    f.write('0               0=Thermal grain\n')
    f.write('100_micron      Extension of name of dustkappa_***.inp file\n')
    f.write('----------------------------------------------------------------------------\n')

with open('wavelength_micron.inp','w+') as f:
    f.write('%d\n'%(nlam))
    for value in lam:
        f.write('%13.6e\n'%(value))



#camera_wavelength_micron.inp

# wavelength for make image



lam_cam = np.array([7500, 3000, 2000, 1621.62, 1304.35, 869.57, 652.17, 461.54, 344.83])

nlam_cam = len(lam_cam)



with open('camera_wavelength_micron.inp','w+') as f:
    #f.write('1\n')                       # Format number
    f.write('%d\n'%(nlam_cam))
    for value in lam_cam:
        f.write('%13.6e\n'%(value))


# Plot the 2D slice of sigmad
rr_slice = rr[:, 0, 0]
pp_slice = pp[:, 0, :]
sigmad_slice = sigmad[:, 0, :]

plt.figure(figsize=(10, 6))
plt.plot(np.log10(rr_slice / au), np.log10(sigmad_slice[:, 50]), label='phi = 0')
#plt.plot(rr_slice / au, sigmad_slice[:, -1], label='phi = 2pi')
#plt.plot(rr_slice / au, sigmad_slice[:, iphi], label='phi = pi')
plt.legend()
#plt.xlim(np.log10(r_min / au), np.log10(r_max / au))
#plt.xlim(r_min / au, r_max / au)
#plt.ylim(0,0.5)
plt.xlabel('log(r(au))')
plt.ylabel('log(sigmad) (g/cm^3)')
plt.title('2D Plot of sigmad')
plt.grid(True)
plt.show()

