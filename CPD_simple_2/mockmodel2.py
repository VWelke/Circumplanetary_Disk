# Combines run_simple_2(2 dust species)+ rum_simple_3（for Herbig A star)

# Import NumPy for array handling
#
import numpy as np
#
# A simple grid refinement function
#
def grid_refine_inner_edge(x_orig,nlev,nspan):
    x     = x_orig.copy()
    rev   = x[0]>x[1]
    for ilev in range(nlev):
        x_new = 0.5 * ( x[1:nspan+1] + x[:nspan] )
        x_ref = np.hstack((x,x_new))
        x_ref.sort()
        x     = x_ref
        if rev:
            x = x[::-1]
    return x
#
# Some natural constants
#
au  = 1.49598e13     # Astronomical Unit       [cm]
pc  = 3.08572e18     # Parsec                  [cm]
ms  = 1.98892e33     # Solar mass              [g]
ts  = 5.78e3         # Solar temperature       [K]
ls  = 3.8525e33      # Solar luminosity        [erg/s]
rs  = 6.96e10        # Solar radius            [cm]
#
# Monte Carlo parameters
#
nphot    = 1000000
#
# Grid parameters
#
nr       = 100
ntheta   = 32
nphi     = 1
rin      = 0.5*au
rout     = 100*au
thetaup  = np.pi/2-0.7
nlev_rin = 8
nspan_rin= 3
#
# Disk parameters
#
sigmag0  = 1e3               # Sigma gas at 1 AU
sigmad0  = sigmag0 * 0.01    # Sigma dust at 1 AU
fracbig  = 0.99              # Fraction of dust that is the big grain dust
plsig    = -1.0e0            # Powerlaw of the surface density
hr0      = 0.05              # H_p/r at 1 AU
plh      = 0.1               # Powerlaw of flaring
hrbigsett= 0.02              # The big grains are settled a bit more than the small grains
#
# Star parameters
#
mstar    = 1.65*ms
rstar    = 1.6*rs
tstar    = 7650 #K
pstar    = np.array([0.,0.,0.])
#
# Make the coordinates
#
ri       = np.logspace(np.log10(rin),np.log10(rout),nr+1)
ri       = grid_refine_inner_edge(ri,nlev_rin,nspan_rin)
thetai   = np.linspace(thetaup,np.pi/2,ntheta+1)
phii     = np.linspace(0.e0,np.pi*2.e0,nphi+1)
rc       = 0.5 * ( ri[:-1] + ri[1:] )
thetac   = 0.5 * ( thetai[0:ntheta] + thetai[1:ntheta+1] )
phic     = 0.5 * ( phii[0:nphi] + phii[1:nphi+1] )
nr       = len(rc)
#
# Make the grid
#
qq       = np.meshgrid(rc,thetac,phic,indexing='ij')
rr       = qq[0]
tt       = qq[1]
zr       = np.pi/2.e0 - qq[1]
pp       = qq[2]
#
# Make the dust density model
#
sigmad   = sigmad0 * (rr/au)**plsig
sigmadsm = sigmad*(1.-fracbig)
sigmadbg = sigmad*fracbig
hhrsm    = hr0 * (rr/au)**plh
hhrbg    = hrbigsett * (rr/au)**plh
hhsm     = hhrsm * rr
hhbg     = hhrbg * rr
rhodsm   = ( sigmadsm / (np.sqrt(2.e0*np.pi)*hhsm) ) * np.exp(-(zr**2/hhrsm**2)/2.e0)
rhodbg   = ( sigmadbg / (np.sqrt(2.e0*np.pi)*hhbg) ) * np.exp(-(zr**2/hhrbg**2)/2.e0)


# Mask the CPD region 



#CPD parameters
sigmad02 = 0.01*10**3 #(g/cm**2) #g/cm^2 # dust surface density at 1 au
plsig2 = -1.5# power law index for the dust surface density
plh2  = 1.15 # power law index for the dust scale height
hr02 = 0.1
hrbigsett2 =0.05

# Make the dust density model for CPD
#
sigmad2  = sigmad02 * (rr/au)**plsig2
sigmadsm2 = sigmad2*(1.-fracbig)
sigmadbg2 = sigmad2*fracbig
hhrsm2    = hr02 * (rr/au)**plh2
hhrbg2    = hrbigsett2 * (rr/au)**plh2
hhsm2     = hhrsm2 * rr
hhbg2     = hhrbg2 * rr
rhodsm2   = ( sigmadsm2 / (np.sqrt(2.e0*np.pi)*hhsm2) ) * np.exp(-(zr**2/hhrsm2**2)/2.e0)
rhodbg2   = ( sigmadbg2 / (np.sqrt(2.e0*np.pi)*hhbg2) ) * np.exp(-(zr**2/hhrbg2**2)/2.e0)

# Define the radial range for the CPD
r_min = 30 * au  # Inner radius of the CPD region
r_max = 39.7 * au  # Outer radius of the CPD region

# Create a mask for the specified range of rr
mask = (rr >= r_min) & (rr <= r_max)

# Extract the masked radial values
r_masked = rr[mask]

# Define the midpoint of the CPD region
r_mid = 37.2 * au  # Midpoint of the radial range

# Shift the radial array by the midpoint and take the positive values

r_shifted = np.abs(r_masked - r_mid)

# Update the density model within the specified range
# Apply the shifted and positive radial values
sigmad[mask] += sigmad02 * (r_shifted / au) ** plsig2
sigmadsm[mask] = sigmad[mask] * (1. - fracbig)
sigmadbg[mask] = sigmad[mask] * fracbig
hhrsm[mask] += hr02 * (r_shifted / au) ** plh2
hhrbg[mask] += hrbigsett2 * (r_shifted / au) ** plh2
hhsm[mask] = hhrsm[mask] * r_shifted
hhbg[mask] = hhrbg[mask] * r_shifted
rhodsm[mask] = (sigmadsm[mask] / (np.sqrt(2.e0 * np.pi) * hhsm[mask])) * np.exp(-(zr[mask] ** 2 / hhrsm[mask] ** 2) / 2.e0)
rhodbg[mask] = (sigmadbg[mask] / (np.sqrt(2.e0 * np.pi) * hhbg[mask])) * np.exp(-(zr[mask] ** 2 / hhrbg[mask] ** 2) / 2.e0)

# Add a gap between r = 2.2au to r = 26 au 
# Add a gap between r = 2.2au to r = 26 au
#gap_min = 2.2 * au
#gap_max = 26.0 * au
#gap_mask = (rr >= gap_min) & (rr <= gap_max)
#rhodsm[gap_mask] = 0.0
#rhodbg[gap_mask] = 0.0

import matplotlib.pyplot as plt

# Function to plot the small dust density profile
def plot_small_dust_density(rr, rhodsm):
    plt.figure(figsize=(10, 6))
    plt.plot(rr.flatten() / au, rhodsm.flatten(), label='Small Dust Density (rhodsm)')
    plt.xlabel('Radius (au)')
    plt.ylabel('Dust Density')
    plt.title('Small Dust Density Profile')
    plt.legend()
    plt.grid(True)
    plt.show()
plot_small_dust_density(rr, rhodsm)
# Function to plot the big dust density profile
def plot_big_dust_density(rr, rhodbg):
    plt.figure(figsize=(10, 6))
    plt.plot(rr.flatten() / au, rhodbg.flatten(), label='Big Dust Density (rhodbg)')
    plt.xlabel('Radius (au)')
    plt.ylabel('Dust Density')
    plt.title('Big Dust Density Profile')
    plt.legend()
    plt.grid(True)
    plt.show()

#
# Write the wavelength_micron.inp file
#
lam1     = 0.1e0
lam2     = 7.0e0
lam3     = 25.e0
lam4     = 1.0e4
n12      = 20
n23      = 100
n34      = 30
lam12    = np.logspace(np.log10(lam1),np.log10(lam2),n12,endpoint=False)
lam23    = np.logspace(np.log10(lam2),np.log10(lam3),n23,endpoint=False)
lam34    = np.logspace(np.log10(lam3),np.log10(lam4),n34,endpoint=True)
lam      = np.concatenate([lam12,lam23,lam34])
nlam     = lam.size
#
# Write the wavelength file
#
with open('wavelength_micron.inp','w+') as f:
    f.write('%d\n'%(nlam))
    np.savetxt(f,lam.T,fmt=['%13.6e'])
#
#
# Write the stars.inp file
#
with open('stars.inp','w+') as f:
    f.write('2\n')
    f.write('1 %d\n\n'%(nlam))
    f.write('%13.6e %13.6e %13.6e %13.6e %13.6e\n\n'%(rstar,mstar,pstar[0],pstar[1],pstar[2]))
    np.savetxt(f,lam.T,fmt=['%13.6e'])
    f.write('\n%13.6e\n'%(-tstar))
#
# Write the grid file
#
with open('amr_grid.inp','w+') as f:
    f.write('1\n')                       # iformat
    f.write('0\n')                       # AMR grid style  (0=regular grid, no AMR)
    f.write('100\n')                     # Coordinate system: spherical
    f.write('0\n')                       # gridinfo
    f.write('1 1 0\n')                   # Include r,theta coordinates
    f.write('%d %d %d\n'%(nr,ntheta,1))  # Size of grid
    np.savetxt(f,ri.T,fmt=['%21.14e'])    # R coordinates (cell walls)
    np.savetxt(f,thetai.T,fmt=['%21.14e'])# Theta coordinates (cell walls)
    np.savetxt(f,phii.T,fmt=['%21.14e'])  # Phi coordinates (cell walls)
#
# Write the density file
#
with open('dust_density.inp','w+') as f:
    f.write('1\n')                       # Format number
    f.write('%d\n'%(nr*ntheta*nphi))     # Nr of cells
    f.write('2\n')                       # Nr of dust species
    data = rhodsm.ravel(order='F')         # Create a 1-D view, fortran-style indexing
    np.savetxt(f,data.T,fmt=['%13.6e'])  # The data
    data = rhodbg.ravel(order='F')         # Create a 1-D view, fortran-style indexing
    np.savetxt(f,data.T,fmt=['%13.6e'])  # The data
#
# Dust opacity control file
#
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
#
# Write the radmc3d.inp control file
#
with open('radmc3d.inp','w+') as f:
    f.write('nphot = %d\n'%(nphot))
    f.write('scattering_mode_max = 1\n')
    f.write('iranfreqmode = 1\n')