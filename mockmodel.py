import numpy as np
# Reference files from radmc3d-2.0/examples/run_ppdisk_simple_1/problem_setup.py

# Make some plots in between to visualize the grid and density distribution

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





# Define the parameters of the model

    # radmc3d.inp parameter : main settings for RADMC-3D
nphot    = 1000000  #for the thermal monte carto simulation


    # Grid : defines layout of space

nr       = 32 
ntheta   = 32
nphi     = 1 #axisymmetric for each r and theta 
r_in      = 0.004779*au   # disk inner radius
r_out     = 0.93*au # disk outer radius 
theta_up  = np.pi*0.5 - 0.7e0  # mighr need to be adjusted

        # Coordinate array

r_i       = np.logspace(np.log10(r_in),np.log10(r_out),nr+1)  #+1 because it is not counting cell centers, but the walls
theta_i   = np.linspace(theta_up,0.5e0*np.pi,ntheta+1)  # theta goes to pi/2 lets z starts from zero 
phi_i     = np.linspace(0.e0,np.pi*2.e0,nphi+1)

        # Cell center position array
r_c = 0.5*(r_i[0:nr]+r_i[1:nr+1])  # average of the cell wall x-pos and cell wall x+1 pos for each cell in the array
theta_c = 0.5*(theta_i[0:ntheta]+theta_i[1:ntheta+1])
phi_c = 0.5*(phi_i[0:nphi]+phi_i[1:nphi+1])


        # Make the grid
            # takes in the center positions of the cells and returns a 3D matrix of the grid
            # indexing='ij' means that the first two indices of the 3D matrix are the r and theta coordinates
qq       = np.meshgrid(r_c,theta_c,phi_c,indexing='ij')
            # Extract the coordinates (r,theta, z) from the 3D matrix
rr       = qq[0]  # final r coor defined by cell center
tt       = qq[1] # final theta coor defined by cell center, just for defining zr
zr       = np.pi/2.e0 - qq[1]    # z = pi/2 - theta, essentially frame rotated by 90 degrees, and z is from 0 to 0.7 radians


    # Density: dust values (for each species) for each cell in the grid
                                                                                                                                                             
        # number of dust species
ndustspec = 1

        # Disk parameters   (Science) - ( to be referenced/lorna?)
sigma_g0 =  10**3*(g/cm**2)   # gas surface density at 1 au
sigma_D0 = 0.01*10**3*(g/cm**2) #g/cm^2 # dust surface density at 1 au
h_r0 = 0.1*au # dust scale height at 1 au
pl_sig = -1.5# power law index for the dust surface density
pl_h  = 1.15 # power law index for the dust scale height

        # dust density function

sigma_D   = sigma_D0 * (rr/au)**pl_sig         # surface density : power law function of r
hh_r      = h_r0 * (rr/au)**pl_h       # scale height : power law function of r, dimensionless
hh       = hh_r * rr         # physical scale height
rho_D     = ( sigma_D / (np.sqrt(2.e0*np.pi)*hh) ) * np.exp(-(zr**2/hh_r**2)/2.e0)  # vertical Gaussian density profile


    # Star and planet parameters  
# Star parameters
#
mstar    = 1.65*M_sun
rstar    = 1.6*R_sun
tstar    = 7650*K
pstar    = np.array([37.2*au,0.,0.])


# Planet parameters
mplanet  = 3*M_jup
rplanet  = 1.17*R_jup
tplanet  = 1000
pplanet  = np.array([0.,0.,0.])

# Write to .inp files


    # Main setting

with open('radmc3d.inp','w+') as f: # Open the file in write-plus (can both read and write the file) mode
# with: ensures file is properly closed after the function is executed even if error occurs
    f.write('nphot = %d\n' %(nphot))  # means nphot = <value of nphot>, %d is placeholder later replaced by value of nphot, \n is new line 
    f.write('scattering_mode_max = 1\n') # if scattering opacity is included anywhere
    f.write('iranfreqmode = 1\n') # frequency grid is based on input files


    # Grid

with open('amr_grid.inp','w+') as f:
    f.write('1\n')                       # iformat
    f.write('0\n')                       # AMR grid style  (0=regular grid, no AMR)
    f.write('100\n')                     # Coordinate system: spherical
    f.write('0\n')                       # gridinfo
    f.write('1 1 0\n')                   # Include r,theta coordinates, phi not included 
    f.write('%d %d %d\n'%(nr,ntheta,1))  # Size of grid
    for value in r_i:
        f.write('%13.6e\n'%(value))      # X coordinates (cell walls position)
    for value in theta_i:
        f.write('%13.6e\n'%(value))      # Y coordinates (cell walls)
    for value in phi_i:
        f.write('%13.6e\n'%(value))      # Z coordinates (cell walls)



    # Density

with open('dust_density.inp','w+') as f:
    f.write('1\n')                       # Format number -> it has format so dont need labels like nphot
    f.write('%d\n'%(nr*ntheta*nphi))     # number of cells in the grid
    f.write('%d\n'%(ndustspec))          # how many species of dust 

    # WHy flatten even if grid is 3D? 
    data = rho_D.ravel(order='F')         # Fortran-like index ordering, flatten density structure to 1D
    data.tofile(f, sep='\n', format="%13.6e")  # 13 characters long, include integer and decider, 6 decimal long
    # sep = '\n' means each value is separated by a new line
    f.write('\n') # new line at the end of the file
        

    # Stars and Planets
with open('stars.inp','w+') as f:
    f.write('2\n')
    f.write('1 %d\n\n'%(nlam))
    f.write('%13.6e %13.6e %13.6e %13.6e %13.6e\n\n'%(rstar,mstar,pstar[0],pstar[1],pstar[2]))
    for value in lam:
        f.write('%13.6e\n'%(value))
    f.write('\n%13.6e\n'%(-tstar))