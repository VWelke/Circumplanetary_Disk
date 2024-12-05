import numpy as np
# Reference files from radmc3d-2.0/examples/run_ppdisk_simple_1/problem_setup.py


# Import some astronomical constants
from astropy.constants import M_sun, L_sun, R_sun,T_sun, au, pc, G

# Convert constants to cgs units
au = au.cgs.value     # Astronomical Unit       [cm]
pc = pc.cgs.value     # Parsec                  [cm]
M_sun = M_sun.cgs.value  # Solar mass              [g]
T_sun = T_sun.value      # Solar temperature       [K]
L_sun = L_sun.cgs.value  # Solar luminosity        [erg/s]
R_sun = R_sun.cgs.value  # Solar radius            [cm]





# Define the parameters of the model

    # radmc3d.inp parameter : main settings for RADMC-3D
nphot    = 1000000  #for the thermal monte carto simulation


    # Grid : defines layout of space

nr       = 32 
ntheta   = 32
nphi     = 1
r_in      = 10*au   # disk inner radius
r_out     = 100*au # disk outer radius 
theta_up  = 2*np.pi

        # Coordinate array

r_i       = np.logspace(np.log10(r_in),np.log10(r_out),nr+1)  #+1 because it is not counting cell centers, but the walls
theta_i   = np.linspace(0.e0,np.pi*2.e0,ntheta+1)
phi_i     = np.linspace(0.e0,np.pi*2.e0,nphi+1)

        # Cell center position array
r_c = 0.5*(r_i[0:nr]+r_i[1:nr+1])  # average of the cell wall x-pos and cell wall x+1 pos for each cell in the array
theta_c = 0.5*(theta_i[0:ntheta]+theta_i[1:ntheta+1])
phi_c = 0.5*(phi_i[0:nphi]+phi_i[1:nphi+1])


        # Make the grid
            # takes in the center positions of the cells and returns a 3D matrix of the grid
            # indexing='ij' means that the first two indices of the 3D matrix are the r and theta coordinates
qq       = np.meshgrid(rc,thetac,phic,indexing='ij')
            # Extract the coordinates (r,theta, z) from the 3D matrix
rr       = qq[0]  # final r coor defined by cell center
tt       = qq[1] # final theta coor defined by cell center
zr       = np.pi/2.e0 - qq[1]    # z = pi/2 - theta, z is the vertical coordinate


    # Density: dust values (for each species) for each cell in the grid
                                                                                                                                                             
        # number of dust species
ndustspec = 1
        # dust density function





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
    data = rhod.ravel(order='F')         # Fortran-like index ordering, flatten density structure to 1D
    data.tofile(f, sep='\n', format="%13.6e")  # 13 characters long, include integer and decider, 6 decimal long
    # sep = '\n' means each value is separated by a new line
    f.write('\n') # new line at the end of the file
        

