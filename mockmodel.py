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

# Density: dust values (for each species) for each cell in the grid

# Number of cells in the grid is defined by




# Write to .inp files


# Main setting

with open('radmc3d.inp','w+') as f: # Open the file in write-plus (can both read and write the file) mode
# with: ensures file is properly closed after the function is executed even if error occurs
    f.write('nphot = %d\n' %(nphot))  # means nphot = <value of nphot>, %d is placeholder later replaced by value of nphot, \n is new line 
    f.write('scattering_mode_max = 1\n') # if scattering opacity is included anywhere
    f.write('iranfreqmode = 1\n') # frequency grid is based on input files


# Grid


# Density

with open('dust_density.inp','w+') as f:
    f.write('1\n')                       # Format number
    f.write('%d\n'%(nr*ntheta*nphi))     # number of cells in the grid
    f.write('ndust%d\n'%(ndustspec))                       # how many species of dust 
    data = rhod.ravel(order='F')         # Create a 1-D view, fortran-style indexing ?
    data.tofile(f, sep='\n', format="%13.6e")  #?
    f.write('\n')
        


