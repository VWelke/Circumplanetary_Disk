# This file converts RADMC3d output in the form of .out files to FITS files.
# Reference:
# Author: charango
# Repository: https://github.com/charango/fargo2radmc3d/tree/master
# File: https://github.com/charango/fargo2radmc3d/blob/master/radmc_to_fits.py



# Read the data into fits file
from astropy.io import fits

output_directory = 'C:/Users/LHEM/Desktop/MSci Project/Images_RADMC_3D_fits/'  # Specify the directory to save the fits file

output_filename = f'image_out.fits'

output_path_fits = os.path.join(output_directory, output_filename)

# Save the reshaped data as a FITS file
hdu = fits.PrimaryHDU(reshaped_data) # Create a new HDU object to store the data, where PrimaryHDU is the main data structure with header and data
hdul = fits.HDUList([hdu]) # Create a new HDUList object to store the HDU
hdul.writeto(output_path_fits, overwrite=True)  