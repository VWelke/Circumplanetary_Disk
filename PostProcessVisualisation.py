# The file aims to plot post-processed .out files by calling image function (2D intensity maps ) from simpleplot.py. 
# Given different amount of .out files to be plotted, this file might simply serve as the general template for plotting the data.
# Modifications will be made to make plots presentable


# Call the image function from simpleplot.py
import numpy as np
import sys 
import os
import matplotlib.pyplot as plt

sys.path.append('C:/Users/LHEM/radmc3d-2.0-master/radmc3d-2.0-master/python/radmc3d_tools')  # Add the path to the directory containing simpleplot.py

# Check if the file exists in the specified directory
if os.path.isfile('C:/Users/LHEM/radmc3d-2.0-master/radmc3d-2.0-master/python/radmc3d_tools/simpleplot.py'):
    print("simpleplot.py found in the specified directory.")
else:
    print("simpleplot.py not found in the specified directory.")

# Import the image function from simpleplot.py
try:
    from simpleplot import image
    print("Import successful.")
except ImportError as e:
    print("Import failed:", e)

# Load the data from the .out files
# If there are multiple .out files, load them all in a list

out_path = 'C:/Users/LHEM/Desktop/Van_Code_Projects/Circumplanetary_Disk/CPD_simple_1/image.out'  # Specify the path to the directory containing the .out files


# Check if the file exists in the specified directory
if os.path.isfile(out_path):
    print("File found in the specified directory.")
else:
    print("File not found in the specified directory.")

# Plot it
# As a typical image.out have have its actual data starting from the 6th row   [unit: erg/s/cm^2/Hz/ster]
# The first 5 rows are the header information: iformat, number of pixels(x,y), number of wavlength images, pixel size(cm,x,y), wavelength (micron)


with open(out_path, 'r') as file: # opens the file in read mode
    lines = file.readlines() #returns a list of all lines in the file

# Parse header information
grid_size = list(map(int, lines[1].strip().split())) # Grid size in x and y directions
pixel_size_cm = float(lines[3].strip().split()[0])  # Pixel size in cm

# Convert pixel size to au
pixel_size_au = pixel_size_cm / 1.496e13

# Generate x and y coordinates in au
x_coords = np.linspace(0, grid_size[0] * pixel_size_au, grid_size[0])
y_coords = np.linspace(0, grid_size[1] * pixel_size_au, grid_size[1])

data_values = lines[5:]  # Numerical data

# Convert the data into a numpy array and reshape
data_values = np.array([float(value.strip()) for value in data_values if value.strip()]) # Convert to float and remove empty strings
reshaped_data = data_values.reshape(grid_size)  # Reshape the 1d array into a 2d array based

# Print the reshaped data's size and dimensions
print(f"Reshaped data dimensions: {reshaped_data.shape}")
print(reshaped_data)


# Plot the image
fig, ax = plt.subplots()
cax = ax.imshow(reshaped_data, extent=[x_coords.min(), x_coords.max(), y_coords.min(), y_coords.max()], cmap='inferno', aspect='equal')

# Add labels and title
plt.colorbar(cax, label='Intensity')  # Add colorbar for context
plt.xlabel('X (au)')
plt.ylabel('Y (au)')
plt.title('RADMC-3D Image')

# Save the image to a directory as png file
# plt.savefig('C:/Users/LHEM/Desktop/MSci Project/Images_RADMC_3D_png/image.png')  # Save the plot as a .png file

plt.show()