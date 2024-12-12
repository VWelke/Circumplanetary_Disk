# The file aims to plot post-processed .out files by calling image function (2D intensity maps ) from simpleplot.py. 
# Given different amount of .out files to be plotted, this file might simply serve as the general template for plotting the data.
# Modifications will be made to make plots presentable


# Call the image function from simpleplot.py
import numpy as np
import sys 
import os
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

out_path = 'C:/Users/LHEM/radmc3d-2.0-master/radmc3d-2.0-master/examples/run_simple_1/image.out'  # Specify the path to the directory containing the .out files


# Check if the file exists in the specified directory
if os.path.isfile(out_path):
    print("File found in the specified directory.")
else:
    print("File not found in the specified directory.")

# Plot it
f = np.genfromtxt(out_path, delimiter=None)
print("Data loaded successfully.")

# Check the dimensions and structure of the loaded data
print("Data shape:", f.shape)
print("Data type:", f.dtype)
print("Data sample:\n", f[:5])  # Print the first 5 rows of the data

image(f)  # Call the image function to plot the data


