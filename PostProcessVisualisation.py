# The file aims to plot post-processed .out files by calling image function (2D intensity maps ) from simpleplot.py. 
# Given different amount of .out files to be plotted, this file might simply serve as the general template for plotting the data.
# Modifications will be made to make plots presentable


# Call the image function from simpleplot.py
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