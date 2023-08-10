"""
Created on Tue Aug  8 00:05:41 2023

@author: Valerio

This script is an example of a rough visualization of an array of data,
representing an ensemble of speres in 3D space 

"""

# Import necessary libraries
import pandas as pd
import matplotlib.pyplot as plt
from crysty.objects import  make_random_3Dspheres
                    

# Generate a pandas DataFrame with random 3D spheres data.
# The function will generate 5 random spheres with:
# - x, y, z coordinates ranging from -100 to 100
# - diameters ranging from 20 to 40
arr = make_random_3Dspheres(5, 100, (20, 40))

# Display the generated DataFrame to the console
print(arr)

# Extract X, Z, and Y columns and store them as a separate DataFrame called 'positions'
positions = arr[['X','Z','Y']]

# Separate each column from the dataframe for easier referencing
x = arr.X
y = arr.Y
z = arr.Z
diameters = arr.diameter

# Print the average diameter of the spheres in pixels
print("Average diameter is", diameters.mean(), "pixels")

# Creating figure
fig = plt.figure(figsize = (10, 7)) # specify size for plotting
ax = plt.axes(projection ="3d") # Add 3D axes

scale_factor = 5  # This factor may need to be adjusted based on the appearance

# The area of each point (s) is determined by 
# the square of its diameter scaled by the scale factor.
area = scale_factor * diameters ** 2

# Plot the data as 3D scatter points. 
ax.scatter3D(x,y,z, 
             s = area,
             color = "green", )
plt.title("simple 3D scatter plot") # plot title
 
# show plot
plt.show()

