#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 00:05:41 2023

@author: Valerio
"""

# Import necessary libraries
import pandas as pd
import matplotlib.pyplot as plt
from objects import rotate_points3D


#%%
"""
This script import an array of data representing an ensemble of speres in 3D space.
Then, it rotates them to simulate a detector positioned over a certain direction.
Finally, it gives a rough visualization. 
"""

# Load data from the specified text file into a pandas DataFrame. 
# The data is tab-separated.
arr = pd.read_csv('3D-spheres/3D-spheres.txt', sep='\t')

# Display the loaded data
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


# Given parameters for the rotation
axis = [0, 0, 1]  # rotation axis
rot_angle = 35.1  # rotation angle, in degrees

# Position of the spheres
v = positions  # for this example, assuming 'positions' is predefined

# Compute rotations
vrot = rotate_points3D(v, axis, 2 * rot_angle)
vrot_plus = rotate_points3D(v, axis, +rot_angle)
vrot_minus = rotate_points3D(v, axis, -rot_angle)

# Create a figure and plot the three sets of points
fig = plt.figure(figsize=(15, 5))

scale_factor = 5  # This factor may need to be adjusted based on the appearance
# The area of each point (s) is determined by 
# the square of its diameter scaled by the scale factor.
area = scale_factor * diameters ** 2

# Plotting vrot_plus
ax1 = fig.add_subplot(131, projection="3d")  # 1x3 grid, 1st subplot
ax1.scatter3D(vrot_plus[:,0], 
              vrot_plus[:,1], 
              vrot_plus[:,2], 
              s = area,
              color="green")
ax1.set_title(f"Rotation Angle: +{rot_angle}°")
ax1.view_init(elev=0., azim=0.)

# Plotting vrot
ax2 = fig.add_subplot(132, projection="3d")  # 1x3 grid, 2nd subplot
ax2.scatter3D(vrot[:,0], 
              vrot[:,1], 
              vrot[:,2], 
              s = area,
              color="blue")
ax2.set_title(f"Rotation Angle: {2*rot_angle}°")
ax2.view_init(elev=0., azim=0.)

# Plotting vrot_minus
ax3 = fig.add_subplot(133, projection="3d")  # 1x3 grid, 3rd subplot
ax3.scatter3D(vrot_minus[:,0], 
              vrot_minus[:,1], 
              vrot_minus[:,2], 
              s = area,
              color="red")
ax3.set_title(f"Rotation Angle: -{rot_angle}°")
ax3.view_init(elev=0., azim=0.)

# Display the plots
plt.tight_layout()  # To ensure that the subplots don't overlap
plt.show()
