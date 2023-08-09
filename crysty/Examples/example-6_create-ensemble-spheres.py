#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 00:05:41 2023

@author: Valerio

Imports data about an ensemble of spheres and creates an image of the thickness.

"""

# Import necessary libraries
import pandas as pd
import matplotlib.pyplot as plt
from objects import (EllipseImageEnsemble,rotate_points3D)


#%%
"""
Import and elaborate an array of data representing an ensemble of speres in 3D space 
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


#%%
"""
Rotate the array of points to simulate a detector positioned at a certain direction.
"""

# Given parameters for the rotation
axis = [0, 0, 1]  # rotation axis
rot_angle = 35.1  # rotation angle, in degrees

# Position of the spheres
v = positions

# Compute rotation
vrot_plus = rotate_points3D(v, axis, +rot_angle)


#%%
"""
Creates an image of the thickness of the ensemple of spheres 
and shows some visualization options.
"""

# Camera settings and parameters:

# Define the physical size of one camera pixel (4 µm in meters).
pix_real = 4.0e-6 

# Super-resolution rate is the ratio between the size of the camera pixels 
# and the size of the pixels used during computation. A higher super-resolution 
# rate can improve image resolution during post-processing.
sup_rate = 2

# Calculate the radius of each sphere using its diameter.
radius = diameters / 2 

# Create an ensemble of spheres for visualization:
ensemble = EllipseImageEnsemble()

# Filter the sphere positions to only X and Z dimensions (essentially projecting 
# the 3D positions of the spheres onto a 2D plane by discarding the Y (depth) coordinate).
v2D_plus = vrot_plus[:,[0,2]]

# Generate an image of the ensemble without any special visualization features.
img = ensemble(v2D_plus, radius, pixel_size=pix_real, sup_rate=1)

# Display the generated image.
plt.figure(figsize=(5,5))
img_show = img[::-1, ::-1]  # Invert the image for visualization purposes.
plt.imshow(img_show)
plt.show()

# Generate an image of the ensemble highlighting spheres by their indices.
ensemble.find_with_index = True
ensemble.index_find = [4, 6]        # Spheres to highlight by index.
ensemble.highlight_index = 15       # Define the highlight intensity.
img = ensemble(v2D_plus, radius, pixel_size=pix_real, sup_rate=sup_rate)

# Display the highlighted image.
plt.figure(figsize=(5,5))
img_show = img[::-1, ::-1]  # Invert the image for visualization purposes.
plt.imshow(img_show)
plt.show()

# Reset the "find by index" flag.
ensemble.find_with_index = False

# Generate an image of the ensemble highlighting spheres by their positions.
ensemble.find_with_position = True
ensemble.xy_find = v2D_plus[1:4,:2]        # Positions of spheres to highlight.
ensemble.highlight_position = 10            # Define the highlight intensity.
img = ensemble(v2D_plus, radius, pixel_size=pix_real, sup_rate=sup_rate)

# Display the highlighted image.
plt.figure(figsize=(5,5))
img_show = img[::-1, ::-1]  # Invert the image for visualization purposes.
plt.imshow(img_show)
plt.show()
