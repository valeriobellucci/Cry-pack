#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 00:05:41 2023

@author: Valerio
"""

# Import necessary libraries
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from utility import transmitted_wavefield, intensity
from objects import (EllipseImageEnsemble, MaterialRefractiveIndex, 
                     rotate_points3D)


#%%
"""
Load a dataset representing an ensemble of speres in 3D space.
"""

# Load data from the specified text file into a pandas DataFrame. 
# The data is tab-separated.
arr = pd.read_csv('3D-spheres/3D-spheres.txt', sep='\t')

# Extract X, Z, and Y columns and store them as a separate DataFrame called 'positions'
positions = arr[['X','Z','Y']]

# Separate each column from the dataframe for easier referencing
x = arr.X
y = arr.Y
z = arr.Z
diameters = arr.diameter


#%%
"""
Rotate the array of points to simulate the detector positioned at a certain angle.
"""
axis = [0, 0, 1]  # rotation axis
rot_angle = 35.1  # rotation angle, in degrees
vrot_plus = rotate_points3D(positions, axis, +rot_angle) # Compute rotations


#%%
"""
Creates an image of the thickness of the ensemple of spheres. 
"""

pix_real = 4.0e-6 # Camera pixel (4 µm in meters).
oversampling = 2 # Oversampling improves image resolution during post-processing.
radius = diameters / 2 # Calculate the radius of each sphere using its diameter.

ensemble = EllipseImageEnsemble() # Create an ensemble of spheres for visualization
v2D_plus = vrot_plus[:,[0,2]] # Discard the Y (depth) coordinate

# Generate an image of the ensemble.
img = ensemble(v2D_plus, radius, pixel_size=pix_real, oversampling=oversampling)

# Display the generated image.
plt.figure(figsize=(5,5))
plt.imshow(img)
plt.show()

#%%

"""
This script calculates and prints the refractive index of soda-lime glass 
at a given energy (in eV). It then calculates the wavefield of an image 
after being transmitted through an object of the material and computes 
the corresponding intensity.
"""

# Photon energy, given in electron volts (eV).
energy_eV = 12_000

# Initialize a material using its refractive index data from the given file path.
material = MaterialRefractiveIndex('Materials/refractive-index_soda-lime-glass.txt')

# Obtain the refractive index of the material for the specified energy.
refractive_index = material.get_refractive_index(energy_eV)

# Print the refractive index and the corresponding energy in keV.
print(f"refractive_index={refractive_index:.3} for Energy={energy_eV/1000} keV ")

# Calculate the wavefield of the image after being transmitted through the object.
img_wavefield = transmitted_wavefield(img, refractive_index, energy_eV)

# Compute the intensity of the transmitted wavefield.
transmitted_intensity = intensity(img_wavefield)

# Display the transmitted intensity image.
plt.figure(figsize = (7,7))
plt.imshow(transmitted_intensity)
plt.show()

# Extract a smaller section from the transmitted intensity for detailed viewing.
zoom = transmitted_intensity[200:300, 300:400]
plt.figure(figsize = (7,7))
plt.imshow(zoom)

# Define the slice of the image over which to compute the section.
yc = 20
yh = 10
xc = 50
xw = 100

# Draw a rectangle on the zoom image to represent the section computed for the line plot.
# Rectangle coordinates: (x_start, y_start), width, height
rectangle = plt.Rectangle((xc-xw/2, yc-yh/2), xw, yh, edgecolor='red', facecolor='none')
plt.gca().add_patch(rectangle)

plt.show()

# Extract the band from the zoomed image to compute the line profile.
band = zoom[int(yc-yh/2) : int(yc+yh/2), int(xc-xw/2) : int(xc+xw/2)]

# Compute the mean of the band along the y-axis to generate a line profile.
line = np.mean(band, axis=0)

# Plot the line profile.
plt.plot(line)
plt.show()


