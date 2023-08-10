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
from crysty.utility import (transmitted_wavefield, 
                    intensity, section, Propagate)
from crysty.objects import (EllipseImageEnsemble, MaterialRefractiveIndex, 
                     rotate_points3D, make_random_3Dspheres)


#%%
"""
Load a dataset representing an ensemble of speres in 3D space.
"""

# Generate a pandas DataFrame with random 3D spheres data.
# The function will generate 20 random spheres with:
# - x, y, z coordinates ranging from -100 to 100
# - diameters ranging from 20 to 30
arr = make_random_3Dspheres(20, 100, (20, 30))

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

pixel_real = 4.0e-6 # Camera pixel (4 µm in meters).
oversampling = 10 # Oversampling improves image resolution during post-processing.
radius = diameters / 2 # Calculate the radius of each sphere using its diameter.

ensemble = EllipseImageEnsemble() # Create an ensemble of spheres for visualization
v2D_plus = vrot_plus[:,[0,2]] # Discard the Y (depth) coordinate

# Generate an image of the ensemble.
img = ensemble(v2D_plus, radius, pixel_size=pixel_real, oversampling=oversampling)

# Display the generated image.
plt.figure(figsize=(5,5))
plt.imshow(img)
plt.show()

#%%

"""
This section calculates the refractive index of the obhect at a given energy (in eV). 
It then calculates the transmitted wavefield through an object.
"""

energy_eV = 12_000 # Photon energy (eV)

# Initialize a material using its refractive index data from the given file path.
material = MaterialRefractiveIndex('refractive-index_soda-lime-glass.txt')

# Obtain the refractive index of the material for the specified energy.
refractive_index = material.get_refractive_index(energy_eV)

# Calculate the wavefield of the image after being transmitted through the object.
img_wavefield = transmitted_wavefield(img, refractive_index, energy_eV)

# Compute the intensity of the transmitted wavefield.
transmitted_intensity = intensity(img_wavefield)

# Display the transmitted intensity image.
plt.figure(figsize = (7,7))
plt.imshow(transmitted_intensity)
plt.show()

# Extract a smaller section from the transmitted intensity for detailed viewing.
#zoom = transmitted_intensity[200:300, 300:400] # oversampling=2
zoom = transmitted_intensity[2000:4000, 2000:4000] # oversampling=10

# Show a section of the image
#section(zoom,25,10,50,100) # oversampling=2
section(zoom, [1000,1000], [2000,50]) # oversampling=10


#%%
"""
-----------------------------
Propagation of the Wavefield
-----------------------------
This section of the code performs the propagation of the wavefield 
using the Fourier method, and then visualizes the result.
"""
distance = 0.13 # propagation distance in meters
experiment = Propagate(img_wavefield, distance, energy_eV, pixel_real,
             oversampling=10, fresnel=False, 
             blur=True, sigma=None, downsample=True,
             anti_aliasing=True)

img_result = experiment.get_propagated_intensity()

# Display the downsampled image
plt.figure(figsize=(5,5))
plt.imshow(img_result, vmin=0., vmax=1., cmap='gray')
plt.show()

# Extract a smaller section from the downsampled image for detailed viewing
zoom = img_result[200:400, 200:400]
section(zoom, [100,100], [200,10])  # Display a section of the zoomed-in image
