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
from scipy.ndimage import gaussian_filter
from skimage.transform import resize
from utility import (transmitted_wavefield, free_space_propagator, 
                    intensity, section)
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
oversampling = 10 # Oversampling improves image resolution during post-processing.
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
This section calculates the refractive index of the obhect at a given energy (in eV). 
It then calculates the transmitted wavefield through an object.
"""

energy_eV = 12_000 # Photon energy (eV)

# Initialize a material using its refractive index data from the given file path.
material = MaterialRefractiveIndex('Materials/refractive-index_soda-lime-glass.txt')

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
zoom = transmitted_intensity[1000:1500, 1500:2000] # oversampling=10

# Show a section of the image
#section(zoom,25,10,50,100) # oversampling=2
section(zoom,125,50,250,500) # oversampling=10


#%%
"""
-----------------------------
Propagation of the Wavefield
-----------------------------
This section of the code performs the propagation of the wavefield 
using the Fourier method, and then visualizes the result.
"""

# Parameters
distance = 0.13  # Distance in meters between sample and detector
size = int(img.shape[0] / oversampling)  # Calculate the size of the original image

# Calculate the propagator for free space based on the given parameters
space_propagator = free_space_propagator(size, energy_eV, 
                                         distance, pix_real, 
                                         oversampling=oversampling,
                                         fresnel=False)

# Propagate the wavefield using FFT
result_fft = np.fft.fft2(img_wavefield) * space_propagator
propagated_wavefield = np.fft.ifft2(result_fft)

# Calculate the intensity of the propagated wavefield
propagated_intensity = intensity(propagated_wavefield)

# Display the propagated intensity
plt.figure(figsize=(5,5))
plt.imshow(propagated_intensity, vmin=0., vmax=1., cmap='gray')
plt.show()

# Extract a smaller section from the transmitted intensity for detailed viewing
zoom = propagated_intensity[1000:1500, 1500:2000] 
section(zoom, 125, 50, 250, 500)  # Display a section of the zoomed-in image

#%%
"""
-------------------------
Blurring and Downsampling
-------------------------
This section smoothes the propagated intensity by applying a Gaussian blur, 
then downsamples the image to the original pixel size, and visualizes the result.
"""

# Define the sigma for Gaussian blurring based on oversampling
sigma = oversampling / 10.0

# Blur the propagated intensity
img_blurred = gaussian_filter(propagated_intensity, sigma=sigma)

# Downsample the blurred image
img_downsampled = resize(img_blurred, (size, size), anti_aliasing=True)

# Display the downsampled image
plt.figure(figsize=(5,5))
plt.imshow(img_downsampled, vmin=0., vmax=1., cmap='gray')
plt.show()

# Extract a smaller section from the downsampled image for detailed viewing
zoom = img_downsampled[100:150, 150:200]
section(zoom, 12, 5, 25, 50)  # Display a section of the zoomed-in image
