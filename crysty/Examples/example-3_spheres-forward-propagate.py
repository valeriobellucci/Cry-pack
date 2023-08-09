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
from utility import (wavelength, free_space_propagator, 
                     make_gauss2D, make_gauss2D_propagator)
from objects import (EllipseImageEnsemble,
                     create_ellipse_image as create_ellipse,
                     MaterialRefractiveIndex, rotate_points3D
                     )


#%%
"""
This section is an example of a rough visualization of an array of data,
representing an ensemble of speres in 3D space 
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



#%%
"""
This section rotates the array of points to simulate 
a detector positioned at a certain direction.
Finally, it gives a rough visualization. 
"""

# Given parameters for the rotation
axis = [0, 0, 1]  # rotation axis
rot_angle = 35.1  # rotation angle, in degrees

# Position of the spheres
v = positions  # for this example, assuming 'positions' is predefined

# Compute rotations
vrot_plus = rotate_points3D(v, axis, +rot_angle)
vrot_minus = rotate_points3D(v, axis, -rot_angle)

# Create a figure and plot the three sets of points
fig = plt.figure(figsize=(5, 5))

scale_factor = 5  # This factor may need to be adjusted based on the appearance
# The area of each point (s) is determined by 
# the square of its diameter scaled by the scale factor.
area = scale_factor * diameters ** 2

# Plotting vrot_plus. You can add more rotations as subplots.
ax1 = fig.add_subplot(111, projection="3d")  # 1x3 grid, 1st subplot
ax1.scatter3D(vrot_plus[:,0], 
              vrot_plus[:,1], 
              vrot_plus[:,2], 
              s = area,
              color="green")
ax1.set_title(f"Rotation Angle: +{rot_angle}°")
ax1.view_init(elev=0., azim=0.) # esplicitly set the angle of view

# Display the plots
plt.tight_layout()  # To ensure that the subplots don't overlap
plt.show()


#%%
"""
Creates an image of the thickness of the ensemple of spheres. 
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

# Generate an image of the ensemble.
img = ensemble(v2D_plus, radius, pixel_size=pix_real, sup_rate=sup_rate)

# Display the generated image.
plt.figure(figsize=(5,5))
img_show = img[::-1, ::-1]  # Invert the image for visualization purposes.
plt.imshow(img_show)
plt.show()

#%%

energy_eV = 12_000
energy_keV = energy_eV / 1000
wave_length = wavelength(energy_keV)


#%%

material = MaterialRefractiveIndex('Materials/refractive-index_soda-lime-glass.txt')
delta = material.get_beta(energy_eV)
beta = material.get_delta(energy_eV)
print(delta, beta)

#%%
ensemble = EllipseImageEnsemble()
img = ensemble(vrot_plus, radius, pixel_size=pix_real, sup_rate=2)


refractive_index = 1 - delta + 1j * beta

k_material = 2 * np.pi * refractive_index / wave_length # material wave number
transmitted_wavefield = np.exp( -1j * img * k_material)
transmitted_intensity = np.absolute(transmitted_wavefield)**2

plt.figure(figsize = (10,10))
plt.imshow(transmitted_intensity)
plt.show() 


#%%

distance = 1.30 # m 0.13
Kz = 2 * np.pi * 1.0 / wave_length # "vacuum" wave number
propagator = np.exp( 1j * Kz * distance)

result_fft = np.fft.fft2(transmitted_wavefield) * propagator
result = np.absolute( np.fft.ifft2(result_fft) )**2

plt.figure(figsize = (10,10))
plt.imshow(result)
plt.show() 

#%%

""" Example for creating a 2D mesh """
_, x_mesh, y_mesh = create_ellipse(500, radius, pixel_size=4.0e-6)

plt.imshow(y_mesh)
plt.show() 

plt.imshow(x_mesh)
plt.show()

#%%

pixel_size = pix_real
space_propagator = free_space_propagator(326, wave_length, distance, pixel_size, fresnel=False)

result_fft = np.fft.fft2(transmitted_wavefield) * space_propagator
result_fft = np.fft.fftshift(result_fft)
result_fft = result_fft[10:300,10:300]
propagated_wavefield = np.fft.ifft2(np.fft.fftshift(result_fft))
propagated_intensity = np.absolute( propagated_wavefield )**2

plt.figure(figsize = (5,5))
plt.imshow(propagated_intensity, vmin=0., vmax=1., cmap='gray')
plt.show() 

#%%

# Usage example:
gauss, _, _ = make_gauss2D(326, sigma=20, pixel_size = pix_real)
# Plotting the generated Gaussian
plt.figure(figsize = (5,5))
plt.imshow(gauss)
plt.show() 

gauss, _, _ = make_gauss2D(326, sigma=1, pixel_size = pix_real)

#%%

pixel_size = pix_real
resolution_propagator = make_gauss2D_propagator(326, sigma=1, pixel_size = pix_real)
space_propagator = free_space_propagator(326, wave_length, distance, pixel_size, fresnel=False)

result_fft = np.fft.fft2(transmitted_wavefield) * space_propagator * resolution_propagator
result_fft = np.fft.fftshift(result_fft)
result_fft = result_fft #[100:4900,100:4900]
propagated_wavefield = np.fft.ifft2(np.fft.fftshift(result_fft))
propagated_intensity = np.absolute( propagated_wavefield )**2



plt.figure(figsize = (5,5))
plt.imshow(propagated_intensity, vmin=0., vmax=np.max(result)*1, cmap='gray')
plt.show()

#%%

detail = propagated_intensity[150:200,100:150]
plt.figure(figsize = (5,5))
plt.imshow(detail, vmin=0., vmax=np.max(result)*1, cmap='gray')
plt.show() 

