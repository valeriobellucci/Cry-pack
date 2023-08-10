"""
This script propagate the image of an object through a crystal
"""

#%% 
# ---------------------------------
# Load packages 
# ---------------------------------

import numpy as np
import quantities as q
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
from skimage.transform import resize
from crysty.utility import (screen, intensity, transmitted_wavefield, 
                            section, free_space_propagator)
from crysty.crystal_class import Crystal
from crysty.objects import create_ellipse_image, MaterialRefractiveIndex


#%%

# ---------------------------------
# Parameters of the setup 
# ---------------------------------
 
n = 512 # detector pixels ( n x n )
pixel_real = 1.0e-6 # Camera pixel size in meters.
energy_eV = 10_000 # photon energy in eV
oversampling = 1 # Oversampling improves image resolution during post-processing.

#%%

# ---------------------------------
# Object 
# ---------------------------------

# Create a sample object, a sphere.
radius = 50 # radius of the sample object (sphere) in pixels
xc, yc = (0, 0) # center of the sample object

# Create the projected thickness image of the sample
img_sample, _, _ = create_ellipse_image(n, radius, pixel_size=pixel_real,
                                        x_centre=xc, y_centre=yc, 
                                        oversampling=oversampling)

# Display the projected thickness image
section(img_sample, [250,250], [500,20])  # Display a section 

material = MaterialRefractiveIndex('refractive-index_soda-lime-glass.txt')

# Obtain the refractive index of the material for the specified energy.
refractive_index = material.get_refractive_index(energy_eV)

# Calculate the wavefield of the image after being transmitted through the object.
img_wavefield = transmitted_wavefield(img_sample, refractive_index, energy_eV)
img = intensity(img_wavefield)
section(img, [250,250], [500,20])  # Display a section 

#%%

# ---------------------------------
# Crystal 
# ---------------------------------

# Initialize a Silicon crystal with specific parameters:
# - Material: Silicon (Si)
# - Crystal orientation: (220)
# - Asymmetry angle: 0 degrees
# - Photon energy: `energy_eV`
crystal = Crystal('Si', [2, 2, 0], 0. * q.deg, energy_eV * q.eV)

# Calculate the effective pixel size by dividing the real pixel size by the oversampling factor.
pixel_size_effective = pixel_real / oversampling 

# Calculate the effective size of the detector based on the effective pixel size and number of pixels.
detector_size_effective = pixel_size_effective * n * oversampling

# Calculate the maximum q-value (spatial frequency) and step size based on the pixel size.
qx_max = (1. / (2. * pixel_size_effective))
qx_step = (1. / (2. * detector_size_effective))

# Create an array of q-values (spatial frequencies) centered around zero.
qx = np.arange(- qx_max / 2., + qx_max / 2., qx_step ) 

# Set qy values to be the same as qx (since the spatial frequencies are the same in both directions).
qy = qx

# Create a 2D mesh grid of QX and QY values.
QX, QY = np.meshgrid(qx, qy) * (1. / q.m)

# Compute the Fourier transform of the crystal's transfer function for given q-values and energy.
crystal_transfer_fft = crystal.QX_or_QY_transfer_fourier(QY, energy_eV * q.eV)

# Display the computed intensity of the crystal's transfer function.
screen(intensity(crystal_transfer_fft))


#%%

# -------------------------------------------------------
# Free- space between Object - Crystal - Detector 
# -------------------------------------------------------

distance_object_to_crystal = 0.2  # Distance in meters between sample and crystal
distance_crystal_to_detector = 0.2  # Distance in meters between crystal and detector

# Calculate the propagator for free space between object and crystal
space_obj_to_cry = free_space_propagator(n, energy_eV, 
                                         distance_object_to_crystal, 
                                         pixel_real, 
                                         oversampling=oversampling,
                                         fresnel=False)

# Calculate the propagator for free space between crystal and detector
space_cry_to_det = free_space_propagator(n, energy_eV, 
                                         distance_crystal_to_detector, 
                                         pixel_real, 
                                         oversampling=oversampling,
                                         fresnel=False)

#%%
# ---------------------------------
# PROPAGATION 
# ---------------------------------

# Propagate the wavefield between object and crystal using FFT
result_obj_to_cry_fft = np.fft.fft2(img_wavefield) * space_obj_to_cry
section(intensity(np.fft.ifft2(result_obj_to_cry_fft)), [250,250], [500,20])  # Display a section 


# Propagate the wavefield after the crystal using FFT            
result_obj_after_cry_fft = result_obj_to_cry_fft * np.fft.ifftshift(crystal_transfer_fft)
section(intensity(np.fft.ifft2(result_obj_after_cry_fft)), [250,250], [500,20])  # Display a section

# Propagate the wavefield between the crystal and detector using FFT            
result_detector_fft = result_obj_after_cry_fft * space_cry_to_det
result_detector = np.fft.ifft2(result_detector_fft)
propagated_intensity = intensity(result_detector)
section(propagated_intensity, [250,250], [500,20])  # Display a section


#%%

# ---------------------------------
# Smoothing and downsampling 
# ---------------------------------

# Define the sigma for Gaussian blurring based on oversampling
sigma = oversampling / 2.0

# Blur the propagated intensity
img_blurred = gaussian_filter(propagated_intensity, sigma=sigma)

# Downsample the blurred image
img_downsampled = resize(img_blurred, (n, n), anti_aliasing=True)

# Display the downsampled image
plt.figure(figsize=(5,5))
plt.imshow(img_downsampled, cmap='gray')
plt.show()

# Extract a smaller section from the downsampled image for detailed viewing
zoom = img_downsampled[150:350, 150:350]
section(zoom, [100,100], [200,20])  # Display a section of the zoomed-in image


