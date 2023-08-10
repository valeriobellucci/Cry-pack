"""
Crystal Transfer Function Visualization Script
==============================================

This script is designed to visualize the Fourier transfer function's intensity 
for a given crystal structure, specifically for Silicon (Si) in this case. 

Key Features:
-------------
- Initialization of a Silicon (Si) crystal with specific orientation, asymmetry angle, and photon energy.
- Calculation of spatial frequencies (qx and qy) based on defined pixel and detector sizes.
- Generation of a 2D mesh grid of these spatial frequencies.
- Computation and visualization of the intensity of the crystal's Fourier transfer function.

Usage:
------
Run the script directly to compute and visualize the Fourier transfer function's intensity.

"""

# Import necessary libraries
import numpy as np
import quantities as q
from crysty.utility import screen, intensity
from crysty.crystal_class import Crystal

# Initialize a Silicon (Si) crystal with the following specifications:
# - Crystal orientation: [2, 2, 0]
# - Asymmetry angle: 0 degrees
# - Photon energy: 29.85 keV
crystal_1 = Crystal('Si', [2, 2, 0], 0. * q.deg, 29.85 * q.keV)

# Define the effective pixel size (simplified for now as 1 um).
pixel_size_effective = 55. * q.um / 55.

# Define the effective size of the detector, given as 512 times the pixel size.
detector_size_effective = pixel_size_effective * 512

# Calculate the spatial frequencies (qx) for image reconstruction.
# The maximum qx is determined by the Nyquist frequency.
qx_max = (1. / (2. * pixel_size_effective)).simplified.magnitude
qx_step = (1. / (2. * detector_size_effective)).simplified.magnitude
qx = np.arange(- qx_max / 2., + qx_max / 2., qx_step ) * (1. / q.m)

# As the grid is symmetric, qy values are the same as qx values.
qy = qx

# Create a 2D mesh grid of QX and QY values for spatial frequencies.
QX, QY = np.meshgrid(qx.simplified.magnitude, qy.simplified.magnitude) * (1. / q.m)

# Compute the intensity of the crystal's Fourier transfer function and display it.
screen(intensity(crystal_1.QX_or_QY_transfer_fourier(QY, 29.85 * q.keV)))
