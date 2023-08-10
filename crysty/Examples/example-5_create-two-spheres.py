#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 16:49:29 2023

@author: Valerio

Test image composed of just two spheres

"""

# Import necessary libraries
import matplotlib.pyplot as plt
from crysty.objects import create_ellipse_image as create_ellipse

#%%

# Create an ellipse image with a size of 500x500 pixels and a radius of 25 pixels.
# The center of the ellipse is at the center of the image (0,0).
img_1, _, _ = create_ellipse(500, 25, x_centre=0., y_centre=0.)

# Create another ellipse image with the same size and radius.
# The center of this ellipse is at (50,150) in the image.
img_2, _, _ = create_ellipse(500, 25, x_centre=50, y_centre=150)

# Combine (overlay) the two ellipse images.
# Pixel values from both images are added together at each position.
img_test = img_1 + img_2

# Display the combined image using matplotlib.
plt.imshow(img_test)
plt.show()
