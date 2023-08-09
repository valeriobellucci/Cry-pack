# -*- coding: utf-8 -*-
"""
Created on Wed May 16 18:07:38 2018

@author: xx9905
"""

import quantities as q
from numpy.fft import fftshift
from syris.physics import energy_to_wavelength
from syris.util import make_tuple
from concert.storage import read_image
from syris.bodies.simple import *
import numpy as np
import cry_functions as cry
import matplotlib.pyplot as plt



def make_4_spheres(n, radius, separation=0.75, pixel_size=1 * q.m, material=None, queue=None):
    """Make a sphere with image shape (*n*, *n*), *radius* and *pixel_size*. Sphere center is in n /
    2 + 0.5, which means between two adjacent pixels. *pixel_size*, *material* and *queue*, which is
    an OpenCL command queue, are used to create :class:`.StaticBody`.
    xcentre, ycentre control the position of the centre of the body in units of half of the field of
        view. 0 is the centre of the image, while +-1 are its borders.
    xdil, ydil control the dilatation of the body, making the sphere an ellipse.
    """
    radius = radius.simplified.magnitude
    image = np.zeros((n, n), dtype=cfg.PRECISION.np_float)
    image1 = image
    image2 = image
    image3 = image
    image4 = image
    
    # sphere 1
    x1_centre, y1_centre = separation, 0.0
    y1, x1 = (np.mgrid[(-n / 2 - y1_centre * n / 2) : (n / 2 - y1_centre * n / 2) : 1 , 
                     (-n / 2 - x1_centre * n / 2) : (n / 2 - x1_centre * n / 2) : 1 ])                     
    x1 = (x1 + .5) * pixel_size.simplified.magnitude
    y1 = (y1 + .5) * pixel_size.simplified.magnitude
    valid1 = np.where(x1 ** 2 + y1 ** 2 < radius ** 2)
    image1[valid1] = 2 * np.sqrt(radius ** 2 - x1[valid1] ** 2 - y1[valid1] ** 2)
    
    # sphere 2
    x2_centre, y2_centre = -separation, 0.0
    y2, x2 = (np.mgrid[(-n / 2 - y2_centre * n / 2) : (n / 2 - y2_centre * n / 2) : 1 , 
                     (-n / 2 - x2_centre * n / 2) : (n / 2 - x2_centre * n / 2) : 1 ])                     
    x2 = (x2 + .5) * pixel_size.simplified.magnitude
    y2 = (y2 + .5) * pixel_size.simplified.magnitude
    valid2 = np.where(x2 ** 2 + y2 ** 2 < radius ** 2)
    image2[valid2] = 2 * np.sqrt(radius ** 2 - x2[valid2] ** 2 - y2[valid2] ** 2)
    
    # sphere 3
    x3_centre, y3_centre = 0.0, separation
    y3, x3 = (np.mgrid[(-n / 2 - y3_centre * n / 2) : (n / 2 - y3_centre * n / 2) : 1 , 
                     (-n / 2 - x3_centre * n / 2) : (n / 2 - x3_centre * n / 2) : 1 ])                     
    x3 = (x3 + .5) * pixel_size.simplified.magnitude
    y3 = (y3 + .5) * pixel_size.simplified.magnitude
    valid3 = np.where(x3 ** 2 + y3 ** 2 < radius ** 2)
    image3[valid3] = 2 * np.sqrt(radius ** 2 - x3[valid3] ** 2 - y3[valid3] ** 2)
    
    # sphere 1
    x4_centre, y4_centre = 0.0, -separation
    y4, x4 = (np.mgrid[(-n / 2 - y4_centre * n / 2) : (n / 2 - y4_centre * n / 2) : 1 , 
                     (-n / 2 - x4_centre * n / 2) : (n / 2 - x4_centre * n / 2) : 1 ])                     
    x4 = (x4 + .5) * pixel_size.simplified.magnitude
    y4 = (y4 + .5) * pixel_size.simplified.magnitude
    valid4 = np.where(x4 ** 2 + y4 ** 2 < radius ** 2)
    image4[valid4] = 2 * np.sqrt(radius ** 2 - x4[valid4] ** 2 - y4[valid4] ** 2)
    
    image = image1 + image2 + image3 + image4

    return StaticBody(image * q.m, pixel_size, material=material, queue=queue)


def make_sphere_generic(n, radius, pixel_size=1 * q.m, material=None, queue=None, x_centre=0, y_centre=0, xdil=1, ydil=1):
    """Make a sphere with image shape (*n*, *n*), *radius* and *pixel_size*. Sphere center is in n /
    2 + 0.5, which means between two adjacent pixels. *pixel_size*, *material* and *queue*, which is
    an OpenCL command queue, are used to create :class:`.StaticBody`.
    xcentre, ycentre control the position of the centre of the body in units of half of the field of
        view. 0 is the centre of the image, while +-1 are its borders.
    xdil, ydil control the dilatation of the body, making the sphere an ellipse.
    """
    xdil = float(xdil)
    ydil = float(ydil)
    image = np.zeros((n, n), dtype=cfg.PRECISION.np_float)
    y, x = (np.mgrid[(-n / 2 - y_centre * n / 2) / ydil : (n / 2 - y_centre * n / 2) / ydil : 1 / ydil, 
                     (-n / 2 - x_centre * n / 2) / xdil : (n / 2 - x_centre * n / 2) / xdil : 1 / xdil])                     
    x = (x + .5) * pixel_size.simplified.magnitude
    y = (y + .5) * pixel_size.simplified.magnitude
    radius = radius.simplified.magnitude
    valid = np.where(x ** 2 + y ** 2 < radius ** 2)
    image[valid] = 2 * np.sqrt(radius ** 2 - x[valid] ** 2 - y[valid] ** 2)

    return StaticBody(image * q.m, pixel_size, material=material, queue=queue)


def make_disk_generic(n, radius, thickness, pixel_size=1 * q.m, material=None, queue=None, x_centre=0, y_centre=0, xdil=1, ydil=1):
    """Make a disk with image shape (*n*, *n*), *radius*, *thickness* and *pixel_size*. Sphere center is in n /
    2 + 0.5, which means between two adjacent pixels. *pixel_size*, *material* and *queue*, which is
    an OpenCL command queue, are used to create :class:`.StaticBody`.
    xcentre, ycentre control the position of the centre of the body in units of half of the field of
        view. 0 is the centre of the image, while +-1 are its borders.
    xdil, ydil control the dilatation of the body, making the disk an ellipse.
    """
    xdil = float(xdil)
    ydil = float(ydil)
    image = np.zeros((n, n), dtype=cfg.PRECISION.np_float)
    y, x = (np.mgrid[(-n / 2 - y_centre * n / 2) / ydil : (n / 2 - y_centre * n / 2) / ydil : 1 / ydil, 
                     (-n / 2 - x_centre * n / 2) / xdil : (n / 2 - x_centre * n / 2) / xdil : 1 / xdil])                     
    x = (x + .5) * pixel_size.simplified.magnitude
    y = (y + .5) * pixel_size.simplified.magnitude
    radius = radius.simplified.magnitude
    valid = np.where(x ** 2 + y ** 2 < radius ** 2)
    thickness = thickness.simplified.magnitude
    image[valid] = thickness

    return StaticBody(image * q.m, pixel_size, material=material, queue=queue)



def make_4_disks(n, radius, thickness, separation=0.75, pixel_size=1 * q.m, material=None, queue=None):
    """Make a sphere with image shape (*n*, *n*), *radius* and *pixel_size*. Sphere center is in n /
    2 + 0.5, which means between two adjacent pixels. *pixel_size*, *material* and *queue*, which is
    an OpenCL command queue, are used to create :class:`.StaticBody`.
    xcentre, ycentre control the position of the centre of the body in units of half of the field of
        view. 0 is the centre of the image, while +-1 are its borders.
    xdil, ydil control the dilatation of the body, making the sphere an ellipse.
    """
    radius = radius.simplified.magnitude
    thickness = thickness.simplified.magnitude
    image = np.zeros((n, n), dtype=cfg.PRECISION.np_float)
    image1 = image
    image2 = image
    image3 = image
    image4 = image
    # sphere 1
    x1_centre, y1_centre = separation, 0.0
    y1, x1 = (np.mgrid[(-n / 2 - y1_centre * n / 2) : (n / 2 - y1_centre * n / 2) : 1 , 
                     (-n / 2 - x1_centre * n / 2) : (n / 2 - x1_centre * n / 2) : 1 ])                     
    x1 = (x1 + .5) * pixel_size.simplified.magnitude
    y1 = (y1 + .5) * pixel_size.simplified.magnitude
    valid1 = np.where(x1 ** 2 + y1 ** 2 < radius ** 2)
    image1[valid1] = thickness
    
    # sphere 2
    x2_centre, y2_centre = -separation, 0.0
    y2, x2 = (np.mgrid[(-n / 2 - y2_centre * n / 2) : (n / 2 - y2_centre * n / 2) : 1 , 
                     (-n / 2 - x2_centre * n / 2) : (n / 2 - x2_centre * n / 2) : 1 ])                     
    x2 = (x2 + .5) * pixel_size.simplified.magnitude
    y2 = (y2 + .5) * pixel_size.simplified.magnitude
    valid2 = np.where(x2 ** 2 + y2 ** 2 < radius ** 2)
    image2[valid2] = thickness
    
    # sphere 3
    x3_centre, y3_centre = 0.0, separation
    y3, x3 = (np.mgrid[(-n / 2 - y3_centre * n / 2) : (n / 2 - y3_centre * n / 2) : 1 , 
                     (-n / 2 - x3_centre * n / 2) : (n / 2 - x3_centre * n / 2) : 1 ])                     
    x3 = (x3 + .5) * pixel_size.simplified.magnitude
    y3 = (y3 + .5) * pixel_size.simplified.magnitude
    valid3 = np.where(x3 ** 2 + y3 ** 2 < radius ** 2)
    image3[valid3] = thickness
    
    # sphere 1
    x4_centre, y4_centre = 0.0, -separation
    y4, x4 = (np.mgrid[(-n / 2 - y4_centre * n / 2) : (n / 2 - y4_centre * n / 2) : 1 , 
                     (-n / 2 - x4_centre * n / 2) : (n / 2 - x4_centre * n / 2) : 1 ])                     
    x4 = (x4 + .5) * pixel_size.simplified.magnitude
    y4 = (y4 + .5) * pixel_size.simplified.magnitude
    valid4 = np.where(x4 ** 2 + y4 ** 2 < radius ** 2)
    image4[valid4] = thickness
    
    image = image1 + image2 + image3 + image4

    return StaticBody(image * q.m, pixel_size, material=material, queue=queue)


def make_rectangle(intensity, pixel_number, x_pixels, y_pixels, x_centre=0, y_centre=0):
    """Make a rectangle image with intensity intensity_max image shape (*pixel_number*, *pixel_number*), 
    x pixel number *x_pixels* and y pixel number *y_pixels*. xcentre, ycentre control the position of 
    the centre of the body in units of half of the field of view. 0 is the centre of the image, 
    while +-1 are its borders.
    """
    rectangle = np.zeros((pixel_number, pixel_number))
    (rectangle[pixel_number * (1 + y_centre)/2 - y_pixels/2 : pixel_number * (1 + y_centre)/2 + y_pixels/2, 
               pixel_number * (1 + x_centre)/2 - x_pixels/2 : pixel_number * (1 + x_centre)/2 + x_pixels/2]) = intensity
    return rectangle






