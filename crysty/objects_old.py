#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  7 15:54:23 2023

@author: Valerio
"""

import numpy as np
import quantities as q
from utility import make_couple, wavelength
import pyopencl.array as cl_array

def pyopencl_array_from_numpy(data, queue=None):
    """ Turn a 2D numpy float array into a pyopencl.array.Array"""
    
    if isinstance(data, np.ndarray):
        if len(data.shape) != 2 or data.dtype.kind != 'f':
            raise TypeError("Only 2D arrays of real numbers are supported")
        
        return cl_array.to_device(queue, data.astype(np.float32))
    
    else:
        raise TypeError('Unsupported data type {}'.format(type(data)))


class Object():
    """ Class representing a transmission/phase object """
    
    def __init__(self, thickness, pixel_size, material=None):
        self.material = material
        self.thickness = pyopencl_array_from_numpy(thickness.simplified.magnitude)
        self.pixel_size = make_couple(pixel_size)
        
    def _transfer(self, shape, pixel_size, energy, offset, t=None):
        """Transfer function implementation based on a refractive index."""
        refractive_index = self.material.get_refractive_index(energy)
        wlength = wavelength(energy/1000) # energy in eV
        projected_thickness = self.thickness

        return transfer(proj, refractive_index, wlength)

class Sphere(Object):
    """ 
    Represents a spherical object.

    Creates an image of a sphere with dimensions (n, n) based on the specified radius 
    and pixel size. The center of the sphere is positioned at n / 2 + 0.5, implying 
    that it is located between two adjacent pixels.
    """
    
    def __init__(self, n, radius, pixel_size=1 * q.m, material=None):
        super(Object, self).__init__(pixel_size, material)
        self.n = n
        self.radius = radius
    
        self.image = np.zeros((self.n, self.n), dtype=np.float32)
        y, x = np.mgrid[-self.n / 2:self.n / 2, -self.n / 2:self.n / 2]
        x = (x + .5) * self.pixel_size.simplified.magnitude
        y = (y + .5) * self.pixel_size.simplified.magnitude
        self.radius = self.radius.simplified.magnitude
        valid = np.where(x ** 2 + y ** 2 < self.radius ** 2)
        self.image[valid] = 2 * np.sqrt(self.radius ** 2 - x[valid] ** 2 - y[valid] ** 2)
    
        self.thickness = pyopencl_array_from_numpy(self.image)


pass