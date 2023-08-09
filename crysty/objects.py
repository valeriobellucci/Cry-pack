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
import math
from scipy.interpolate import interp1d
import pandas as pd

def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis / math.sqrt(np.dot(axis, axis))
    a = math.cos(theta / 2.0)
    b, c, d = -axis * math.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])

def rotate_points3D(points, axis, rot_angle):
    """
    Rotate a list of 3D points around a given axis by a specified angle.

    Parameters:
    - points (list or array-like): A list of 3D points to be rotated.
    - axis (list or array-like): A 3-element list or array indicating the axis of rotation.
    - rot_angle (float): The rotation angle in degrees.

    Returns:
    - numpy.ndarray: An array of rotated 3D points.
    """
    theta = np.deg2rad(rot_angle)
    rotated_points = np.dot(rotation_matrix(axis, theta), np.transpose(points))
    return np.transpose(rotated_points)

def create_ellipse_image(n, radius, 
                         pixel_size=1, x_centre=0, y_centre=0, 
                         xdil=1., ydil=1., oversampling=1):
    """
    Generate an image of a sphere or ellipse based on given parameters.
    
    Parameters:
        n (int): Dimensions of the input coordinates (n, n).
        radius (float): Radius of the sphere/ellipse, in pixels.
        pixel_size (float, optional): Size of each pixel in the image. Default is 1.
        x_centre (float, optional): X-coordinate of the ellipse's center, in pixels. 
                                    Default is 0 (center of the image).
        y_centre (float, optional): Y-coordinate of the ellipse's center, in pixels. 
                                    Default is 0 (center of the image).
        xdil (float, optional): X-axis dilation of the ellipse. Default is 1 (no dilation).
        ydil (float, optional): Y-axis dilation of the ellipse. Default is 1 (no dilation).
        oversampling (int, optional): oversampling rate. It is the ratio between 
                                    the size of the camera pixels over
                                    the size of the pixels used during computation 
    Returns:
        tuple: A tuple containing the generated image with size (n, n) * oversampling, 
        and X and Y coordinate arrays.
    """
    
    # Adjust for the oversampling
    n = n * oversampling
    x_centre = x_centre * oversampling
    y_centre = y_centre * oversampling
    radius = radius * oversampling
    pixel_size = pixel_size / oversampling
    
    # Initialize a blank image
    image = np.zeros((n, n), dtype="float64")
    
    # Generate X and Y coordinate grids
    y, x = np.mgrid[
        (-n/2 - y_centre)/ydil : (n/2 - y_centre)/ydil : 1/ydil, 
        (-n/2 - x_centre)/xdil : (n/2 - x_centre)/xdil : 1/xdil
    ]
    
    # Adjust the coordinates based on the pixel size
    x = (x + 0.5) * pixel_size
    y = (y + 0.5) * pixel_size
    radius = radius * pixel_size
    
    # Identify the valid coordinates that lie inside the sphere/ellipse boundary
    valid = np.where(x**2 + y**2 < radius**2)
    
    # Update the image values for the valid coordinates
    image[valid] = 2 * np.sqrt(radius**2 - x[valid]**2 - y[valid]**2)

    return image, x, y


class EllipseImageEnsemble:
    """
    Generates an image by overlaying multiple ellipses. 

    The ellipses can be highlighted based on their coordinates or index.

    Parameters:
        find_with_position (bool): Flag to decide if ellipses are to be highlighted based on their position.
        xy_find (1-d or 2-d list): List of [x,y] positions to be highlighted. Pixels units.
        highlight_position (list): Intensity multiplier for ellipses found by position.
        find_with_index (bool): Flag to decide if ellipses are to be highlighted based on their index.
        index_find (int or list): List of indices to be highlighted.
        highlight_index (int): Intensity multiplier for ellipses found by index.
    """
    
    def __init__(self, 
                 find_with_position=False, 
                 xy_find=[0., 0.], 
                 highlight_position=[10],
                 find_with_index=False, 
                 index_find=0,
                 highlight_index=5
                 ):
        self.find_with_position = find_with_position
        self.xy_find = [xy_find] if np.ndim(xy_find) == 1 else xy_find
        self.highlight_position = highlight_position
        self.find_with_index = find_with_index
        self.index_find = [index_find] if np.ndim(index_find) == 0 else index_find
        self.highlight_index = highlight_index
        
    def __call__(self, array, radius, n=None, pixel_size=1, oversampling=1):
        """
        Generates an image by overlaying ellipses based on the input parameters.

        Parameters:
            array (np.array): Array with shape (N,2) containing (x,y) coordinates of the centers of ellipses.
            radius (float or np.array): Radius of ellipses, in pixels. Can be a scalar or an array of length N.
            n (int, optional): Size of the resulting image. If not given, computed based on input array and radius.
            pixel_size (float): Pixel size for the image. Default is 1.
            oversampling (int, optional): oversampling rate. It is the ratio between 
                                        the size of the camera pixels over
                                        the size of the pixels used during computation
        Returns:
            np.array: Image with superimposed ellipses, size (n, n) * oversampling.
        """
        
        # Compute 'n' if not given
        if n is None:
            n = 2 * int(
                        np.amax(array) 
                         - np.amin(array) 
                         + np.max(radius)
                         ) 
        
        # Ensure 'radius' is an array of appropriate length
        if np.ndim(radius) == 0:
            radius = np.full((len(array),), radius)    
        
        # Ensure 'index_find' and 'xy_find' are lists
        if np.ndim(self.index_find) == 0:
            self.index_find = [self.index_find] 
        
        if np.ndim(self.xy_find) == 1:
            self.xy_find = [self.xy_find]
        
        # Initialize the image
        img = np.zeros((n*oversampling, n*oversampling), dtype="float64")

        # Generate ellipses and superimpose them
        for i_point, point in enumerate(array):
            xc = point[0]
            yc = point[1]
            r = radius[i_point]
            image_now, _, _ = create_ellipse_image(n, r, 
                                                   x_centre=xc, y_centre=yc, 
                                                   pixel_size=pixel_size,
                                                   oversampling=oversampling)
            
            # Highlight based on position
            if self.find_with_position and [xc, yc] in self.xy_find:
                image_now *= self.highlight_position
                print("Found: object with center ", (xc, yc))
            
            # Highlight based on index
            if self.find_with_index and i_point in self.index_find:
                image_now *= self.highlight_index
                print("Found: object with index ", i_point)

            img += image_now

        return img


class MaterialRefractiveIndex:
    def __init__(self, file_path):
        """
        Initialize the MaterialRefractiveIndex with data from a file.
        Get the table of refractive index (Delta, Beta) vs Energy (eV) 
        in a .txt file from https://henke.lbl.gov/optical_constants/getdb2.html
        :param file_path: Path to the .txt file containing energy, delta, and beta values.
        """
        # Read the data from the file using pandas
        df = pd.read_csv(file_path, sep='\s+', skiprows=2, names=['Energy', 'Delta', 'Beta'])

        # Extract the data
        self.energy = df['Energy'].values
        self.delta = df['Delta'].values
        self.beta = df['Beta'].values

        # Interpolate the data
        self.delta_func = interp1d(self.energy, self.delta, kind='cubic')
        self.beta_func = interp1d(self.energy, self.beta, kind='cubic')

    def get_delta(self, energy):
        """
        Returns the interpolated value of delta for a given energy.

        :param energy: Energy (eV) value to find the corresponding delta.
        :return: Interpolated delta value.
        """
        return self.delta_func(energy)

    def get_beta(self, energy):
        """
        Returns the interpolated value of beta for a given energy.

        :param energy: Energy (eV) value to find the corresponding beta.
        :return: Interpolated beta value.
        """
        return self.beta_func(energy)
    
    def get_refractive_index(self, energy):
        """
        Returns the interpolated value of the refractive index for a given energy.
            refractive_index = 1 - delta - 1j * beta

        :param energy: Energy (eV) value to find the corresponding refractive index.
        :return: Interpolated refractive index value.
        """
        refractive_index = 1 - self.delta_func(energy) + 1j * self.beta_func(energy)

        return refractive_index



def pyopencl_array_from_numpy(data, queue=None):
    """ Turn a 2D numpy float array into a pyopencl.array.Array"""
    
    if isinstance(data, np.ndarray):
        if len(data.shape) != 2 or data.dtype.kind != 'f':
            raise TypeError("Only 2D arrays of real numbers are supported")
        
        return cl_array.to_device(queue, data.astype(np.float32))
    
    else:
        raise TypeError('Unsupported data type {}'.format(type(data)))


pass