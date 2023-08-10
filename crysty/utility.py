# -*- coding: utf-8 -*-
"""
Created on Wed Jun  6 14:00:50 2018

@author: Valerio
"""

import numpy as np
import quantities as q
import scipy
import scipy.constants as sci # Collection of constants
import matplotlib.pyplot as plt
from skimage.transform import resize
from scipy.ndimage import gaussian_filter

""" Basic physical constants """
electron_radius = sci.physical_constants['classical electron radius'][0]  # 2.81794*10**-15 m
h_planck = sci.physical_constants['Planck constant in eV s'][0]  # 4.1356692*10**-15 eV s
polariz = 1.0  # Polarization factor now
# sci.c = c = speed of light in vacuum = 2.99792458*10**8 m/s
# electronic density
# sci.m_e = me = electron mass = 9.109*10**-31
# sci.epsilon_0 = Epo = the electric constant (vacuum permittivity) = 8.854*10**-12
# sci.e = elementary charge = 1.602*10**-19

def wavelength(energy_keV):
    """ Retuns the wavelength in metres  """
    return (h_planck * sci.speed_of_light) / (energy_keV * 1000.0)


def fft2c(data):
    return np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(data)))

def ifft2c(data):
    return np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(data)))
    
def intensity(data):
    return np.absolute(data) ** 2

def phase(data):
    return np.angle(data)

def fresnel_number(radius, distance, wavelength):
    fresnel_number = radius ** 2 / (distance * wavelength)
    return fresnel_number.simplified.magnitude
   
def none_to_zero_m(value):    
    if value == None:
        value = 0. * q.m
    return value

def screen(image, vmin=None, vmax=None, extent=None, plot_title = '', fontsize = 20):
    """ Shows the image 2D with a colorbar. you can set the minimum and maximum values vmin, vmax 
        to show in the colorbar. 
        
        vmin, vmax: min and max range of the colorbar
        extent: array containing the xmax, xmin, ymax, ymin
    """
    plt.figure()
    if extent == None:
        initial = plt.imshow(image, interpolation = 'none', vmin=vmin, vmax=vmax, origin='lower')
    else:
        initial = plt.imshow(image, interpolation = 'none', vmin=vmin, vmax=vmax, origin='lower'
                            , extent=extent)
    plt.title(plot_title, fontsize = fontsize)
    plt.colorbar(initial)


def show():
    plt.show()
    
def make_couple(data):
    """
    Convert input data to a tuple of two numbers.

    Parameters:
    - data (int/float/list/tuple): A single number or a list/tuple/array of two numbers.

    Returns:
    - tuple: A tuple of two numbers.

    If the input is a single number:
        - Returns a tuple with both elements equal to the input.
    If the input is a list/array/tuple of two numbers:
        - Returns a tuple of those two numbers.

    Raises:
    - ValueError: If input is not a single number or a list/array/tuple of exactly 2 numbers.
    """
    
    # If data is a single number
    if isinstance(data, (int, float)):
        return (data, data)
    
    # If data is a list, array, or tuple
    elif isinstance(data, (list, tuple)) and len(data) == 2:
        return tuple(data)
    
    else:
        raise ValueError("The input must be a single number or a list/array/tuple of exactly 2 numbers.")



def interp2(X, Y, V, Xq, Yq):
    """
    Returns interpolated values of a function of two variables at specific query points using linear 
    interpolation. The results always pass through the original sampling of the function. 
        X and Y contain the coordinates of the sample points. 
        V contains the corresponding function values at each sample point. 
        Xq and Yq contain the coordinates of the query points.    
    """
    
    ''' Definition of the interpolating function '''
    V_real = scipy.interpolate.RectBivariateSpline(X[1,:], Y[:,1], V.real)
    V_imag = scipy.interpolate.RectBivariateSpline(X[1,:], Y[:,1], V.imag)
    n = len(X[1,:])
    ''' Interpolation is made point-by-point, so matrices have to be unrolled in vectors and then 
        rolled back into matrices. '''   

    ''' The coordinate matrices are unrolled in vectors. '''
    Xq_unrol = Xq.reshape((1,n*n))
    Yq_unrol = Yq.reshape((1,n*n))

    ''' The interpolating functions are evaluated in the interpolated image. '''
    V_unrol_real = V_real.ev(Xq_unrol, Yq_unrol) 
    V_unrol_imag = V_imag.ev(Xq_unrol, Yq_unrol)
    V_unrol = V_unrol_real + 1j * V_unrol_imag

    ''' The coordinate matrices and the interpolated image are rolled back into matrices. '''
    Xq_rol = Xq_unrol.reshape((n,n))
    Yq_rol = Yq_unrol.reshape((n,n))
    V_temp = V_unrol.reshape((n,n))
    Vout = np.array(V_temp)    
    
    return Vout


def interp2python(X, Y, V, Xq, Yq, kind='cubic', fill_value=0.):
    """
    Returns interpolated values of a function of two variables at specific query points using linear 
    interpolation. The results always pass through the original sampling of the function. 
        X and Y contain the coordinates of the sample points. 
        V contains the corresponding function values at each sample point. 
        Xq and Yq contain the coordinates of the query points.    
    """
    
    ''' Definition of the interpolating function '''
    V_function_real = scipy.interpolate.interp2d(X, Y, V.real, kind=kind, fill_value=fill_value)
    V_interpolated_real = V_function_real(Xq[1,:], Yq[:,1])
    V_function_imag = scipy.interpolate.interp2d(X, Y, V.imag, kind=kind, fill_value=fill_value)
    V_interpolated_imag = V_function_imag(Xq[1,:], Yq[:,1])
    V_interpolated = V_interpolated_real + 1j * V_interpolated_imag
    
    return V_interpolated    


def wave_number(refractive_index, energy_eV):
    """
    Calculate the wave number of a material based on its refractive index and energy.

    Parameters:
    - refractive_index (float): Refractive index of the material.
    - energy_eV (float): Energy in electron volts.

    Returns:
    - float: Wave number of the material.
    """
    
    # Convert the energy from eV to keV for wavelength calculation.
    energy_keV = energy_eV / 1000
    
    # Calculate the wavelength of the material.
    wave_length = wavelength(energy_keV)
    
    # Calculate and return the wave number of the material.
    return 2 * np.pi * refractive_index / wave_length 


def transmitted_wavefield(thickness_image, refractive_index, energy_eV):
    """
    Calculate the transmitted wavefield of a material.

    Parameters:
    - thickness_image (2D numpy array): Image representing the thickness of the material.
    - refractive_index (float): Refractive index of the material.
    - energy_eV (float): Energy in electron volts.

    Returns:
    - 2D numpy array: Transmitted wavefield of the material.
    """

    # Calculate the wave number of the material.
    k_material = wave_number(refractive_index, energy_eV)
    
    # Calculate and return the transmitted wavefield of the material.
    return np.exp(1j * thickness_image * k_material)



def free_space_propagator(n, energy_eV, propagation_distance, 
                          pixel_size, oversampling=1, fresnel=True):
    """
    Computes the Fourier space propagator, either using the Fresnel approximation 
    or a more general form that takes into account wave_length.

    Parameters:
    - n (int): The size of the input square matrix (n x n).
    - wave_length (int or float): Photon energy of the wave in consideration, 
                                    in eV.
    - propagation_distance (float): The distance over which propagation 
                                    is considered, in meters.
    - pixel_size (float): The size of a pixel in spatial coordinates, in meters.
    - oversampling (int, optional): oversampling rate. It is the ratio between 
                                the size of the camera pixels over
                                the size of the pixels used during computation 
    - fresnel (bool, optional): Whether to use the Fresnel approximation. 
                                Default is True.

    Returns:
    - numpy.ndarray: The Fourier space propagator as a 2D matrix of shape (n, n).
    """
    
    # Calculate the wavelength from the photon energy
    energy_keV = energy_eV / 1000
    wave_length = wavelength(energy_keV)
    
    # Create frequency arrays for Fourier space
    freqs = np.fft.fftfreq(n*oversampling) / (pixel_size / oversampling)
    fx = np.tile(freqs, [n*oversampling, 1])
    fy = fx.transpose()

    # Check if we're using the Fresnel approximation
    if fresnel:
        propagator = np.exp(1j * 2 * np.pi / wave_length * propagation_distance) * \
                     np.exp(-1j * np.pi * wave_length * propagation_distance * (fx**2 + fy**2))
    else:
        propagator = np.exp(1j * 2 * np.pi / wave_length * propagation_distance * \
                            np.sqrt(1 - (fx * wave_length)**2 - (fy * wave_length)**2))

    return propagator

def make_gauss2D(n, sigma=1, pixel_size=1, x_centre=0, y_centre=0, oversampling=1):
    """
    Generate a 2D Gaussian image.

    Parameters:
    - n (int): Size of the resulting square image matrix (n x n).
    - sigma (float): Standard deviation of the Gaussian, in pixels.
    - pixel_size (float, optional): Size of a pixel in the resulting image,
                                    in meters . Default is 1.
    - x_centre, y_centre (float, optional): Center of the Gaussian 
                                            in pixel coordinates. 
                                            Default is (0, 0).
    - oversampling (int, optional): oversampling rate. It is the ratio between 
                                the size of the camera pixels over
                                the size of the pixels used during computation

    Returns:
    - image (numpy.ndarray): A 2D array representing the Gaussian image.
    - x, y (numpy.ndarray): Meshgrids representing the x and y pixel coordinates.
    """
    
    # Adjust for the oversampling
    n = n * oversampling
    x_centre = x_centre * oversampling
    y_centre = y_centre * oversampling
    sigma = sigma * oversampling
    pixel_size = pixel_size / oversampling
    
    # Initialize an empty image of the desired size
    image = np.zeros((n, n), dtype="float64")
    
    # Create meshgrid with coordinates adjusted by the center values
    y, x = np.mgrid[-n/2 - y_centre : n/2 - y_centre:1, 
                    -n/2 - x_centre : n/2 - x_centre:1]
    
    # Adjust the coordinates based on the pixel size
    x = (x + 0.5) * pixel_size
    y = (y + 0.5) * pixel_size
    sigma = sigma * pixel_size
    sigmaX = sigma 
    sigmaY = sigma 

    # Compute the 2D Gaussian values for each coordinate in the image
    # image = (1 / (2 * np.pi * sigmaX**2 * sigmaX**2)) * \
             # np.exp(- (x**2 / (2 * sigmaX**2) + y**2 / (2 * sigmaY**2)))
    image =  np.exp(- (x**2 / (2 * sigmaX**2) + y**2 / (2 * sigmaY**2)))
    image = image / np.sum(image)
    return image, x, y

def make_gauss2D_propagator(n, sigma=1, pixel_size=1, oversampling=1):
    """
    Generate the 2D Fourier transform of a Gaussian image.
    
    This function creates a 2D Gaussian image and then returns its 
    2D Fourier transform after shifting the zero frequency component to the center.
    
    Parameters:
    - n (int): Size of the resulting square image matrix (n x n).
    - sigma (float): Standard deviation of the Gaussian, in pixels.
    - pixel_size (float, optional): Size of a pixel in the resulting image,
                                    in meters. Default is 1.
    - oversampling (int, optional): oversampling rate. It is the ratio between 
                                the size of the camera pixels over
                                the size of the pixels used during computation

    Returns:
    - numpy.ndarray: Fourier transform of the 2D Gaussian image.
    """
    
    # Generate the Gaussian image using the provided make_gauss2D function
    gauss, _, _ = make_gauss2D(n, sigma=sigma, pixel_size=pixel_size, 
                               x_centre=0, y_centre=0, oversampling=oversampling)
    
    # Shift the zero frequency component to the center, take the 2D Fourier transform, and return
    return np.fft.fft2(np.fft.fftshift(gauss))


def section(img, center, size, vertical=False, cmap='gray', 
            vmin=None, vmax=None):
    """
    Plot a section of an image with a rectangle indicating the sliced region 
    and then plots the averaged line profile of that region.
    
    Parameters:
    - img (numpy array): 2D array representing the image.
    - center (tuple or array): (xc, yc) coordinate of the center of the section.
    - size (tuple or array): (xw, yh) width and height of the section.
    - yc (float): Center y-coordinate of the slice.
    - yh (float): Height of the slice.
    - xc (float): Center x-coordinate of the slice.
    - xw (float): Width of the slice.
    - vertical(bool, optional): compute a vertical section, default is horizontal.
    - cmap (string, optional): color map of the image
        
    Returns:
    None
    """
    # unpacking data about the section.
    xc, yc = center
    xw, yh = size
    
    # make sure that the coordinates are integers.
    xc = int(xc)
    yc = int(yc)
    xw = int(xw)
    yh = int(yh)
    
    # Display the provided image.
    plt.figure(figsize=(7,7))
    plt.imshow(img, cmap=cmap, vmin=vmin, vmax=vmax)
    
    # Draw a rectangle on the image to represent the section computed for the line plot.
    # Rectangle coordinates: (x_start, y_start), width, height
    rectangle = plt.Rectangle((xc-xw/2, yc-yh/2), xw, yh, edgecolor='red', facecolor='none')
    plt.gca().add_patch(rectangle)
    plt.show()
    
    # Extract the band from the image to compute the line profile.
    band = img[int(yc-yh/2):int(yc+yh/2), int(xc-xw/2):int(xc+xw/2)]
    
    # Compute the mean of the band along the y-axis to generate a line profile.
    if vertical is True:
        line = np.mean(band, axis=1)
    else:
        line = np.mean(band, axis=0)
    
    # Plot the line profile.
    plt.figure(figsize=(7,5))
    plt.plot(line)
    plt.show()



class Propagate:
    """
    A class for propagating an image wavefield.

    Attributes:
    - space_propagator: The propagator calculated based on input parameters.
    - propagated_wavefield: The wavefield after propagation.
    - propagated_intensity: The intensity of the propagated wavefield.
    """

    def __init__(self, img_wavefield, distance, energy_eV, pixel_size,
                 oversampling=10, fresnel=False, 
                 blur=True, sigma=None, downsample=True,
                 anti_aliasing=True):
        """
        Initialize the Propagate class.

        Parameters:
        - img_wavefield (ndarray): The input image wavefield.
        - distance (float): Propagation distance in meters.
        - energy_eV (float): Energy in electron volts.
        - pixel_size (float): The pixel size of the image.
        - oversampling (int, optional): Oversampling factor. Defaults to 10.
        - fresnel (bool, optional): If True, use Fresnel propagation. Defaults to False.
        - blur (bool, optional): If True, blur the image. Defaults to True.
        - sigma (float, optional): Sigma for the Gaussian blur. If None, calculated from oversampling.
        - downsample (bool, optional): If True, downsample the image. Defaults to True.
        - anti_aliasing (bool, optional): If True, apply anti-aliasing when downsampling. Defaults to True.
        """

        self.img_wavefield = img_wavefield
        self.oversampling = oversampling
        self.size = int(img_wavefield.shape[0] / oversampling)  # Calculate the size of the original image
        self.blur = blur
        self.sigma = sigma
        self.downsample = downsample
        self.anti_aliasing = anti_aliasing

        # Calculate the propagator for free space based on the given parameters
        self.space_propagator = free_space_propagator(self.size, energy_eV, 
                                                        distance, pixel_size, 
                                                        oversampling=self.oversampling,
                                                        fresnel=fresnel)

        self.propagated_wavefield = None
        self.propagated_intensity = None

    @staticmethod
    def intensity(wavefield):
        # Intensity is calculated as the absolute square of the wavefield
        return np.abs(wavefield)**2

    def get_propagated_wavefield(self):
        """Computes and returns the propagated wavefield."""
        result_fft = np.fft.fft2(self.img_wavefield) * self.space_propagator
        self.propagated_wavefield = np.fft.ifft2(result_fft)
        return self.propagated_wavefield

    def get_propagated_intensity(self):
        """Computes and returns the propagated intensity, with optional blurring and downsampling."""
        
        # Calculate the intensity of the propagated wavefield
        self.propagated_intensity = Propagate.intensity(self.get_propagated_wavefield())

        if self.blur:
            # Define the sigma for Gaussian blurring based on oversampling
            sigma = self.sigma if self.sigma is not None else self.oversampling / 10.0
            # Blur the propagated intensity
            self.propagated_intensity = gaussian_filter(self.propagated_intensity, sigma=sigma)
        
        if self.downsample:
            # Downsample the blurred image
            self.propagated_intensity = resize(self.propagated_intensity, (self.size, self.size), anti_aliasing=self.anti_aliasing)
        
        return self.propagated_intensity


pass