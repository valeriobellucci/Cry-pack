# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 10:41:13 2018

@author: xx9905
"""

from syris.geometry import Trajectory
from syris.devices.sources import *
from bragg_magnifier_class import Bragg_Magnifier
from monochromator_class import Monochromator
from crystal_class import Crystal
from utility import screen
from construction import *
from bragg_magnifier_class import Bragg_Magnifier
from monochromator_class import Monochromator
from syris.bodies.simple import *
from newe_objects import make_4_spheres, make_sphere_generic, make_disk_generic, make_4_disks, make_rectangle



def make_triangle(n=128):
    lin = np.linspace(0, 2, n, endpoint=False)
    y = np.abs(lin - 1)
    z = np.zeros(n)

    return zip(z, y, z) * q.mm


' I need a source. '
class Source:
    """ Source. 
        Makes the source image at the object position.
    
    Choose a source.
    Can be 'bending_magnet' or 'wiggler'. Example: Source('wiggler') 
    """
    
    from syris.devices.sources import BendingMagnet, Wiggler 
    syris.init()
    
    def __init__(self, type_of_source, where = 'before'):
        
        if (where == 'before') | (where == True):
            self.ps_now = ps
        elif (where == 'no_BM') | (where == False):
            self.ps_now = ps_exit
        else:
            print 'Choose a position for the Source. Before the BM or without the BM.'          
        
        self.pixel_size = self.ps_now
        self.trajectory = Trajectory([(n / 2, n / 2, 0)] * self.pixel_size) # Position of the centre of the source
        
        if type_of_source == 'bending_magnet':
            self.source = BendingMagnet(electron_energy, el_current, magnetic_field, 
                                        sample_distance, dE, size, self.pixel_size, 
                                        self.trajectory)
        elif type_of_source == 'wiggler':                           
            self.source = Wiggler(electron_energy, el_current, magnetic_field, 
                                  sample_distance, dE, size, self.pixel_size, 
                                  self.trajectory, num_periods) 
        else:
            print 'Choose a source type'   
        self.type_of_source = type_of_source
                             
    
    def transfer(self, energy):
        return self.source.transfer((n, n), self.pixel_size, energy).get()


    def apply_blur(self, intensity, distance, pixel_size, queue=None, block=False):
        """Apply source blur based on van Cittert-Zernike theorem at *distance*."""
        return self.source.apply_blur(intensity, distance, pixel_size, queue=None, block=False)


class Monochromator_simple:
    """ Simplified monochromator, double Crystals in non-dispersive configuration
    
    """
    def __init__(self, shift='vertical', where = 'before'):
        
        ### Old code, possibly useless, it used to work with the below QY / self.M_now and QX / self.M_now
        if (where == 'before') | (where == True):
            self.M_now = M
        elif (where == 'after') | (where == False) | (where == 'no_BM'):
            self.M_now = 1
        else:
            print 'Choose a position for the Monochromator. Before of after the BM, or without the BM.'
            
        self.shift = shift
        if shift == 'vertical':
            self.Qnow = QY_exit * self.M_now #QY / self.M_now
        elif shift == 'horizontal':
            self.Qnow = QX_exit * self.M_now #QX / self.M_now
            
        self.crystal_1 = Crystal(material_mono, miller_index_mono, miscut_mono[0], central_energy)
        self.crystal_2 = Crystal(material_mono, miller_index_mono, miscut_mono[1], central_energy)
        
    def transfer(self, energy):
        """ Crystalline transfer function """
        frequency_transfer_1 = self.crystal_1.QX_or_QY_transfer_fourier(Qnow, energy)
        frequency_transfer_2 = self.crystal_2.QX_or_QY_transfer_fourier(Qnow, energy) 
        frequency_transfer_12 = frequency_transfer_1 * frequency_transfer_2
        return frequency_transfer_12
    
    def crystals_reflectivity(self, energy , crystals=0):
        """ Reflectivity of the monochromator """
        if crystals == 0:
            ref = intensity(self.crystal_1.QX_or_QY_transfer_fourier(Qnow, energy)
                            * self.crystal_2.QX_or_QY_transfer_fourier(Qnow, energy) )       
        elif crystals == 1:
            ref = intensity(self.crystal_1.QX_or_QY_transfer_fourier(Qnow, energy)) 
        elif crystals == 2:
            ref = intensity(self.crystal_2.QX_or_QY_transfer_fourier(Qnow, energy) ) 
        else:
            print 'Choose a correct crystal, nothing or 0 for both crystals, 1 for crystal 1 and 2 for crystal 2.'
        return ref


    def crystals_phaseshift(self, energy, crystals=0):
        """ Phase shift produced by the monochromator """
        if crystals == 0:
            pha = phase(self.crystal_1.QX_or_QY_transfer_fourier(Qnow, energy)
                        * self.crystal_2.QX_or_QY_transfer_fourier(Qnow, energy) )       
        elif crystals == 1:
            pha = phase(self.crystal_1.QX_or_QY_transfer_fourier(Qnow, energy)) 
        elif crystals == 2:
            pha = phase(self.crystal_2.QX_or_QY_transfer_fourier(Qnow, energy) ) 
        else:
            print 'Choose a correct crystal, nothing or 0 for both crystals, 1 for crystal 1 and 2 for crystal 2.'
        return pha
        

    def propagate(self, input_image, energy):
        """ Takes an input_image in real space, propagates it in fourier space through the monochromator, 
            and returns the output_image in real space.
        """        
        input_image_fs = fft2c(input_image) 
        output_image_fs = input_image_fs * self.transfer(energy)
        output_image = ifft2c(output_image_fs)
        return output_image
    
    def intensity(self, input_image, energy):
        """ Takes an input_image in real space, propagates it in fourier space through the monochromator, 
            and returns the intensity_image in real space.
        """ 
        return intensity(self.propagate(input_image, energy))
    
    
    def phase(self, input_image):
        """ Takes an input_image in real space, propagates it in fourier space through the monochromator, 
            and returns the phase_image in real space.
        """ 
        return phase(self.propagate(input_image, energy))
    

# The monochromator transfer function is already in fourier space   
#mono_0_fs = Monochromator_2C.transfer()

# Turn the monochromator transfer function in real space
#mono_0 = ifft2c(mono_0_fs)

' I need air (or helium ??) '
class Space:
    """
    Create space where the wave propagates.
    
    Inputs: 
        space --> has to be a quantity input, in meters, mm or etc...    
    
    Possible computation ways are 'syris' or 'analytic'.
    propagator methods produce the propagators before or after the BM.
    propagate method take an input_image in real space, propagates it in fourier space for the distance, 
            and returns the output_image in real space.
    """
    syris.init()
    
    def __init__(self, space):
        self.space = space
    
    def propagator_before_BM(self, energy, computation='syris'):
        """ Returns the propagator if the space is positioned before the Bragg Magnifier """
        if computation == 'syris':
            propagator = np.fft.ifftshift( compute_propagator(n, self.space,
                                                              energy_to_wavelength(energy), ps,
                                                                mollified=mollified).get() )
        elif computation == 'analytic':
            propagator = np.exp(-1j * (2 * np.pi) * (self.space / (2. / energy_to_wavelength(energy))
                                * (QX ** 2 + QY ** 2) ).simplified.magnitude) 
        else:
            print 'Choose a propagation method'
        return propagator
        
    def propagator_after_BM(self, energy, computation='syris'):
        """ Returns the propagator if the space is positioned after the Bragg Magnifier """
        if computation == 'syris':
            propagator = np.fft.ifftshift( compute_propagator(n, self.space,
                                                              energy_to_wavelength(energy), ps_exit, 
                                                              mollified=mollified).get() )
        elif computation == 'analytic':
            propagator = np.exp(-1j * (2*np.pi) * (self.space / (2. / energy_to_wavelength(energy))
                                * (QX_exit ** 2 + QY_exit ** 2) ).simplified.magnitude) 
        else:
            print 'Choose a propagation method'
        return propagator
        
    def propagate_before_BM(self, input_image, energy, computation='syris'):
        """ Takes an input_image in real space, propagates it in fourier space for the distance, 
            and returns the output_image in real space.
            
            The distance is positioned before the Bragg Magnifier.
        """        
        input_image_fs = fft2c(input_image) 
        output_image_fs = input_image_fs * self.propagator_before_BM(energy, computation=computation)
        output_image = ifft2c(output_image_fs)
        return output_image
    
    def propagate_after_BM(self, input_image, energy, computation='syris'):
        """ Takes an input_image in real space, propagates it in fourier space for the distance, 
            and returns the output_image in real space.
            
            The distance is positioned after the Bragg Magnifier.
        """        
        input_image_fs = fft2c(input_image) 
        output_image_fs = input_image_fs * self.propagator_after_BM(energy, computation=computation)
        output_image = ifft2c(output_image_fs)
        return output_image



' I need a sample. '
class Sample:
    """ Sample. 
        Makes the object image.
    
    Choose a sample. Can be:
        'sphere_scaled' : uniform plexiglass sphere with radius rescaled for the FoV. 
            It has attribute self.ratio --> the radius of the sphere is 1/ratio.
            Default ratio is 4, which means the radius of the sphere is 1/4 of the FoV.
        'sphere_fix : uniform plexiglass sphere with fixed radius.
            Default radius is 50 * q.um .
    """
    syris.init()
    
    def __init__(self, sample_name, where = 'before', ratio = ratio, radius = radius,
                 thickness = thickness, x_centre=x_centre, y_centre=y_centre, xdil=xdil, ydil=ydil,
                 separation=separation, material=sample_material):
        
        if (where == 'before') | (where == True):
            self.ps_now = ps
        elif (where == 'after') | (where == False) | (where == 'no_BM'):
            self.ps_now = ps_exit
        else:
            print 'Choose a position for the Sample. Before of after the BM, or without the BM.'        
        
        if library_sample == False:
            if type(sample) == object:
                self.sample = sample
            else:
                raise ValueError('Sample should be an object')
        elif library_sample == True:
            if sample_name == 'sphere_scaled':        
                self.ratio = ratio
                self.radius = n / self.ratio * self.ps_now
                self.sample = make_sphere(n, self.radius, pixel_size=self.ps_now, material=material)
                
            elif sample_name == '4_spheres_scaled':        
                self.ratio = ratio * 2.
                self.separation = separation
                self.radius = n / self.ratio * self.ps_now
                self.sample = make_4_spheres(n, self.radius, separation = self.separation, 
                                             pixel_size=self.ps_now, material=material)
                                             
            elif sample_name == 'sphere_fix':  
                self.radius = radius
                self.sample = make_sphere(n, self.radius, pixel_size=self.ps_now, material=material)
                
            elif sample_name == 'sphere_generic':  
                self.radius = radius
                self.x_centre = x_centre
                self.y_centre = y_centre
                self.xdil=xdil
                self.ydil=ydil
                self.sample = make_sphere_generic(n, self.radius, pixel_size=self.ps_now,
                                                  material=material, x_centre=self.x_centre,
                                                  y_centre=self.y_centre, xdil=self.xdil, ydil=self.ydil)
                                                  
            elif sample_name == 'disk_generic':  
                self.radius = radius
                self.thickness = thickness
                self.x_centre = x_centre
                self.y_centre = y_centre
                self.xdil=xdil
                self.ydil=ydil
                self.sample = make_disk_generic(n, self.radius, self.thickness, pixel_size=self.ps_now,
                                                material=material, x_centre=self.x_centre,
                                                y_centre=self.y_centre, xdil=self.xdil, ydil=self.ydil)
                                                
            elif sample_name == '4_disks_scaled':        
                self.ratio = ratio * 2.
                self.separation = separation
                self.thickness = thickness
                self.radius = n / self.ratio * self.ps_now
                self.sample = make_4_disks(n, self.radius, self.thickness, separation = self.separation, 
                                             pixel_size=self.ps_now, material=material)
                                         
            else:
                raise ValueError('Choose a sample')
        else:
            raise ValueError('library_sample should be True or False')
        self.sample_name = sample_name    
    
    def transfer(self, energy):
        if self.sample_name == 'rectangle':
            sample_transfer = self.sample
        else:
            sample_transfer = self.sample.transfer(shape, self.ps_now, energy).get()
        return sample_transfer
    
    def transfer_before_BM(self, energy):
        if self.sample_name == 'rectangle':
            sample_transfer = self.sample
        else:
            sample_transfer = self.sample.transfer(shape, ps, energy).get()
        return sample_transfer
        
    def transfer_after_BM(self, energy):
        if self.sample_name == 'rectangle':
            sample_transfer = self.sample
        else:
            sample_transfer = self.sample.transfer(shape, ps_exit, energy).get()
        return sample_transfer


class BraggMagnifier:
    """Propagates the object image through the Bragg magnifier
    
    Possible propagation methods are:
        'ED' --> effective distance method. 
        'RCT' --> reciprocal coordinate transformation method.
        'Syris' --> effective distance method using the propagators from Syris
    """
    
    def __init__(self):
    
    
    def propagation(self, fourier_wavefield, energy, crystalline_transfer_functions=True,
                    magnify=True, mollified=False, method='ED'):
        ''' Returns an object containing all the quantities of the propagation. '''
        propagation = Bragg_Magnifier(crystals_BM, distances_BM).propagate_magnifier_Q(fourier_wavefield,
                                            Q_magn, energy,
                                            crystalline_transfer_functions=
                                            crystalline_transfer_functions,
                                            magnify=magnify, mollified=mollified, method=method)
        return propagation

    def wavefield(self, fourier_wavefield, energy, crystalline_transfer_functions=True,
                    magnify=True, mollified=False, method='ED'):
        ''' Exit wavefield in Fourier space. '''  
        propagation = self.propagation(fourier_wavefield, energy,
                                       crystalline_transfer_functions=
                                       crystalline_transfer_functions,
                                       magnify=magnify, mollified=mollified, method=method) 
        return propagation.wavefield  
    

class Monochromator_now:
    """Propagates the object image through the Monochromator
    
    Possible propagation methods are:
        'ED' --> effective distance method. 
        'RCT' --> reciprocal coordinate transformation method.
        'Syris' --> effective distance method using the propagators from Syris
    """
    
    def __init__(self):
    
    
    def propagation(self, fourier_wavefield, energy, crystalline_transfer_functions=True,
                    magnify=True, mollified=False, method='ED'):
        ''' Returns an object containing all the quantities of the propagation. '''
        propagation = Monochromator(crystals_mono, distances_mono).propagate_magnifier_Q(fourier_wavefield,
                                            Q_magn, energy,
                                            crystalline_transfer_functions=
                                            crystalline_transfer_functions,
                                            magnify=magnify, mollified=mollified, method=method)
        return propagation

    def wavefield(self, fourier_wavefield, energy, crystalline_transfer_functions=True,
                    magnify=True, mollified=False, method='ED'):
        ''' Exit wavefield in Fourier space. '''  
        propagation = self.propagation(fourier_wavefield, energy,
                                       crystalline_transfer_functions=
                                       crystalline_transfer_functions,
                                       magnify=magnify, mollified=mollified, method=method) 
        return propagation.wavefield 

    def transfer(self, Q_in, energy, crystalline_transfer_functions=True,
                                  magnify=True, mollified=False, method='ED'):
        return propagator_Q(self, Q_in, energy, crystalline_transfer_functions=
                            crystalline_transfer_functions, magnify=magnify, mollified=mollified,
                            method=method).propagator


class Medipix:
    """ Temptative Medipix detector class. """
    from syris.devices.cameras import Camera
    syris.init()
    
    def __init__(self):    
        self.pixel_size = 55 * q.um
        self.gain =  1 # Med = 1, Example = 0.1
        self.dark_current =  0 # Med = 0, Example = 500
        self.amplifier_sigma =  0 # Med = 0 ??, Example = 23
        self.bits_per_pixel =  32 # Med = 12 or 24 , Example = 32
        self.shape = (512, 512) # Med = (256, 256) single, Hexa = 2x3 module = (768, 512) , Example = (n, n)

        
    def camera(self, exp_time=1 * q.s):
        camera_now = syris.devices.cameras.Camera(self.pixel_size, self.gain, self.dark_current, 
                                              self.amplifier_sigma, self.bits_per_pixel, self.shape
                                              , exp_time = exp_time, dtype=np.float32)
        return camera_now                                    
    
    def get_image(self, intensity_image, exp_time=1 * q.s):       
        camera_image = self.camera(exp_time = exp_time).get_image(intensity_image)
        return camera_image
    
    def print_image(self, intensity_image, exp_time=1 * q.s):
        camera_image = self.get_image(intensity_image, exp_time = exp_time)
        print camera_image
    
    def show_image(self, intensity_image, exp_time=1 * q.s):     
        camera_image = self.get_image(intensity_image, exp_time = exp_time)
        plt.figure()
        initial = plt.imshow(camera_image, interpolation = 'none')
        plt.colorbar(initial)
        plt.show()

    def phantom_show_image(self, intensity_image, exp_time=1 * q.s):     
        camera_image = self.get_image(intensity_image, exp_time = exp_time)
        plt.figure()
        initial = plt.imshow(camera_image, interpolation = 'none')
        plt.colorbar(initial)
        

class MedipixHexaGaAs:
    """ Medipix Hexa made of GaAs 500 um thick, shape = (512, 768) pixels """
    
    def camera(self, exp_time=1 * q.s):
        from pcd import TimepixHexa
        camera_now = TimepixHexa(exp_time=exp_time)
        return camera_now                                    
    
    def get_image(self, intensity_image, wavelength, photons_pixel_size, exp_time=1 * q.s): 
        from syris.imageprocessing import crop
        x0 = round(intensity_image.shape[0] / 2
                    - (28.16*q.mm/(photons_pixel_size)).simplified.magnitude / 2)
        y0 = round(intensity_image.shape[1] / 2
                    - (42.24*q.mm/(photons_pixel_size)).simplified.magnitude / 2)
        x1 = round(intensity_image.shape[0] / 2
                    + (28.16*q.mm/(photons_pixel_size)).simplified.magnitude / 2)
        y1 = round(intensity_image.shape[1] / 2
                    + (42.24*q.mm/(photons_pixel_size)).simplified.magnitude / 2)         
        photons = intensity_image[x1:x0:-1,y0:y1]
        photons = photons * exp_time.simplified.magnitude
        camera_image = self.camera(exp_time = exp_time).get_image(photons, wavelength, photons_pixel_size)
        return camera_image    
    
    def print_image(self, intensity_image, wavelength, photons_pixel_size, exp_time=1 * q.s):
        camera_image = self.get_image(intensity_image, wavelength, photons_pixel_size, exp_time = exp_time)
        print camera_image
    
    def show_image(self, intensity_image, wavelength, photons_pixel_size, exp_time=1 * q.s):     
        camera_image = self.get_image(intensity_image, wavelength, photons_pixel_size, exp_time = exp_time)
        plt.figure()
        initial = plt.imshow(camera_image, interpolation = 'none')
        plt.colorbar(initial)
        plt.show()

    def phantom_show_image(self, intensity_image, wavelength, photons_pixel_size, exp_time=1 * q.s):     
        camera_image = self.get_image(intensity_image, wavelength, photons_pixel_size, exp_time = exp_time)
        plt.figure()
        initial = plt.imshow(camera_image, interpolation = 'none')
        plt.colorbar(initial)


