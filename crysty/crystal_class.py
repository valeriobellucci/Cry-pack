# -*- coding: utf-8 -*-
"""
Created on Tue Jun  5 09:36:45 2018

@author: Valerio
"""

from cry_functions import (compute_crystal_frequency_entrance, compute_crystal_frequency_exit,
                           bragg_angle, width_angle_o, magnification, asymmetry_best)
import numpy as np
import quantities as q
import pyopencl.array as cl_array
from utility import (intensity, make_couple)



class Crystal():
    """
    Class representing a crystal in Bragg diffraction.
    """
    
    def __init__(self, material, miller_indices, miscut, positioning_energy, **params):
        """
        INPUT:
            material: string representing the material, for now is 'Si' or 'Ge'.
            miller_indices: list representing the miller indices of the diffracting plane,
                            for example [2,2,0] for the Bragg Magnifier.
            miscut: miscut angle of the crystal with the quantity units, usually degrees or radians
                    --> multiply the number for q.deg or q.rad . The lattice plane parallel to the
                    physical surface has miscut zero (symmetric Bragg diffraction). If the miscut
                    makes the incidence angle of x-rays on the crystal more grazing, the miscut
                    angle is negative.If the miscut makes the incidence angle of x-rays on the
                    crystal less grazing, the miscut angle is positive.
            positioning_energy: central energy of the x-rays, for which the crystal diffracts
                    exactly at the Bragg Angle. This energy is used to define the angle between
                    beam and diffracting planes in space, so to position the crystal.
                    Energy should be a quantity object,usually expressed in keV or eV -->
                    multiply the number for q.keV or q.eV.
        """

        ''' Basic properties of the crystal. '''
        self.material = material  # 'Si' or 'Ge' for now
        self.miller_indices = miller_indices
        self.miscut = miscut
        self.asymmetry = - np.pi / 2 + self.miscut.rescale('rad').magnitude
        self.energy = positioning_energy
        self.energy_keV = positioning_energy.rescale('keV').magnitude
        positioning_energy = self.energy

        ''' Positioning of the crystal on the beam and some features of the beam. '''
        ''' misalignment angle of the crystal respect the Bragg angle, in arcsec
            (non a quantity object) '''
        self.misalignment_angle = params.get('misalignment_angle', 0.)
        ''' rocking angle of the crystal respect the Bragg angle,in arcsec
        (non a quantity object) '''
        self.rocking_angle = params.get('rocking_angle', 0.)
        '''Centering of the rocking curve: True|False
                                           True: the rocking curve is centred at the zero,
                                           False: the refraction shift is visible '''
        self.centred = params.get('centred', True)
        ''' I sthe crystal parallel or antiparallel to the previous crystal?
            True|False
            True: this crystal is in the parallel, non-dispersive geometry
            False: this crystal is in the antiparallel, dispersive geometry
        '''
        self.parallel = params.get('parallel', True)
        ''' Beam direction as seen by the crystal, 'entrance' or 'exit' '''
        self.beam_direction = params.get('beam_direction', 'entrance')
        ''' Diffraction direction, horizontal 'H' or vertical 'V' '''
        self.diffraction_direction = params.get('diffraction_direction', 'V')
 

    def description(self):
        """ Print a description of the crystal """
        desc_str = ("""The crystal is a %s %s with miscut of %s degrees. It is positioned at the Bragg
                    angle for energy %s keV"""
                    % (self.material, self.miller_indices, self.miscut.rescale('deg').magnitude,
                       self.energy_keV))
        print(desc_str)
    

    def _process_energy_arg(self, energy=None):
        # If no energy is provided, default to self.energy
        if not energy:
            return self.energy
    
        # Check if energy is a quantity object
        if not isinstance(energy, q.Quantity):
            raise ValueError("Provided energy is not a quantity object.")
    
        # Check if energy is a scalar (i.e., not an array or list)
        if hasattr(energy, "shape") and energy.shape:
            raise ValueError("Only a single value for energy is expected, not an array or list.")
    
        return energy


    def bragg_angle(self, energy = None):
        """
        Calculates the Bragg angle.

        INPUT:
            energy: energy of the x-rays, should be a quantity object, usually expressed in keV
                    or eV --> multiply the number for q.keV or q.eV. Default value is the central
                    energy at which the crystal is positioned self.energy
        RETURNS: Bragg angle in radians (quantity package).
        """
        energy = self._process_energy_arg(energy)
        return bragg_angle(
                            energy.rescale('keV').magnitude, 
                           self.material,
                           self.miller_indices
                           ) * q.rad

    def alfa(self, energy = None):
        """
        Calculates the incidence angle on the crystal.

        INPUT:
            energy: energy of the x-rays, should be a quantity object, usually expressed in keV
                    or eV --> multiply the number for q.keV or q.eV. Default value is the central
                    energy at which the crystal is positioned self.energy
        RETURNS: Incidence angle on the crystal in radians (quantity package).
        """
        energy = self._process_energy_arg(energy)
        return (self.bragg_angle(energy) + self.miscut).rescale('rad')

    def beta(self, energy = None):
        """
        Calculates the output angle coming out of the crystal.

        INPUT:
            energy: energy of the x-rays, should be a quantity object, usually expressed in keV
                    or eV --> multiply the number for q.keV or q.eV. Default value is the central
                    energy at which the crystal is positioned self.energy
        RETURNS: Output angle coming out of the crystal in radians (quantity package).
        """
        energy = self._process_energy_arg(energy)    
        return (self.bragg_angle(energy) - self.miscut).rescale('rad')        

    def fwhm_angle_entrance(self, energy = None):
        """
        Calculates the fwhm of the crystal rocking curve at a given energy.

        INPUT:
            energy: energy of the x-rays, should be a quantity object, usually expressed in keV
                    or eV --> multiply the number for q.keV or q.eV. Default value is the central
                    energy at which the crystal is positioned self.energy
        RETURNS: fwhm of the crystal rocking curve, in radians (quantity package).
        """
        energy = self._process_energy_arg(energy)    
        return width_angle_o(energy.rescale('keV').magnitude, self.asymmetry, self.material,
                             self.miller_indices) * q.rad

    def magnification(self, energy = None):
        """
        Calculates the magnification provided by the crystal at a given energy.

        INPUT:
            energy: energy of the x-rays, should be a quantity object, usually expressed in keV
                    or eV --> multiply the number for q.keV or q.eV. Default value is the central
                    energy at which the crystal is positioned self.energy
        RETURNS: magnification, float number.
        """
        energy = self._process_energy_arg(energy)    
        return magnification(energy.rescale('keV').magnitude, self.asymmetry, self.material,
                             self.miller_indices)

    def asymmetry_best(self, energy = None):
        """
        Calculates the best miscut that the crystal should have to maximize the resolution at
        a given energy.

        INPUT:
            energy: energy of the x-rays, should be a quantity object, usually expressed in keV or
                    eV --> multiply the number for q.keV or q.eV. Default value is the central
                    energy at which the crystal is positioned self.energy
        RETURNS: miscut angle that maximizes the resolution, in radians (quantity package).
        """
        energy = self._process_energy_arg(energy)    
        return asymmetry_best(energy.rescale('keV').magnitude, self.material,
                              self.miller_indices) * q.rad

    def _transfer_fourier(self, shape, pixel_size, energy, t=None, queue=None, out=None,
                          block=False):
        """
        Compute the crystalline transfer function, WITH fourier shift.

        INPUT:
            shape: Should contain the pixels number. Can be an integer for square shape, or a tuple
                   of two integers for rectangular shape --> example: n for square or (nx, ny) for
                   rectangle.
            pixel_size: Pixel size in quantity units, usually um or m --> multiply the number for
                        q.um or q.m . Can be an integer for square pixels, or a tuple of two
                        integersfor rectangular pixels --> example: p for square pixels or (px, py)
                        forrectangular pixels.
            energy: energy of the x-ray beam. Can be equal to the positioning energy self.energy
                    or not. Should be a quantity object, usually expressed in keV or eV -->
                    multiply the number for q.keV or q.eV.
        RETURNS : Crystal frequency transfer
        """

        # energy corresponding to the Bragg angle
        positioning_energy_keV = self.energy.rescale('keV').magnitude
        # energy of the incident x-rays
        experiment_energy_keV = energy.rescale('keV').magnitude
        # difference between positioning and experiment energy
        shift_poly_energy = positioning_energy_keV - experiment_energy_keV

        (nx, ny) = make_couple(shape)  # number of pixels for x and y direction
        (px, py) = make_couple(pixel_size)  # size of pixels for x and y direction

        # spatial frequencies over x
        qx_max = (1. / (2. * px)).simplified.magnitude
        qx_step = (1. / (2. * px * nx)).simplified.magnitude
        qx = np.arange(- qx_max / 2., + qx_max / 2., qx_step)

        # spatial frequencies over y
        qy_max = (1. / (2. * py)).simplified.magnitude
        qy_step = (1. / (2. * py * ny)).simplified.magnitude
        qy = np.arange(- qy_max / 2., + qy_max / 2., qy_step)
        QX, QY = np.meshgrid(qx, qy)

        ''' choose a diffraction direction, horizontal 'H' or vertical 'V' '''
        if self.diffraction_direction == 'H':
            QX_or_QY = QX
        elif self.diffraction_direction == 'V':
            QX_or_QY = QY
        else:
            print('''Choose diffraction direction, horizontal 'H' or vertical 'V' ''')

        ''' choose a beam direction as seen by the crystal, 'entrance' or 'exit' '''
        if self.beam_direction == 'entrance':
            crystal_frequency_transfer = \
                compute_crystal_frequency_entrance(QX_or_QY,
                                                    self.misalignment_angle,
                                                    self.asymmetry,
                                                    positioning_energy_keV,
                                                    shift_poly_energy,
                                                    self.material,
                                                    self.miller_indices,
                                                    self.centred,
                                                    self.parallel,
                                                    1.
                                                    )[0]
        elif self.beam_direction == 'exit':
            crystal_frequency_transfer = \
                compute_crystal_frequency_exit(QX_or_QY,
                                                self.misalignment_angle,
                                                self.asymmetry,
                                                positioning_energy_keV,
                                                shift_poly_energy,
                                                self.material,
                                                self.miller_indices,
                                                self.centred,
                                                self.parallel,
                                                1.
                                                )[0]
        else:
            print('''Choose beam direction, 'entrance' or 'exit' ''')

        # Fix the singularity issue
        reflectivity = intensity(crystal_frequency_transfer)
        i = 0
        while np.any((reflectivity > 1.) | (reflectivity.round(3) == 0.)):    
            print('Singularity issue found in crystalline transfer function. Try to solve: ', i)
            size = reflectivity.shape[0] 
            correction = np.where((reflectivity > 1.) | (reflectivity.round(3) == 0.)) 
            crystal_frequency_transfer[(reflectivity > 1.) | (reflectivity.round(3) == 0.)] \
                = crystal_frequency_transfer[correction[0][size/2] - i, correction[1][size/2] - i]
                
            reflectivity = intensity(crystal_frequency_transfer)
            i += 1

        crystal_frequency_transfer = np.fft.fftshift(crystal_frequency_transfer)
        crystal_frequency_transfer = crystal_frequency_transfer.astype(np.complex128)

        return cl_array.to_device(queue, crystal_frequency_transfer)

    def crystal_transfer_fourier(self, shape, pixel_size, energy, 
                                 t=None, queue=None, out=None,
                                 block=False):
        """
        Compute the crystalline transfer function, WITHOUT fourier shift.

        INPUT:
            shape: Should contain the pixels number. Can be an integer for square shape, or a tuple
                   of two integers for rectangular shape --> example: n for square or (nx, ny) for
                   rectangle.
            pixel_size: Pixel size in quantity units, usually um or m --> multiply the number for
                        q.um or q.m . Can be an integer for square pixels, or a tuple of two
                        integersfor rectangular pixels --> example: p for square pixels or (px, py)
                        forrectangular pixels.
            energy: energy of the x-ray beam. Can be equal to the positioning energy self.energy
                    or not. Should be a quantity object, usually expressed in keV or eV -->
                    multiply the number for q.keV or q.eV.
        RETURNS : Crystal frequency transfer
        """

        return np.fft.ifftshift(self.transfer_fourier(shape, pixel_size, energy, t=t, queue=queue,
                                                      out=out, block=block).get())

    def QX_or_QY_transfer_fourier(self, QX_or_QY, energy, singularity_shift=0.01):
        """
        Compute the crystalline transfer functions.

        input:
        QX_or_QY : mesh of frequencies over X or Y direction. It should be a quantity number.
                    (depending whether the magnification is in the X or Y directio)
        energy: energy of the x-ray beam. Can be equal to the positioning energy self.energy
                    or not. Should be a quantity object, usually expressed in keV or eV -->
                    multiply the number for q.keV or q.eV. 
        singularity_shift = shift in the rocking angle to exclude the singularity in zero. This shift is
                        expressed in units of the Darwin width (half of the FWHM).
        
        returns: Crystal frequency transfer
        """

        # energy corresponding to the Bragg angle
        positioning_energy_keV = self.energy.rescale('keV').magnitude
        # energy of the incident x-rays
        experiment_energy_keV = energy.rescale('keV').magnitude
        # difference between positioning and experiment energy
        shift_poly_energy = positioning_energy_keV - experiment_energy_keV

        # Rescales spatial frequency in (1/m) and removes dimension
        QX_or_QY = QX_or_QY.simplified.magnitude

        if self.beam_direction == 'entrance':
            crystal_frequency_transfer = \
                compute_crystal_frequency_entrance(QX_or_QY,
                                                    self.misalignment_angle,
                                                    self.asymmetry,
                                                    positioning_energy_keV,
                                                    shift_poly_energy,
                                                    self.material,
                                                    self.miller_indices,
                                                    self.centred,
                                                    self.parallel,
                                                    1.,
                                                    singularity_shift=
                                                    singularity_shift
                                                    )[0]
        elif self.beam_direction == 'exit':
            crystal_frequency_transfer = \
                compute_crystal_frequency_exit(QX_or_QY,
                                                self.misalignment_angle,
                                                self.asymmetry,
                                                positioning_energy_keV,
                                                shift_poly_energy,
                                                self.material,
                                                self.miller_indices,
                                                self.centred,
                                                self.parallel,
                                                1.,
                                                singularity_shift=
                                                singularity_shift
                                                )[0]
        else: 
            print ('Choose direction')
      
        # Fix the singularity issue
        reflectivity = intensity(crystal_frequency_transfer)
        i = 0
        while np.any((reflectivity > 1.) | (reflectivity.round(3) == 0.)):    
            print ('Singularity issue found in crystalline transfer function. Try to solve: ', i)
            size = reflectivity.shape[0] 
            correction = np.where((reflectivity > 1.) | (reflectivity.round(3) == 0.)) 
            crystal_frequency_transfer[(reflectivity > 1.) | (reflectivity.round(3) == 0.)] \
                = crystal_frequency_transfer[correction[0][size/2] - i, correction[1][size/2] - i]
                
            reflectivity = intensity(crystal_frequency_transfer)
            i += 1

        return crystal_frequency_transfer



pass