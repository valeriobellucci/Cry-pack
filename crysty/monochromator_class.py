# -*- coding: utf-8 -*-
"""
Created on Mon Jul  2 12:28:25 2018

@author: xx9905
"""

import syris
import pyopencl.array as cl_array
from monochromator_class import Monochromator
from crystal_class import Crystal
from bragg_magnifier_class import Bragg_Magnifier


' I need a monochromator. '

class Monochromator(Bragg_Magnifier):
    """ Monochromator double Crystals in non-dispersive configuration
    
    """

    def propagate_Q(self, D0, Q_in, energy, crystalline_transfer_functions=True,
                                  magnify=True, mollified=False, method='ED'):
        """ Propagate through space and the crystals of the Monochromator with the different methods. 
        
        Input: (D0, Q, energy)
        D0 : Fourier transform of the input wavefield at the source position. It has to be an nxn matrix.
        Q_in : has to be a tuple (QX, QY). Reciprocal space coordinates for X and Y direction. 
            Each QX and QY has to be an nxn matrix of quantity numbers.
        energy: energy of the x-ray beam. Can be equal to the positioning energy self.energy
            or not. Should be a quantity object, usually expressed in keV or eV --> multiply the
            number for q.keV or q.eV.
        crystalline_transfer_functions = True / [True,True,True,True, ...] for the number of
            crystals composing the Bragg Magnifier. If crystalline_transfer_functions=False,
            all the crystal functions are removed from propagation. If one of the entry of
            crystalline_transfer_functions is False, that crystal function is removed from
            propagation. If one of the entry is 0, that 0 means False. If one of the entry is 1,
            that 1 means True.
        magnify = True. If false, the magnification is not applied.
        mollified = False or True. Apply mollification to the propagators.
        method = 'ED', 'RCT', 'Syris', 'ED Stano' -> use the Effective Distance method (ED),
            or Reciprocal Space Transformation method (RCT), using the Syris propagators or the
            Effective Distance Method from the Stano's notebooks.
        
        Output: object
        .wavefield : Fourier transform of the otput wavefield at the detector position 
            (or any arbitrary point after the magnifier). It is an nxn matrix.
        .Qout : tuple (QX_out, QY_out). Reciprocal space coordinates for X and Y direction after
                magnification. It is an nxn matrix of quantity numbers.
        .crystal_frequency_transfer : crystal frequency transfer of the crystal.
        .Q_s : array composed of arrays (QX_out, QY_out) of all the reciprocal space coordinates used for propagation.
                The number of coordinates is equal to the number of crystals + 1.
        .propagator : total propagator.
        """
        
        return self.propagate_magnifier_Q(D0, Q_in, energy, crystalline_transfer_functions=
                                            crystalline_transfer_functions, magnify=magnify,
                                            mollified=mollified, method=method)                          


    def propagator_Q(self, Q_in, energy, crystalline_transfer_functions=True,
                                  magnify=True, mollified=False, method='ED'):
        """ Propagator through space and the crystals of the Monochromator with the different methods. 
        
        Input: (Q, energy)
        Q_in : has to be a tuple (QX, QY). Reciprocal space coordinates for X and Y direction. 
            Each QX and QY has to be an nxn matrix of quantity numbers.
        energy: energy of the x-ray beam. Can be equal to the positioning energy self.energy
            or not. Should be a quantity object, usually expressed in keV or eV --> multiply the
            number for q.keV or q.eV.
        crystalline_transfer_functions = True / [True,True,True,True, ...] for the number of
            crystals composing the Bragg Magnifier. If crystalline_transfer_functions=False,
            all the crystal functions are removed from propagation. If one of the entry of
            crystalline_transfer_functions is False, that crystal function is removed from
            propagation. If one of the entry is 0, that 0 means False. If one of the entry is 1,
            that 1 means True.
        magnify = True. If false, the magnification is not applied.
        mollified = False or True. Apply mollification to the propagators.
        method = 'ED', 'RCT', 'Syris', 'ED Stano' -> use the Effective Distance method (ED),
            or Reciprocal Space Transformation method (RCT), using the Syris propagators or the
            Effective Distance Method from the Stano's notebooks.
        
        Output: object
        .propagator : total propagator.
        """                              
        return self.propagator_magnifier_Q(Q_in, energy, crystalline_transfer_functions=
                                            crystalline_transfer_functions, magnify=magnify,
                                            mollified=mollified, method=method).propagator                              


    def propagate_ps(self, D0, ps_in, energy, crystalline_transfer_functions=True,
                                  magnify=True, mollified=False, method='ED'):
        """ Propagate through space and the crystals of the Monochromator with the different methods. 
        
        Input: (D0, Q, energy)
        D0 : Fourier transform of the input wavefield at the source position. It has to be an nxn matrix.
        ps_in : entry pixel size. Has to de a quantity value in the form ps_in = make_tuple(ps)
        energy: energy of the x-ray beam. Can be equal to the positioning energy self.energy
            or not. Should be a quantity object, usually expressed in keV or eV --> multiply the
            number for q.keV or q.eV.
        crystalline_transfer_functions = True / [True,True,True,True, ...] for the number of
            crystals composing the Bragg Magnifier. If crystalline_transfer_functions=False,
            all the crystal functions are removed from propagation. If one of the entry of
            crystalline_transfer_functions is False, that crystal function is removed from
            propagation. If one of the entry is 0, that 0 means False. If one of the entry is 1,
            that 1 means True.
        magnify = True. If false, the magnification is not applied.
        mollified = False or True. Apply mollification to the propagators.
        method = 'ED', 'RCT', 'Syris', 'ED Stano' -> use the Effective Distance method (ED),
            or Reciprocal Space Transformation method (RCT), using the Syris propagators or the
            Effective Distance Method from the Stano's notebooks.
        
        Output: object
        .wavefield : Fourier transform of the otput wavefield at the detector position 
            (or any arbitrary point after the magnifier). It is an nxn matrix.
        .Qout : tuple (QX_out, QY_out). Reciprocal space coordinates for X and Y direction after
                magnification. It is an nxn matrix of quantity numbers.
        .crystal_frequency_transfer : crystal frequency transfer of the crystal.
        .Q_s : array composed of arrays (QX_out, QY_out) of all the reciprocal space coordinates used for propagation.
                The number of coordinates is equal to the number of crystals + 1.
        .propagator : total propagator.
        """
        
        return self.propagate_magnifier_ps(D0, ps_in, energy, crystalline_transfer_functions=
                                            crystalline_transfer_functions, magnify=magnify,
                                            mollified=mollified, method='method')    


    def propagator_ps(self, ps_in, energy, crystalline_transfer_functions=True,
                                  magnify=True, mollified=False, method='ED'):
        """ Propagator through space and the crystals of the Monochromator with the different methods. 
        
        Input: (Q, energy)
        ps_in : entry pixel size. Has to de a quantity value in the form ps_in = make_tuple(ps)
        energy: energy of the x-ray beam. Can be equal to the positioning energy self.energy
            or not. Should be a quantity object, usually expressed in keV or eV --> multiply the
            number for q.keV or q.eV.
        crystalline_transfer_functions = True / [True,True,True,True, ...] for the number of
            crystals composing the Bragg Magnifier. If crystalline_transfer_functions=False,
            all the crystal functions are removed from propagation. If one of the entry of
            crystalline_transfer_functions is False, that crystal function is removed from
            propagation. If one of the entry is 0, that 0 means False. If one of the entry is 1,
            that 1 means True.
        magnify = True. If false, the magnification is not applied.
        mollified = False or True. Apply mollification to the propagators.
        method = 'ED', 'RCT', 'Syris', 'ED Stano' -> use the Effective Distance method (ED),
            or Reciprocal Space Transformation method (RCT), using the Syris propagators or the
            Effective Distance Method from the Stano's notebooks.
        
        Output: object
        .Qout : tuple (QX_out, QY_out). Reciprocal space coordinates for X and Y direction after
                magnification. It is an nxn matrix of quantity numbers.
        .crystal_frequency_transfer : crystal frequency transfer of the crystal.
        .Q_s : array composed of arrays (QX_out, QY_out) of all the reciprocal space coordinates used for propagation.
                The number of coordinates is equal to the number of crystals + 1.
        .propagator : total propagator.
        """
        
        return self.propagator_magnifier_ps(ps_in, energy, crystalline_transfer_functions=
                                            crystalline_transfer_functions, magnify=magnify,
                                            mollified=mollified, method='method').propagator
