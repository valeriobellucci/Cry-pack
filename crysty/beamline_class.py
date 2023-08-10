# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 09:43:33 2018

@author: xx9905
"""

from bragg_magnifier_class import Bragg_Magnifier
from monochromator_class import Monochromator
from beamline_components import (Source, Monochromator_simple, Space, Sample, BraggMagnifier,
                                 Monochromator_now, Medipix, MedipixHexaGaAs)
from syris.util import save_image


class Beamline:
    """ Beamline class. It has to be initialized with the beamline components.
    
    __init__(*beamline_components, **components_specifications): builds the beamline with the beamline
            components and the components specifications. Load the components instances in the 
            Beamline instance.
    run(): runs the beam through the components of the beamline producing propagators and images.
            loads them in the Beamline instance.
    phantom_show_recorded_image(): Show the image recroded by the detector without plt.show(), 
            so allowing to show multiple images.
    show_recorded_image(): Show the image recorded by the detector.
    """    
    
    def __init__(self, *beamline_components, **components_specifications):
        """ All the beamline components that compose the beamline have to be added to the initialization
            of the Beamline instance with their keywords. The possible keywords are:
                Source -> keyword: 'source'
                Sample -> keyword: 'sample'
                Monochromator -> keyword: 'mono' or 'monochromator'
                Bragg Magnifier -> keyword: 'BM' or 'Bragg Magnifier'
                Detector -> 'detector'
                Space between all the components (equivalent to all the spaces below) -> 'space'
                Space between monochromator and sample -> 'space_mono_sample'
                Space between sample and Bragg Magnifier -> 'space_sample_BM'
                Space between Bragg Magnifier and detector -> 'space_BM_detector'
                Space inside the Bragg Magnifier when there is no Bragg Magnifier -> 'space_inside_BM' 
                    (also activated by 'space')
                No Space inside the Bragg Magnifier -> 'no_space_inside_BM' 
                    when there is no Bragg Magnifier, specifies that there is no space in place of the BM
                    when there is a Bragg Magnifier, puts the spaces inside the BM to zero z_01, z_12, z_23, z_34, z_45
                ### BUT 'no_space_inside_BM' CURRENTLY DOES NOT WORK SINCE THE FUNCTION THAT PROPAGATES
                    INSIDE THE BM IS IN ANOTHER SCRIPT forward_propagation AND TAKES THE z SPACE DIRECTLY FROM params
                
        Example: The beamline instance Beamline('source', 'sample', 'detector') will be composed of just 
                    a source, a sample and a detector.
            
            The details of the components composing the beamline can be specified by the keyworded 
            variables composing the components specifications. All of them have a default value.
                source_type = 'bending_magnet' or 'wiggler'
                monochromator = Monochromator_2C() or other monochromator instance
                distance_mono_sample = z_mono_sample or other quantity distance (example: distance * q.m)
                distance_sample_BM = z_sample_BM or other quantity distance (example: distance * q.m)
                sample_type = 'sphere_scaled' (sphere with radius scaled at a fraction of the Field of View), 
                                'sphere_fix' (sphere with fixed radius), 
                                'sphere_generic' (sphere or ellipse with generic size and position)
                                'disk_generic' (disk or elliptic disk with generic size and position)
                                'uniform_fourier'/'uniform fourier' (uniform ones matrix in Fourier Space)
                    sample_ratio : in the case of 'sphere_scaled', additional specification -> ratio 
                        between image field of view and sample radius. Default sample_ratio = 10.
                    sample_radius : in the case of 'sphere_fix', additional specification ->  radius 
                        of the plexiglass sphere. Default sample_radius = 100 * q.um
                distance_BM_detector = z_BM_detector or other quantity distance (example: distance * q.m)
                distance_inside_BM = z_01 + z_12 + z_23 + z_34 + z_45 from parameters script or other 
                    quantity distance (example: distance * q.m)
                detector = MedipixHexaGaAs() or Medipix() other detector instance
        """
        
        """ Selects the beamline components to add to the beamline. """
        # Gets the beamline components from __init__
        self.beamline_components = np.array([ components for components in beamline_components ])
        # List of keywords to call the components of the beamline.
        self.keywords = np.array(['source', 'mono', 'monochromator', 'mono_simple',
                                  'monochromator_simple', 'sample','BM', 'Bragg Magnifier',
                                  'detector', 'space', 'space_mono_sample', 'space_sample_BM',
                                   'space_BM_detector', 'space_inside_BM', 'no_space_inside_BM'])
        # Creates a boolean mask for filtering the beamline components that correspond to certain keywords.        
        self.logic_mask = np.sum([self.beamline_components == key for key in self.keywords],
                                 axis=0, dtype=bool)
        # Filters the beamline components that correspond to certain keywords.
        self.beamline_components = self.beamline_components[self.logic_mask]
        print('The beamline is composed of: ', self.beamline_components)

        """ Build the beamline according to the beamline components and the components specifications 
            in __init__ or default components specifications.
            Loads the Beamline instance with instances of the components of the beamline."""
        # Bragg Magnifier (keyword: BM or Bragg Magnifier)       
        if np.any(self.beamline_components == 'BM') or np.any(self.beamline_components == 'Bragg Magnifier'):
            self.there_is_BM = True
        else:
            self.there_is_BM = False        
        # Source (keyword: source)
        if np.any(self.beamline_components == 'source'):
            self.source_type = components_specifications.get('source_type', 'bending_magnet')
            self.source = Source(self.source_type, where = self.there_is_BM)
        else:
            self.source_type = None
        # Monochromator (keyword: monochromator or mono)
        if np.any(self.beamline_components == 'mono' ) or np.any(self.beamline_components == 'monochromator'):
            self.monochromator = components_specifications.get('monochromator', 'monochromator')
            self.monochromator = Monochromator_now
        elif np.any(self.beamline_components == 'mono_simple' ) or np.any(self.beamline_components == 'monochromator_simple'):
            self.monochromator = components_specifications.get('monochromator', 'monochromator_simple')
            self.monochromator = Monochromator_simple(where = self.there_is_BM)
        else:
            self.monochromator = None
        # Sample (keyword: sample)   
        if np.any(self.beamline_components == 'sample'):
            self.sample_type = components_specifications.get('sample_type', 'sphere_scaled')
            self.sample_ratio = components_specifications.get('sample_ratio', 10.)
            self.sample_radius = components_specifications.get('sample_radius', 100 * q.um)
            self.sample_x_centre = components_specifications.get('x_centre', 0)
            self.sample_y_centre = components_specifications.get('y_centre', 0)
            self.sample_xdil = components_specifications.get('xdil', 1)
            self.sample_ydil = components_specifications.get('ydil', 1)
            self.separation = components_specifications.get('separation', 0.75)
            self.sample = Sample(self.sample_type, where = self.there_is_BM, 
                                 ratio = self.sample_ratio, radius = self.sample_radius, 
                                 x_centre=self.sample_x_centre, y_centre=self.sample_y_centre, 
                                 xdil=self.sample_xdil, ydil=self.sample_ydil, 
                                 separation=self.separation)
        else:
            self.sample_type = None
        # Detector (keyword: detector)    
        if np.any(self.beamline_components == 'detector'):
            self.detector = components_specifications.get('detector', MedipixHexaGaAs())
            self.detector = self.detector
        else:
            self.detector = None
        # Spaces (keyword: space) add all the spaces
        # Space between Monochromator and Sample (keyword: space_mono_sample)    
        if np.any(self.beamline_components == 'space_mono_sample') or np.any(self.beamline_components == 'space'):
            self.distance_mono_sample = components_specifications.get('distance_mono_sample', z_mono_sample)
            self.space_mono_sample = Space(self.distance_mono_sample)
        else:
            self.distance_mono_sample = None
        # Space between Sample and Bragg Magnifier (keyword: space_sample_BM)    
        if np.any(self.beamline_components == 'space_sample_BM') or np.any(self.beamline_components == 'space'):
            self.distance_sample_BM = components_specifications.get('distance_sample_BM', z_sample_BM)
            self.space_sample_BM = Space(self.distance_sample_BM)
        else:
            self.distance_sample_BM = None
        # Space inside the Bragg Magnifier (keyword: no_space_inside_BM, space_inside_BM)
        if np.any(self.beamline_components == 'no_space_inside_BM'):
            self.distance_inside_BM = None
        elif (self.there_is_BM == True) & np.any(self.beamline_components == 'no_space_inside_BM'):
            distances_BM = np.zeros(number_of_crystals_BM) * q.m
            self.distance_inside_BM = np.sum(distances_BM.simplified.magnitude) * q.m
        elif ( (self.there_is_BM == False) & np.any(self.beamline_components == 'space') | 
                (self.there_is_BM == False) & np.any(self.beamline_components == 'space_inside_BM') ): 
            self.default_distance_inside_BM = np.sum(distances_BM.simplified.magnitude) * q.m
            self.distance_inside_BM = components_specifications.get('distance_inside_BM', self.default_distance_inside_BM )
            self.space_inside_BM = Space(self.distance_inside_BM)
        elif self.there_is_BM == True:
            self.distance_inside_BM = np.sum(distances_BM.simplified.magnitude) * q.m
        else:
            self.distance_inside_BM = None
        # Space between Bragg Magnifier and Detector (keyword: space_BM_detector)    
        if np.any(self.beamline_components == 'space_BM_detector') or np.any(self.beamline_components == 'space'):
            self.distance_BM_detector = components_specifications.get('distance_BM_detector', z_BM_detector)
            self.space_BM_detector = Space(self.distance_BM_detector) 
        else:
            self.distance_BM_detector = None
        
        if ( (self.sample_type != None) & (self.sample_type != 'uniform_fourier') 
            & (self.sample_type != 'uniform fourier') ): 
            # Calculates the expected Fresnel number for the object in the beamline    
            self.fresnel_distance = ( none_to_zero_m(self.distance_sample_BM) 
                                        + none_to_zero_m(self.distance_inside_BM) 
                                        + none_to_zero_m(self.distance_BM_detector) )
            if ( (self.sample_type == 'sphere_fix') | (self.sample_type == 'sphere_scaled') 
                | (self.sample_type == '4_spheres_scaled') | (self.sample_type == '4_disks_scaled') 
                | (self.sample_type == 'sphere_generic') | (self.sample_type == 'disk_generic') ):
                self.sample_size = self.sample.radius
            print('The Fresnel number is: ', self.fresnel_number)
            if self.fresnel_number < 1.:
                print('    Far-field regime')
            elif self.fresnel_number > 1.:
                print('    Near-field regime')
    
    
    def run(self, save = False, save_path = '', save_format = '.tiff', computation = 'syris', 
            BM_propagation_method = 'ED', exp_time = 1 * q.s, 
            stop_after = None, stop_before = None, test = False, 
            sample_interaction = 'absorption and phase object', 
            crystalline_transfer_functions=True, magnify=True, mollified=False):
        """ Runs the beam through the beamline producing the propagators and the images. 
            Loads the Beamline instance with them. 
            
        Possible keyworded variables in running the beam are:
            save = False or True (anything different from True is considered False): saves the real 
                and imaginary parts of transfers,propagators and intermediate images to files, in the 
                path specified by save_path and with format specified by save_format.
            save_path : specifies the path where to save the files if save = True. Default is '' so it
                saves in the running folder.
            save_format : specifies in what format to save the arrays if save = True. Default is '.tiff'
            computatio = 'syris' or 'analytic' -> the free space propagators are produced by syris or
                                                    with their analytic formula.
            BM_propagation_method = 'ED' or 'RCT' or 'Syris'
                                    'ED' --> effective distance method. 
                                    'RCT' --> reciprocal coordinate transformation method.
                                    'Syris' --> effective distance method using the propagators from Syris
            stop_before/stop_after : stop before or after a certain element
                stop_before: 'mono', 'monochromator', 'sample', 'BM', 'Bragg Magnifier', 'detector'
                stop_after: 'source', 'mono', 'monochromator', 'sample', 'BM', 'Bragg Magnifier'
            test = False or True --> if True activates the test mode: it stores the BM_fs transfer function.
            sample_interaction = 'absorption and phase object' -> default, the sample object is an
                                        absorption and phase mixed object
                                     'pure phase object' -> the sample is a pure phase object
                                     'pure absorption object' -> the sample is a pure absorption object
            there_are_crystals=[True,True,True,True]. If there_are_crystals=False, all the crystals are 
                removed from propagation. If one of the entrances of there_are_crystals=[1st,2nd,3rd,4th] is
                False, that nth crystal is removed from propagation. If one of the entrances is 0, that 0 
                means False. If one of the entrances is 1, that one means True.
        """
        print('Runs beamline composed of: ', self.beamline_components)
        
        
        def save_specifications():
            """ Saves the beamline and running specifications in a txt file when keyworded variable 
                save in the run methos is save = True. 
            """
            specs = np.array(['BEAMLINE SPECIFICATIONS\n',
                        'save path: ' + save_path,
                        'save format of the arrays: ' + save_format,
                        
                        '\nBeamline components: ' + str(self.beamline_components),

                        '\nComponents specifications',
                        'Source type: ' + str(self.source_type),
                        'Monochromator type: ' + str(self.monochromator),
                        'Sample type: ' + str(self.sample_type), 
                        '    Sample ratio: ' + str(self.sample_ratio) + ''' in the case of sphere_scaled ''',
                        '    Sample radius: ' + str(self.sample_radius) + ''' in the case of sphere_fix ''',
                        'Detector type: ' + str(self.detector),
                        'Distance between monochromator and sample = ' + str(self.distance_mono_sample),
                        'Distance between sample and Bragg Magnifier = ' + str(self.distance_sample_BM),
                        'Distance between Bragg Magnifier and detector = ' + str(self.distance_BM_detector),
                        
                        '\nRun specifications',
                        'computation of the free space propagators: ' + computation,
                        'Bragg Magnifier propagation method: ' + BM_propagation_method, 
                        'Detector exposure time = ' + str(exp_time),
                        'Stop before ' + str(stop_before) ,
                        'Stop after ' + str(stop_after)
                        ])
            np.savetxt(save_path + 'Beamline specifications.txt', specs, fmt="%s")
        
        if save == True: save_specifications()
        
        def save_run(array, array_name):
            """ If the keyworded variable save in the run methos is save = True, the save_run function
                Saves the complex array with name array_name to the path save_path (keyworded in run) 
                and to the format save_format (keyworded in run)           
            """
            if (save == True) & (type(array) == np.ndarray):
                save_image(save_path + array_name + '_real' + save_format, array.real)
                save_image(save_path + array_name + '_imag' + save_format, array.imag)


        # RUNS THE BEAMLINE
        print('CALCULATING: ')        
        # Source        
        if self.source_type == None:
            self.s_0 = np.ones((n,n))
            self.s_0_fs = fft2c(self.s_0)
        else:
            print('source transfer function')
            self.s_0 = self.source.transfer(energy)
            print('source transfer function in fourier space')
            self.s_0_fs = fft2c(self.s_0)
        save_run(self.s_0, 'source_transfer')
        save_run(self.s_0_fs, 'source_transfer_fourier_space')
        if stop_after == 'source': return        
        if (stop_before == 'mono') | (stop_before == 'monochromator'): return                   
        # Monochromator
        if self.monochromator == None:
            self.mono_0_fs = 1.
            self.mono_0 = 1.
        else:
            if np.any(self.beamline_components == 'mono' ) or np.any(self.beamline_components == 'monochromator'):
                if self.there_is_BM == False:
                    print('monochromator transfer function in fourier space ') 
                    self.mono_0_fs = self.monochromator.transfer(Q_exit, energy)
                elif self.there_is_BM == True:
                    print('monochromator transfer function in fourier space ')  
                    self.mono_0_fs = self.monochromator.transfer(Q_magn, energy)
                print('monochromator transfer function in real space')
                self.mono_0 = ifft2c(self.mono_0_fs)
            elif np.any(self.beamline_components == 'mono_simple' ) or np.any(self.beamline_components == 'monochromator_simple'):
                print('monochromator transfer function in fourier space ')  
                self.mono_0_fs = self.monochromator.transfer(energy)
                print('monochromator transfer function in real space')
                self.mono_0 = ifft2c(self.mono_0_fs)
            
        save_run(self.mono_0_fs, 'monochromator_transfer_fourier_space')
        save_run(self.mono_0, 'monochromator_transfer')
        if (stop_after == 'mono') | (stop_after == 'monochromator'): return           
        # Space between Monochromator and Sample
        if self.distance_mono_sample == None:
            self.propagator_mono_sample_fs = 1.
        elif self.there_is_BM == False:
            print('propagator of the space between BM and detector (without Magnifier), in Fourier space')
            self.propagator_mono_sample_fs = self.space_mono_sample.propagator_after_BM(energy,
                                                                                        computation=computation)
        else:
            print('propagator of the space between the monochromator and the sample')
            self.propagator_mono_sample_fs = self.space_mono_sample.propagator_before_BM(energy,
                                                                                         computation=computation)        
        save_run(self.propagator_mono_sample_fs, 'space_propagator_between_mono_and_sample_fs')
        if (stop_before == 'sample'): return                      
        # Sample
        if self.sample_type == None:
            self.u_0 = np.ones((n,n), dtype=np.complex64)
            self.u_0_fs = fft2c(self.u_0)
        elif (self.sample_type == 'uniform_fourier') | (self.sample_type == 'uniform fourier'):
            self.u_0_fs = np.ones((n,n), dtype=np.complex64)
            self.u_0 = ifft2c(self.u_0_fs)
        else:
            print('object transfer function')
            self.u_0 = self.sample.transfer(energy)
            self.sample_interaction = sample_interaction
            if self.sample_interaction == 'pure phase object': # Turns the sample in a pure phase object
                self.u_0 = self.u_0 / np.absolute(self.u_0)
            elif self.sample_interaction == 'pure absorption object': # Turns the sample in a pure absorption object
                self.u_0 = self.u_0 / phase(self.u_0)
            print('object transfer function in fourier space')
            self.u_0_fs = ifft2c(self.u_0)
        save_run(self.u_0, 'sample_transfer')
        save_run(self.u_0_fs, 'sample_transfer_fourier_space')
        if stop_after == 'sample': return              
        # Space between Sample and Bragg Magnifier
        if self.distance_sample_BM == None:
            self.propagator_sample_BM_fs = 1.
        elif self.there_is_BM == False:
            print('propagator of the space between BM and detector (without Magnifier), in Fourier space')
            self.propagator_sample_BM_fs = self.space_sample_BM.propagator_after_BM(energy,
                                                                                    computation=computation)
        else:
            print('propagator of the space between sample and BM')
            self.propagator_sample_BM_fs = self.space_sample_BM.propagator_before_BM(energy,
                                                                                     computation=computation) 
        save_run(self.propagator_sample_BM_fs, 'space_propagator_between_sample_and_BM_fs')
        # Wavefield at the entrance of the BM, in Fourier space
        print('wavefield at the entrance of the BM, in Fourier space')
        self.u_before_BM_fs = ( fft2c(self.u_0 * 
                                        ifft2c(self.s_0_fs * self.mono_0_fs * self.propagator_mono_sample_fs)
                                        ) 
                                * self.propagator_sample_BM_fs )       
        save_run(self.u_before_BM_fs, 'Wavefield at the entrance of the BM, in Fourier space')
        if (stop_before == 'BM') | (stop_before == 'Bragg Magnifier'): return              
        # Space inside the Bragg Magnifier
        if (self.there_is_BM == False) & (self.distance_inside_BM != None):
            print('propagator of the space in place of the Bragg Magnifier (without Magnifier), in Fourier space')
            self.propagator_inside_BM_fs = self.space_inside_BM.propagator_after_BM(energy, computation=computation) 
        elif self.distance_inside_BM == None:
            self.propagator_inside_BM_fs = 1.
        else:
            self.propagator_inside_BM_fs = 1.
        save_run(self.propagator_inside_BM_fs, 'space_propagator_in_place_of_the_BM_fs')        
        # Bragg Magnifier
        if self.there_is_BM == False:
            self.u_after_BM_fs = self.u_before_BM_fs * self.propagator_inside_BM_fs
        elif self.there_is_BM == True:
            print('wavefield after the BM')
            self.u_after_BM_fs = BraggMagnifier().wavefield(self.u_before_BM_fs, energy,
                                                            method = BM_propagation_method,
                                                            crystalline_transfer_functions=
                                                            crystalline_transfer_functions,
                                                            magnify=magnify, mollified=mollified) 
            if test == True:
                self.BM_fs = BraggMagnifier().wavefield(np.ones((n,n)), energy,
                                                        method = BM_propagation_method,
                                                            crystalline_transfer_functions=
                                                            crystalline_transfer_functions,
                                                            magnify=magnify, mollified=mollified)     
        save_run(self.u_after_BM_fs, 'wavefield at the exit of the BM, in Fourier space')
        if (stop_after == 'BM') | (stop_after == 'Bragg Magnifier'): 
            self.u_after_BM = ifft2c(self.u_after_BM_fs)
            return              
        # Space between Bragg Magnifier and Detector
        if self.distance_BM_detector == None:
            self.propagator_BM_d_fs = 1.
        else:        
            print('propagator of the space between BM and detector, in Fourier space')
            self.propagator_BM_d_fs = self.space_BM_detector.propagator_after_BM(energy,
                                                                                 computation=computation)
        save_run(self.propagator_BM_d_fs, 'space_propagator_between_BM_and_detector_fs')
        # Wavefield just before the detector, in Fourier space
        print('wavefield just before the detector, in Fourier space')
        self.u_d_fs = self.u_after_BM_fs * self.propagator_BM_d_fs   
        save_run(self.u_d_fs, 'Wavefield just before the detector, in Fourier space')        
        # Wavefield just before the detector, in real space
        print('wavefield just before the detector, in real space')
        self.u_d = ifft2c(self.u_d_fs)
        save_run(self.u_d, 'Wavefield just before the detector, in real space')             
        # Intensity image just before the detector
        print('intensity image just before the detector')
        self.i_d = intensity(self.u_d)
        save_run(self.i_d, 'Intensity image just before the detector')          
        if stop_before == 'detector': return      
        #screen(self.i_d)
        #plt.show()
        # Detector
        if self.detector == None:
            self.recorded_image = self.i_d
        else:
            print('image recorded by the detector')
            self.recorded_image = self.detector.get_image(self.i_d, exp_time=exp_time)
        save_run(self.recorded_image, 'Image recorded by the detector')  
        print('\n')
    
    def phantom_show_recorded_image(self):
        ''' Show the image recroded by the detector without plt.show(), so allowing to show multiple images. '''
        plt.figure()
        initial = plt.imshow(self.recorded_image, interpolation = 'none')
        plt.colorbar(initial)
    
    def show_recorded_image(self):
        ''' Show the image recorded by the detector '''
        self.phantom_show_recorded_image()
        plt.show()



"""
''' Position invariance??? '''
beamline_1 = Beamline( 'sample', 'BM', 'detector', 'space_sample_BM',
                      sample_type = 'rectangle', x_pixels=n/2, y_pixels=n/2, x_centre=0.)
beamline_1.run(exp_time = 1. * q.s, BM_propagation_method = 'RCT', there_are_crystals=[0,0,0,0])

print 'Fresnel number : ', ((1. * q.um * 1024 / 4.)**2 / (lam * 2. * q.m)).simplified

screen(beamline_1.u_0)
screen(beamline_1.u_0[:,::-1])
screen(beamline_1.u_0 - beamline_1.u_0[:,::-1])
image1 = beamline_1.i_d
screen(image1)

image1_horizontal_mirror = beamline_1.i_d[:,::-1]
screen(image1_horizontal_mirror)
image2 = beamline_1.i_d
screen(image2)
horizontal_subtraction = image2 - image1_horizontal_mirror
screen(horizontal_subtraction)

image3 = beamline_1.i_d
image3_vertical_mirror = beamline_1.i_d[::-1,:]
screen(image3_vertical_mirror)
image4 = beamline_1.i_d
vertical_subtraction = image4 - image3_vertical_mirror
screen(vertical_subtraction)

plt.show()
"""
"""
beamline_1 = Beamline( 'sample', 'BM', 'detector', 'space_sample_BM',
                      sample_type = 'sphere_fix', sample_radius = 25. * q.um)
beamline_1.run(exp_time = 1. * q.s, BM_propagation_method = 'RCT', there_are_crystals=[0,0,0,0])
print 'beamline 1'
beamline_2 = Beamline( 'sample', 'BM', 'detector', 'space_sample_BM',
                      sample_type = 'sphere_generic', sample_radius = 25. * q.um)
beamline_2.run(exp_time = 1. * q.s, BM_propagation_method = 'RCT', there_are_crystals=[0,0,0,0])
print 'beamline 2'
beamline_3 = Beamline( 'sample', 'BM', 'detector', 'space_sample_BM',
                      sample_type = 'sphere_generic', sample_radius = 25. * q.um)
beamline_3.run(exp_time = 1. * q.s, BM_propagation_method = 'RCT', there_are_crystals=[0,0,0,0])
print 'beamline 3'
beamline_4 = Beamline( 'sample', 'BM', 'detector', 'space_sample_BM',
                      sample_type = 'sphere_generic', sample_radius = 25. * q.um)
beamline_4.run(exp_time = 1. * q.s, BM_propagation_method = 'RCT', there_are_crystals=[0,0,0,0])
print 'beamline 4'

image1 = beamline_1.i_d
screen(image1)
print image1[512/2 - 1 - 2 : 512/2 + 3, 512/2 - 1 - 2 : 512/2 + 3]
image1_horizontal_mirror = beamline_1.i_d[:,::-1]
screen(image1_horizontal_mirror)
screen(image1_horizontal_mirror - image1)
image1_vertical_mirror = beamline_1.i_d[::-1,:]
screen(image1_vertical_mirror)
screen(image1_vertical_mirror - image1)"""
"""
image2 = beamline_2.i_d
screen(image2)
subtraction12 = image1_horizontal_mirror - image2
screen(subtraction12)

image3 = beamline_3.i_d
image3_horizontal_mirror = beamline_3.i_d[::-1,:]
image4 = beamline_4.i_d
subtraction34 = image3_horizontal_mirror - image4
screen(image3)
screen(image3_horizontal_mirror)
screen(image4)
screen(subtraction34)
              
plt.show()
"""

"""Test the difference between the ED and RCT method while varying the sample to BM propagation space.
Far-field: distance_sample_BM = 100. * q.m, sample_radius = 5. * q.um, Fresnel number = 0.006
Mid-field: distance_sample_BM = 5. * q.m, sample_radius = 15. * q.um, Fresnel number = 1
Near-field: distance_sample_BM = 0. * q.m, sample_radius = 50. * q.um, Fresnel number = 74"""

distance_sample_BM = 0. * q.m
sample_radius = 100. * q.um
# 100-2 um for x55
beamline_1 = Beamline('sample', 'BM', 'space_sample_BM', distance_sample_BM = distance_sample_BM,
                      sample_type = '4_disks_scaled', sample_radius = sample_radius, separation = 0.5)
beamline_1.run(BM_propagation_method = 'ED', sample_interaction = 'pure phase object'
                , save = True, save_path = 'transfer_functions_bm1/')
beamline_2 = Beamline('sample', 'BM', 'space_sample_BM', distance_sample_BM = distance_sample_BM,
                      sample_type = '4_disks_scaled', sample_radius = sample_radius, separation = 0.5)
beamline_2.run(BM_propagation_method = 'RCT', sample_interaction = 'pure phase object'
                , save = True, save_path = 'transfer_functions_bm2/')
print 'distance_inside_BM : ', beamline_1.distance_inside_BM
print 'Fresnel number : ', fresnel_number(sample_radius, distance_sample_BM + z_01, lam)

difference_max_ratio = ( (np.amax(beamline_2.i_d - beamline_1.i_d) - np.amin(beamline_2.i_d - beamline_1.i_d)) 
                    / (np.amax(beamline_1.i_d) - np.amin(beamline_1.i_d)))
print 'difference_max_ratio : ', difference_max_ratio * 100, ' %'
print 'wavelength : ', lam
print 'magnification : ', M

#beamline_1.run(save = True, save_path = 'transfer_functions/', exp_time = 1 * q.s, stop_after = 'sourc')


screen(beamline_1.u_0.real, plot_title = 'Object real part', fontsize = 22)
screen(beamline_1.i_d, plot_title = 'Intensity ED', fontsize = 22)
screen(beamline_2.i_d, plot_title = 'Intensity RCT', fontsize = 22)
screen(beamline_2.i_d - beamline_1.i_d, plot_title = 'Intensity RCT-ED', fontsize = 22)
screen(QX.simplified.magnitude)
screen(QY.simplified.magnitude)

difference_mean_ratio = np.sum(np.absolute(beamline_2.i_d - beamline_1.i_d))  / np.sum(beamline_1.i_d)
print 'difference_mean_ratio : ', difference_mean_ratio * 100, ' %'

plt.figure()
plt.plot(beamline_1.i_d[512,:]) # plot over x
plt.title('X profile of ED')
plt.xlabel('pixels ')
plt.ylabel('Intensity ')

plt.figure()
plt.plot(beamline_1.i_d[:,512]) # plot over y
plt.title('Y profile of ED')
plt.xlabel('pixels ')
plt.ylabel('Intensity ')

plt.figure()
plt.plot(beamline_2.i_d[512,:]) # plot over x
plt.title('X profile of RCT')
plt.xlabel('pixels ')
plt.ylabel('Intensity ')

plt.figure()
plt.plot(beamline_2.i_d[:,512]) # plot over y
plt.title('Y profile of RCT')
plt.xlabel('pixels ')
plt.ylabel('Intensity ')

plt.figure()
plt.plot((beamline_2.i_d - beamline_1.i_d)[512,:]) # plot over x
plt.title('X profile of RCT - ED')
plt.xlabel('pixels ')
plt.ylabel('Intensity ')

plt.figure()
plt.plot((beamline_2.i_d - beamline_1.i_d)[:,512]) # plot over y
plt.title('Y profile of RCT - ED')
plt.xlabel('pixels ')
plt.ylabel('Intensity ')


plt.show()



"""
beamline_1 = Beamline('mono','BM')
beamline_1.run(test = True)
mono = beamline_1.mono_0_fs
BM = beamline_1.BM_fs


screen(intensity(mono))
screen(intensity(BM))
plt.show()
"""
"""    
#beamline_1 = Beamline(distance_BM_detector = 0. * q.m, distance_mono_sample = 0. * q.m, distance_sample_BM = 0. * q.m)
beamline_1 = Beamline('sourc', 'mon', 'sample', 'BM', 'detecto', 
                      'spac', 'space_inside_B', 'space_BM_detecto', 'space_sample_B',
                      source_type = 'wiggler', sample_type = 'sphere_scaled', sample_ratio = 10.,
                      distance_BM_detector = 12000000. * q.mm, distance_mono_sample = 2.0 * q.m, 
                      distance_sample_BM = 12000000. * q.m)
beamline_2 = beamline_1

beamline_1.run( BM_propagation_method = 'RCT', exp_time = 1 * q.s, stop_after = 'BM')
beamline_2.run( BM_propagation_method = 'ED', exp_time = 1 * q.s, stop_after = 'BM')

print '\n phase(beamline_1.u_0) \n', phase(beamline_1.u_0)
print '\n intensity(beamline_1.u_0) \n', intensity(beamline_1.u_0)
print ('\n phase((beamline_1.u_after_BM-beamline_1.u_0) - (beamline_2.u_after_BM-beamline_2.u_0)) \n', 
       phase((beamline_1.u_after_BM-beamline_1.u_0) - (beamline_2.u_after_BM-beamline_2.u_0))
       )
print ('\n intensity((beamline_1.u_after_BM-beamline_1.u_0) - (beamline_2.u_after_BM-beamline_2.u_0)) \n', 
       intensity((beamline_1.u_after_BM-beamline_1.u_0) - (beamline_2.u_after_BM-beamline_2.u_0))
       )

screen(phase(beamline_1.u_0))
screen(intensity(beamline_1.u_0))
screen(phase(beamline_1.u_after_BM))
screen(intensity(beamline_1.u_after_BM))
screen(phase(beamline_2.u_after_BM))
screen(intensity(beamline_2.u_after_BM))
screen( phase((beamline_1.u_after_BM-beamline_1.u_0) - (beamline_2.u_after_BM-beamline_2.u_0)) )
screen( intensity((beamline_1.u_after_BM-beamline_1.u_0) - (beamline_2.u_after_BM-beamline_2.u_0)) )
plt.show()

screen(phase(beamline_1.u_0))
screen(intensity(beamline_1.u_0))
screen(phase(beamline_1.u_after_BM))
screen(intensity(beamline_1.u_after_BM))
screen(phase(beamline_2.u_after_BM))
screen(intensity(beamline_2.u_after_BM))
screen( (phase(beamline_1.u_after_BM)-phase(beamline_1.u_0)) 
        - (phase(beamline_2.u_after_BM)-phase(beamline_2.u_0)) )
screen( (intensity(beamline_1.u_after_BM)-intensity(beamline_1.u_0)) 
        - (intensity(beamline_2.u_after_BM)-intensity(beamline_2.u_0)) )
plt.show()
"""
"""
print beamline_1.source_type
#beamline_1.run(save = True, save_path = 'transfer_functions/', exp_time = 1 * q.s, stop_after = 'sourc')
beamline_1.run( BM_propagation_method = 'RCT', exp_time = 1 * q.s, stop_after = 'BM')


print '\n phase(beamline_1.u_0) \n', phase(beamline_1.u_0)
print '\n intensity(beamline_1.u_0) \n', intensity(beamline_1.u_0)
print '\n phase(beamline_1.u_after_BM) \n', phase(beamline_1.u_after_BM)
print '\n intensity(beamline_1.u_after_BM) \n', intensity(beamline_1.u_after_BM)
print '\n phase(beamline_1.u_after_BM)-phase(beamline_1.u_0) \n', phase(beamline_1.u_after_BM)-phase(beamline_1.u_0)
print '\n intensity(beamline_1.u_after_BM)-intensity(beamline_1.u_0) \n', intensity(beamline_1.u_after_BM)-intensity(beamline_1.u_0)

screen(phase(beamline_1.u_0))
screen(intensity(beamline_1.u_0))
screen(phase(beamline_1.u_after_BM))
screen(intensity(beamline_1.u_after_BM))
screen(phase(beamline_1.u_after_BM)-phase(beamline_1.u_0))
screen(intensity(beamline_1.u_after_BM)-intensity(beamline_1.u_0))
plt.show()"""
"""
beamline_1 = Beamline('source','mono', 'BM','sample','detector', source_type = 'wiggler')
beamline_1.run()

screen(intensity(beamline_1.u_0))
#screen(phase(beamline_1.u_d))
screen(beamline_1.i_d)
beamline_1.phantom_show_recorded_image()
plt.show()"""
"""
#print beamline_1.s_0
print ''
print beamline_1.propagator_BM_d_fs
print beamline_1.distance_BM_detector
print beamline_1.beamline_components
screen(beamline_1.i_d)
screen(np.log10(intensity(beamline_1.u_0_fs)))
#screen(intensity(beamline_1.mono_0_fs))
screen(intensity(beamline_1.u_0))

beamline_1.show_recorded_image()
plt.show()
"""
