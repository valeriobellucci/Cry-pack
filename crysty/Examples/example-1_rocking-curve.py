#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Plots intensity and phase vs rocking angle

"""

# Load these packages
import numpy as np
import matplotlib.pyplot as plt
import cry_functions as crf


mnow = 'Si'
mill_ind = [2,2,0]
energy_keV_now = 31 
asy_best = crf.asymmetry_best(energy_keV_now, mnow, mill_ind)
print (asy_best)
asy_now = asy_best
asy_now = - np.pi/2 - (5.66 * (np.pi)/180)
asy_now = - np.pi/2 - (3.35 * (np.pi)/180)
asy_now = - np.pi/2 
#asy_now = - np.pi / 2 - 16.45
#asy_now = - (np.pi / 2) - 5.92 * np.pi / 180
step = 0.025
maxrange = 2 * crf.arcsec(crf.width_angle_o(energy_keV_now, asy_now, mnow, mill_ind))
minrange = - maxrange - step/2
maxrange_freq = crf.arcsec_to_frequency(maxrange, energy_keV_now)
minrange_freq = crf.arcsec_to_frequency(minrange, energy_keV_now)
step_freq = maxrange_freq / 100
maxrange_keV = 0.002
step_keV = maxrange_keV / 200
minrange_keV = - maxrange_keV - step_keV/2

mis_angle_now = 0
range_kev_now = 0
rock_angl_now = 0.1
cent_now = True
paral_now = True


pre_mag = 1 # Was the beam already magnified with a magnification pre_mag?
vector_rocking_angle = np.arange(minrange,maxrange,step)
vector_cry = crf.wave_crystal_function_exit(vector_rocking_angle, asy_now,
                                   energy_keV_now, mnow, mill_ind, pre_mag)[0]
vector_cry_poly = crf.wave_crystal_poly_exit(vector_rocking_angle, mis_angle_now, asy_now, 
                                    energy_keV_now, range_kev_now, mnow, mill_ind, cent_now, 
                                    paral_now, pre_mag)[0]


fig, ax1 = plt.subplots()

# Plotting the intensity with the left y-axis
ax1.plot(vector_rocking_angle, crf.intensity(vector_cry_poly), 'b-')
ax1.set_xlabel("angle (arcsec)")
ax1.set_ylabel('Intensity', color='b')
ax1.tick_params('y', colors='b')

# Create a second y-axis that shares the same x-axis
ax2 = ax1.twinx()
# Plotting the phase with the right y-axis
ax2.plot(vector_rocking_angle, np.angle(vector_cry_poly), 'r-')
ax2.set_ylabel('Phase change (rad)', color='r')
ax2.tick_params('y', colors='r')

mnow = 'Si'
mill_ind = [2,2,0]
energy_keV_now = 31 #10.76

plt.title(f"intensity and phase vs rocking angle \n \
          {mnow} {tuple(mill_ind)}, {energy_keV_now} keV")
plt.show()
