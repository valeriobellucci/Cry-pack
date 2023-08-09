#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  7 14:43:01 2023

@author: Valerio
"""


from syris.bodies.simple import make_sphere, StaticBody
from syris.materials import make_fromfile
from syris.physics import propagate

import numpy as np
import quantities as q
import matplotlib.pyplot as plt
from utility import screen, intensity
from crystal_class import Crystal



crystal_1 = Crystal('Si', [2, 2, 0], 0. * q.deg, 29.85 * q.keV)

pixel_size_effective = 55. * q.um / 55.
detector_size_effective = pixel_size_effective * 512
qx_max = (1. / (2. * pixel_size_effective)).simplified.magnitude
qx_step = (1. / (2. * detector_size_effective)).simplified.magnitude
qx = np.arange(- qx_max / 2., + qx_max / 2., qx_step ) * (1. / q.m)
qy = qx
QX, QY = np.meshgrid(qx.simplified.magnitude, qy.simplified.magnitude) * (1. / q.m)

screen(intensity(crystal_1.QX_or_QY_transfer_fourier(QY, 29.85 * q.keV)))


n = 512
ratio = 8.
radius = n / ratio * 1. * q.um
sample = make_sphere(n, radius, pixel_size=1. * q.um, material=make_fromfile('pmma_1_31_keV.mat'))

distance_now = 1. * q.m
energies_now = np.arange(29.845, 29.856, 0.001) * q.keV
propagated = propagate([sample, crystal_1], (512, 512), energies_now, distance_now,
                       (1., 1.) * q.um).get()
                       
print(crystal_1.bragg_angle(energies_now))
print(crystal_1.bragg_angle(energies_now)[0])

screen(propagated)

transfer = crystal_1.transfer_fourier(n, 1. * q.um, 29.85 * q.keV).get()
screen(np.absolute(transfer)**2)

cry_transfer = crystal_1.crystal_transfer_fourier(n, 1. * q.um, 29.85 * q.keV)
screen(np.absolute(cry_transfer)**2)

plt.show()



