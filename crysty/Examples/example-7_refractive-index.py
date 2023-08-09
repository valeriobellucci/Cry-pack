#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 00:05:41 2023

@author: Valerio
"""

# Import necessary libraries

from objects import MaterialRefractiveIndex


#%%
"""
Material refractive index as a variable of photon energy.
"""

energy_eV = 12_000

material = MaterialRefractiveIndex('Materials/refractive-index_soda-lime-glass.txt')
refractive_index = material.get_refractive_index(energy_eV)
print(f"refractive_index = 1+{(refractive_index-1):.3} \n for Energy = {energy_eV/1000} keV ")

