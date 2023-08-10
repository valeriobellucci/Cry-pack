#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 00:05:41 2023

@author: Valerio


Initialize the MaterialRefractiveIndex with data from a file.
Get the table of refractive index (Delta, Beta) vs Energy (eV) 
in a .txt file from https://henke.lbl.gov/optical_constants/getdb2.html
:param file_path: Path to the .txt file containing energy, delta, and beta values.

"""

# Import necessary libraries

from crysty.objects import MaterialRefractiveIndex


#%%
"""
Material refractive index as a variable of photon energy.
"""

energy_eV = 12_000

material = MaterialRefractiveIndex('refractive-index_soda-lime-glass.txt')
refractive_index = material.get_refractive_index(energy_eV)
print(f"refractive_index = 1+{(refractive_index-1):.3} \n for Energy = {energy_eV/1000} keV ")

