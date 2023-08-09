#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  9 17:49:16 2023

@author: Valerio
"""

from setuptools import setup, find_packages

setup(
    name="crysty",
    version="0.1",
    packages=find_packages(),
    # List your package's dependencies here
    install_requires=["numpy", "quantities", "scipy", "matplotlib",
                      "scikit-image", "pandas", "pyopencl", "math"

    ],
)
