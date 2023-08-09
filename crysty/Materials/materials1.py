"""Sample material represented by a complex refractive index."""
'''import cPickle
import logging
import os
import sys
import time
import urllib
import urllib2
from distutils.spawn import find_executable
from HTMLParser import HTMLParser
from subprocess import Popen, PIPE
from urlparse import urljoin
import quantities as q
from syris import config as cfg, physics '''
import numpy as np
from scipy import interpolate as interp



# LOG = logging.getLogger(__name__)  ''' What is this?? '''


class Material(object):

    """A material represented by its *name* and *miller_indices* calculated for *energies*."""

    def __init__(self, name, miller_indices, energies):
        """Create material with *name* and store its complex *miller_indices* (delta + ibeta)
        for all given *energies*.
        """
        self._name = name
        self._miller_indices = np.array(miller_indices)
        # To keep track which energies were used.
        self._energies = energies



    @property
    def name(self):
        """Material *name*."""
        return self._name

    @property
    def energies(self):
        """*energies* for which the complex refractive
        index was calculated.
        """
        return self._energies

    @property
    def miller_indices(self):
        """Get complex refractive indices (delta [phase], ibeta [absorption])
        for all energies used to create the material.
        """
        return self._miller_indices


silicon = Material('Si',[2, 2, 0],25)

print silicon.name 
print silicon.miller_indices
print silicon.energies



