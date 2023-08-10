# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 15:45:45 2017

@author: Valerio
"""

""" Importing packages """
import os

# NumPy - the fundamental package for scientific computing with Python.
import numpy as np
# Collection of constants
import scipy.constants as sci
# Optimizing functions
import scipy.optimize as opt
# Deals with various types of I/O (input/output). 
import io
# find and access resources inside packages.
from pkg_resources import resource_string # gets a resource from a package
import pandas as pd


""" Program for calculating crystalline functions.

This is the most updated version of the program!

This program needs to import some .txt tables containing properties of materials

The functions can take as entrance a numpy array

The asymmetry angle is expressed in radians,
symmetric Bragg reflection is obtained for asymmetry = - Pi / 2
asymmetry angle for grazing incidence is negative

"""


""" Basic physical constants """
electron_radius = sci.physical_constants['classical electron radius'][0]  # 2.81794*10**-15 m
h_planck = sci.physical_constants['Planck constant in eV s'][0]  # 4.1356692*10**-15 eV s
polariz = 1.0  # Polarization factor now
# sci.c = c = speed of light in vacuum = 2.99792458*10**8 m/s
# electronic density
# sci.m_e = me = electron mass = 9.109*10**-31
# sci.epsilon_0 = Epo = the electric constant (vacuum permittivity) = 8.854*10**-12
# sci.e = elementary charge = 1.602*10**-19


""" Some useful functions """
def arcsec(radians):
    """ Transforms radians in arcsec """
    return radians * (180.0 * 3600.0) / np.pi


def rad(arcseconds):
    """ Transforms arcsec in radians """
    return arcseconds * np.pi / (180.0 * 3600.0)


def wavelength(energy_kev):
    """ Retuns the wavelength in metres  """
    return (h_planck * sci.speed_of_light) / (energy_kev * 1000.0)


def frequency(energy_kev):
    """ Retuns the frequency  """
    return sci.c / wavelength(energy_kev)


def phase(wave_function):
    """ Returns the phase """
    return np.angle(wave_function)


def modulo(wave_function):
    """ Returns the modulus """
    return np.absolute(wave_function)


def intensity(wave_function):
    """ Returns the intensity of a wave """
    return (np.absolute(wave_function)) ** 2


def frange(start, stop, step):
    """ range function for float numbers """
    i = start
    while i < stop:
        yield i
        i += step


def printmatrix(matrix):
    """ Prints a matrix in a matrix-like form """
    for i in matrix:
        print(i)


def crysty_txt_to_list(crysty_subfolder, file_name, pandas=False):
    """ Import from crysty.materials a txt content to a NumPy array """
    data = resource_string(f'crysty.{crysty_subfolder}', file_name)
    
    # The 'StringIO' class allows us to treat strings as file-like objects.
    buffer = io.StringIO(data.decode('ISO-8859-1'))
    load = np.loadtxt(buffer, dtype=object, encoding='ISO-8859-1')
    
    if pandas is True:
        load = pd.read_csv(buffer, sep='\t')  # Assuming the data is tab-separated

    # Check if the content is just numbers by trying to convert to float64
    try:
        load = load.astype(np.float64)
    except ValueError:
        pass  # If not all values are convertible to float64, don't convert

    return load


def materials_txt_to_list(file_name):
    """ Import from crysty.materials a txt content to a NumPy array """

    return crysty_txt_to_list("materials", file_name)


""" Import txt files containing crystalline constants """
current_directory = os.getcwd()


""" Some materials constants """
costantim = materials_txt_to_list('materials_constants.txt')


""" Define a dictionary of materials """
MATERIALS = {
    # unit_cell describes the position of the atoms in the unit cell 
    'Si': {
        'index': 0,
        'unit_cell_vectors': materials_txt_to_list("unit cell Si.txt")
    },
    'Ge': {
        'index': 1,
        'unit_cell_vectors': materials_txt_to_list("unit cell Ge.txt")
    }
}

""" Form factor of materials - real and imaginary part """
form_factor_real_Si220 = 8.66557
form_factor_real_Si111 = 10.3811
form_factor_real_Ge220 = 23.7972
form_factor_real_Ge111 = 27.3664
form_factor_imag_Si = materials_txt_to_list("f2_Si_50keV.txt")
form_factor_imag_Ge = materials_txt_to_list("f2_Ge_50keV.txt")


""" Polarizability of materials in the forward direction - real and imaginary part """
chio_Si_real = materials_txt_to_list("chio_Si_real.txt")
chio_Si_imag = materials_txt_to_list("chio_Si_imag.txt")
chio_Ge_real = materials_txt_to_list("chio_Ge_real.txt")
chio_Ge_imag = materials_txt_to_list("chio_Ge_imag.txt")


""" Crystalline constants

    Insert the name of the material 'Name' in string_material
    and the vector of the Miller index [xyz] in list_miller
"""

def get_material_property(string_material, property_index):
    """Utility function to fetch a property from costantim based on material."""
    return costantim[MATERIALS[string_material]['index']][property_index]

def density(string_material):
    """Density of a material in g/cm**3."""
    return get_material_property(string_material, 4)

def electron_density(string_material):
    """Electron density of a material."""
    dens = get_material_property(string_material, 4)
    atomic_number = get_material_property(string_material, 2)
    mass_number = get_material_property(string_material, 3)
    return dens * sci.N_A * (atomic_number / mass_number) * 10 ** 6


def get_unit_cell_vector(string_material, vector_index):
    """Utility function to fetch a unit cell vector."""
    unit_cell_vector = MATERIALS[string_material]['unit_cell_vectors'][vector_index][1:]
    return np.asarray(unit_cell_vector, dtype=float) * 10 ** -10

def a1(string_material):
    """1st vector of the unit cell (metres)."""
    return get_unit_cell_vector(string_material, 1)

def a2(string_material):
    """2nd vector of the unit cell (metres)."""
    return get_unit_cell_vector(string_material, 2)

def a3(string_material):
    """3rd vector of the unit cell (metres)."""
    return get_unit_cell_vector(string_material, 3)


def get_reciprocal_vector(string_material, a_vector1, a_vector2):
    """Utility function to fetch a reciprocal vector."""
    cross_product = np.cross(a_vector1, a_vector2)
    return 2 * np.pi * cross_product / volume(string_material)

def volume(string_material):
    """Volume of the unit cell (metres)"""
    a1_vector = a1(string_material)
    a2_vector = a2(string_material)
    a3_vector = a3(string_material)
    return np.dot(np.cross(a1_vector, a2_vector), a3_vector)

def b1(string_material):
    """1st reciprocal vector of the unit cell (metres)"""
    return get_reciprocal_vector(string_material, a2(string_material), a3(string_material))

def b2(string_material):
    """2nd reciprocal vector of the unit cell (metres)"""
    return get_reciprocal_vector(string_material, a3(string_material), a1(string_material))

def b3(string_material):
    """3rd reciprocal vector of the unit cell (metres)"""
    return get_reciprocal_vector(string_material, a1(string_material), a2(string_material))

def d_space(string_material, list_miller):
    """Lattice spacing of a diffraction plane"""
    d_1 = list_miller[0] * b1(string_material)
    d_2 = list_miller[1] * b2(string_material)
    d_3 = list_miller[2] * b3(string_material)
    return 2 * np.pi / np.linalg.norm(d_1 + d_2 + d_3)


def get_chio_data(string_material):
    """Fetch chio (Polarizability) data based on material."""
    if string_material == 'Si':
        return chio_Si_real, chio_Si_imag
    elif string_material == 'Ge':
        return chio_Ge_real, chio_Ge_imag

def chio(energy_kev, string_material):
    """Polarizability in the forward direction (interpolation of the tables)"""
    cho_real_data, cho_imag_data = get_chio_data(string_material)
    
    cho_energy_kev = cho_real_data[:, 0]
    cho_real = cho_real_data[:, 1]
    cho_imag = cho_imag_data[:, 1]

    cho_re = np.interp(energy_kev, cho_energy_kev, cho_real)
    cho_im = np.interp(energy_kev, cho_energy_kev, cho_imag)
    
    return cho_re + 1j * cho_im


def Fo(energy_kev, string_material):
    """ Structure Factor in the forward direction """
    chio_now = chio(energy_kev, string_material)
    wav_len = wavelength(energy_kev)
    numerator = chio_now * np.pi * volume(string_material)
    denominator = electron_radius * wav_len ** 2
    return numerator / denominator


def form_factor_real(string_material, list_miller):
    """ Real part of the form factor.
        (interpolation of the tables)
    """
    mil_if = [abs(item) for item in list_miller]

    if (string_material == 'Si') and (mil_if == [2, 2, 0]):
        return form_factor_real_Si220
    elif (string_material == 'Si') and (mil_if == [1, 1, 1]):
        return form_factor_real_Si111
    elif (string_material == 'Ge') and (mil_if == [2, 2, 0]):
        return form_factor_real_Ge220
    elif (string_material == 'Ge') and (mil_if == [1, 1, 1]):
        return form_factor_real_Ge111
    raise ValueError(f"No real part of the form factor found for material {string_material} and Miller index {list_miller}.")


def form_factor_imag(energy_kev, string_material):
    """ Imaginary part of the form factor
        (interpolation of the tables)
    """
    def interpolate(form_imag):
        return np.interp(energy_kev, form_imag[:, 0], form_imag[:, 1])
        
    if string_material == 'Si':
        return interpolate(form_factor_imag_Si)
    elif string_material == 'Ge':
        return interpolate(form_factor_imag_Ge)
    else: 
        raise ValueError(f"No imaginary part for the form factor found for material {string_material}.")


def form_factor(energy_kev, string_material, list_miller):
    """ Form factor
        (interpolation of the tables)
    """
    f_real = form_factor_real(string_material, list_miller)
    f_imag = form_factor_imag(energy_kev, string_material)
    return f_real + 1j * f_imag


def structure_factor_geometric(string_material, list_miller):
    """ Geometric structure factor """

    # Ensure the material is either Si or Ge
    if string_material not in ['Si', 'Ge']:
        raise ValueError(f"Unsupported material {string_material}. Only Si or Ge is allowed.")

    # Convert the list_miller to a tuple for consistent comparison
    miller_tuple = tuple(list_miller)

    # Check for Miller index and return associated structure factor
    if miller_tuple in [(2,2,0), (-2,-2,0)]:
        return 8
    elif miller_tuple == (1, 1, 1):
        return 4 * (1 - 1j)
    elif miller_tuple == (-1, -1, -1):
        return 4 * (1 + 1j)
    else:
        raise ValueError(f"Unsupported Miller index {miller_tuple} for material {string_material}.")


def Fc(energy_kev, string_material, list_miller):
    """ Structure factor (complex) """
    form_fac = form_factor(energy_kev, string_material, list_miller)
    geom_fac = structure_factor_geometric(string_material, list_miller)
    return form_fac * geom_fac


def Fcplus(energy_kev, string_material, list_miller):
    """ Structure factor (complex) """
    return Fc(energy_kev, string_material, list_miller)


def Fcminus(energy_kev, string_material, list_miller):
    """ Structure factor (complex) """
    mil_neg = [- list_miller[0], - list_miller[1], - list_miller[2]]
    return Fc(energy_kev, string_material, mil_neg)


def chih(energy_kev, string_material, list_miller):
    """ Polarizability in the diffraction direction """
    str_factor = Fc(energy_kev, string_material, list_miller)
    numerator = str_factor * electron_radius * (wavelength(energy_kev)) ** 2
    denominator = np.pi * volume(string_material)
    chih_now = - numerator/denominator
    return chih_now


def chihplus(energy_kev, string_material, list_miller):
    """ Polarizability in the diffraction direction (+ sign) """
    return chih(energy_kev, string_material, list_miller)


def chihminus(energy_kev, string_material, list_miller):
    """ Polarizability in the diffraction direction (- sign) """
    mil_neg = [- list_miller[0], - list_miller[1], - list_miller[2]]
    return chih(energy_kev, string_material, mil_neg)


def beta_crystal(energy_kev, string_material, list_miller):
    """ Crystalline beta function """
    chihpl = chihplus(energy_kev, string_material, list_miller)
    chihmi = chihminus(energy_kev, string_material, list_miller)
    return np.angle(np.sqrt(chihpl * chihmi))


""" HERE START THE CRYSTALLINE FUNCTIONS """


def bragg_angle(energy_kev, string_material, list_miller):
    """  Bragg angle """
    d_now = d_space(string_material, list_miller)
    return np.arcsin(1.0 / (2.0 * d_now) * (h_planck * sci.c) / (energy_kev * 1000.0))


def gamma_h(energy_kev, asymmetry, string_material, list_miller):
    """ Asymmetry exit cosine """
    return np.cos(asymmetry - bragg_angle(energy_kev, string_material, list_miller))


def gamma_o(energy_kev, asymmetry, string_material, list_miller):
    """ Asymmetry entrance cosine """
    return np.cos(asymmetry + bragg_angle(energy_kev, string_material, list_miller))


def gamma(energy_kev, asymmetry, string_material, list_miller):
    """ Asymmetry ratio """
    g_h = gamma_h(energy_kev, asymmetry, string_material, list_miller)
    g_o = gamma_o(energy_kev, asymmetry, string_material, list_miller)
    return g_h/g_o


""" Darwin width and Bragg shift for symmetric and asymmetric reflections
    (entrance and exit)
"""


def delta_o(energy_kev, string_material, list_miller):
    """ Half of the Darwin width for symmetric reflection """
    bragg_ang = bragg_angle(energy_kev, string_material, list_miller)
    numerator = polariz * (electron_radius * (wavelength(energy_kev)) ** 2)
    denominator = np.pi * volume(string_material) * np.sin(2 * bragg_ang)
    Fcpl = Fcplus(energy_kev, string_material, list_miller)
    Fcmi = Fcminus(energy_kev, string_material, list_miller)
    Fc_correction = np.sqrt(Fcpl * Fcmi)
    return numerator / denominator * Fc_correction


def bragg_shift_o(energy_kev, string_material, list_miller):
    """ Bragg shift for symmetric reflection """
    bragg_ang = bragg_angle(energy_kev, string_material, list_miller)
    Fo_now = Fo(energy_kev, string_material)
    numerator = (electron_radius * (wavelength(energy_kev)) ** 2 * abs(Fo_now))
    denominator = 2 * np.pi * volume(string_material) * np.sin(2 * bragg_ang)
    asym_correction = 1 - gamma(energy_kev, -np.pi/2, string_material, list_miller)
    return numerator / denominator * asym_correction


def darwin_o(energy_kev, string_material, list_miller):
    """ Darwin width for symmetric reflection """
    return 2*delta_o(energy_kev, string_material, list_miller)


def delta_os(energy_kev, asymmetry, string_material, list_miller):
    """ Half of the entrance Darwin width for asymmetric reflection -
        no correction for grazing incidence
    """
    bragg_ang = bragg_angle(energy_kev, string_material, list_miller)
    numerator = polariz * (electron_radius * (wavelength(energy_kev)) ** 2)
    denominator = np.pi * volume(string_material) * np.sin(2 * bragg_ang)
    asym_correction = np.sqrt(abs(gamma(energy_kev, asymmetry, string_material, list_miller)))
    Fcpl = Fcplus(energy_kev, string_material, list_miller)
    Fcmi = Fcminus(energy_kev, string_material, list_miller)
    Fc_correction = np.sqrt(Fcpl * Fcmi)
    return numerator / denominator * asym_correction * Fc_correction


def darwin_os(energy_kev, asymmetry, string_material, list_miller):
    """ Entrance Darwin width for asymmetric reflection -
        no correction for grazing incidence
    """
    return 2 * delta_os(energy_kev, asymmetry, string_material, list_miller)


def bragg_shift_os(energy_kev, asymmetry, string_material, list_miller):
    """ Entrance Bragg shift for symmetric reflection -
        no correction for grazing incidence
    """
    bragg_ang = bragg_angle(energy_kev, string_material, list_miller)
    Fo_now = Fo(energy_kev, string_material)
    numerator = electron_radius * (wavelength(energy_kev)) ** 2 * abs(Fo_now)
    denominator = 2 * np.pi * volume(string_material) * np.sin(2 * bragg_ang)
    asym_correction = 1 - gamma(energy_kev, asymmetry, string_material, list_miller)
    return numerator / denominator * asym_correction


def delta_hs(energy_kev, asymmetry, string_material, list_miller):
    """ Half of the exit Darwin width for asymmetric reflection -
        no correction for grazing incidence
    """
    bragg_ang = bragg_angle(energy_kev, string_material, list_miller)
    numerator = polariz * (electron_radius * (wavelength(energy_kev)) ** 2)
    denominator = np.pi * volume(string_material) * np.sin(2 * bragg_ang)
    gam = gamma(energy_kev, asymmetry, string_material, list_miller)
    asym_correction = 1 / (np.sqrt(abs(gam)))
    Fc_correction = np.sqrt(Fcplus * Fcminus)
    return numerator / denominator * asym_correction * Fc_correction


def darwin_hs(energy_kev, asymmetry, string_material, list_miller):
    """ Exit Darwin width for asymmetric reflection -
        no correction for grazing incidence
    """
    return 2 * delta_hs(energy_kev, asymmetry, string_material, list_miller)


def bragg_shift_hs(energy_kev, asymmetry, string_material, list_miller):
    """ Exit Bragg shift for symmetric reflection -
        no correction for grazing incidence
    """
    bragg_ang = bragg_angle(energy_kev, string_material, list_miller)
    Fo_now = Fo(energy_kev, string_material)
    numerator = electron_radius * (wavelength(energy_kev)) ** 2 * abs(Fo_now)
    denominator = 2 * np.pi * volume(string_material) * np.sin(2 * bragg_ang)
    asym_correction = (1 - 1 / gamma(energy_kev, asymmetry, string_material, list_miller))
    return numerator / denominator * asym_correction


"""################ CORRECTIONS FOR GRAZING INCIDENCE #######################"""


def bragg_shift_oc(energy_kev, asymmetry, string_material, list_miller):
    """ Entrance Bragg shift for symmetric reflection -
        WITH correction for grazing incidence
    """
    bragg_ang = bragg_angle(energy_kev, string_material, list_miller)
    g_o = gamma_o(energy_kev, asymmetry, string_material, list_miller)
    g_h = gamma_h(energy_kev, asymmetry, string_material, list_miller)
    chio_now = chio(energy_kev, string_material)
    bragg_corr = chio_now / np.sin(2*bragg_ang)
    chio_corr = (g_o - g_h) * np.sqrt(1 - g_o ** 2) * bragg_corr
    radq = np.sqrt(g_o ** 2 - chio_corr)
    numerator = -g_o + radq
    denominator = np.sqrt(1 - g_o ** 2)
    return numerator / denominator


def bragg_shift_hc(energy_kev, asymmetry, string_material, list_miller):
    """ EXIT Bragg shift for symmetric reflection -
        WITH correction for grazing incidence
    """
    bragg_ang = bragg_angle(energy_kev, string_material, list_miller)
    g_o = gamma_o(energy_kev, asymmetry, string_material, list_miller)
    g_h = gamma_h(energy_kev, asymmetry, string_material, list_miller)
    chio_now = chio(energy_kev, string_material)
    bragg_corr = chio_now / np.sin(2 * bragg_ang)
    chih_corr = (g_o - g_h) * np.sqrt(1 - g_h ** 2) * bragg_corr
    radq = np.sqrt(g_h**2 - chih_corr)
    numerator = (g_h + radq)
    denominator = np.sqrt(1 - g_h ** 2)
    return numerator / denominator


def delta_oc(energy_kev, asymmetry, string_material, list_miller):
    """ Half of the entrance Darwin width for asymmetric reflection -
    WITH correction for grazing incidence.
    """
    bragg_ang = bragg_angle(energy_kev, string_material, list_miller)
    g_h = gamma_h(energy_kev, asymmetry, string_material, list_miller)
    br_shift_oc = bragg_shift_oc(energy_kev, asymmetry, string_material, list_miller)
    delt_o = delta_o(energy_kev, string_material, list_miller)
    bragg_sqrt = (np.sqrt(asymmetry + (np.pi / 2) + bragg_ang))
    numerator = delt_o * np.sqrt(abs(g_h)) * bragg_sqrt
    denominator = (br_shift_oc + (asymmetry + (np.pi / 2) + bragg_ang))
    return numerator / denominator


def delta_hc(energy_kev, asymmetry, string_material, list_miller):
    """ Half of the exit Darwin width for asymmetric reflection -
    WITH correction for grazing incidence.
    """
    bragg_ang = bragg_angle(energy_kev, string_material, list_miller)
    gam = gamma(energy_kev, asymmetry, string_material, list_miller)
    delt_os = delta_os(energy_kev, asymmetry, string_material, list_miller)
    bragg_abs = abs((asymmetry + (np.pi / 2) - bragg_ang))
    numerator = delt_os * ((abs(gam)) ** -1 * bragg_abs)
    asy_bragg_cor = asymmetry + (np.pi / 2) - bragg_ang
    chio_now = chio(energy_kev, string_material)
    tan_cor = asy_bragg_cor * chio_now * np.tan(bragg_ang)
    denominator = np.sqrt(asy_bragg_cor ** 2 - tan_cor - chio_now)
    return numerator / denominator


def darwin_oc(energy_kev, asymmetry, string_material, list_miller):
    """ Darwin width for the entrance beam -
        WITH correction for grazing incidence
    """
    return 2 * delta_oc(energy_kev, asymmetry, string_material, list_miller)


def darwin_hc(energy_kev, asymmetry, string_material, list_miller):
    """ Darwin width for the exit beam -
        WITH correction for grazing incidence
    """
    return 2 * delta_hc(energy_kev, asymmetry, string_material, list_miller)


def b_delta_oc(energy_kev, asymmetry, string_material, list_miller):
    """ Asymmetry factor for the entrance beam -
        WITH correction for grazing incidence
    """
    delta_oc_real = delta_oc(energy_kev, asymmetry, string_material, list_miller).real
    delta_o_real = delta_o(energy_kev, string_material, list_miller).real
    return (delta_oc_real / delta_o_real) ** 2


def width_angle_o(energy_kev, asymmetry, string_material, list_miller):
    """ Width of the rocking curve for the entrance beam -
        WITH correction for grazing incidence
    """
    return 2 * (delta_oc(energy_kev, asymmetry, string_material, list_miller).real)


def width_energy(energy_kev, string_material, list_miller):
    """ Width of the energy curve for the entrance beam
    """
    w_angle = np.absolute(darwin_o(energy_kev, string_material, list_miller))
    br_angle = bragg_angle(energy_kev, string_material, list_miller)
    w_energy = energy_kev * w_angle / np.tan(br_angle)
    return w_energy


def b_delta_hc(energy_kev, asymmetry, string_material, list_miller):
    """ Asymmetry factor for the exit beam -
        WITH correction for grazing incidence
    """
    delta_hc_real = delta_hc(energy_kev, asymmetry, string_material, list_miller).real
    delta_o_real = delta_o(energy_kev, string_material, list_miller).real
    return (delta_hc_real / delta_o_real) ** 2


def width_angle_h(energy_kev, asymmetry, string_material, list_miller):
    """ Width of the rocking curve for the exit beam -
        WITH correction for grazing incidence
    """
    return 2 * (delta_hc(energy_kev, asymmetry, string_material, list_miller).real)



def magnification(energy_kev, asymmetry, string_material, list_miller):
    """ Magnification - WITH correction for grazing incidence """
    delta_oc_real = delta_oc(energy_kev, asymmetry, string_material, list_miller).real
    delta_hc_real = delta_hc(energy_kev, asymmetry, string_material, list_miller).real
    return delta_oc_real / delta_hc_real

"""########################################################################"""


""" A & B crystalline parameters """

def A_crystal(energy_kev, string_material, list_miller):
    """ A crystal parameter """
    return -np.tan(beta_crystal(energy_kev, string_material, list_miller))


def B_crystal(energy_kev, asymmetry, string_material, list_miller):
    """ B crystal parameter """
    gam = gamma(energy_kev, asymmetry, string_material, list_miller)
    chio_now = chio(energy_kev, string_material)
    numerator = (1 - gam) * chio_now.imag
    chihpl = chihplus(energy_kev, string_material, list_miller)
    chihmi = chihminus(energy_kev, string_material, list_miller)
    chi_corr = np.sqrt(chihpl * chihmi)
    beta_now = beta_crystal(energy_kev, string_material, list_miller)
    denominator = 2 * np.cos(beta_now) * np.sqrt(abs(gam)) * polariz * chi_corr
    return numerator / denominator


""" Deviation parameters """


def dev_par_real_o(rocking_angle, asymmetry, energy_kev, string_material, list_miller):
    """ Real part of the deviation parameter for the entrance beam """
    br_shift_oc = bragg_shift_oc(energy_kev, asymmetry, string_material, list_miller)
    wd_angle_o = width_angle_o(energy_kev, asymmetry, string_material, list_miller)
    numerator = rocking_angle * np.pi / (180 * 3600) - br_shift_oc.real
    denominator = wd_angle_o.real / 2
    return numerator / denominator


def dev_par_tot_o(rocking_angle, asymmetry, energy_kev, string_material, list_miller):
    """ Deviation parameter for the entrance beam """
    dev_real_o = dev_par_real_o(rocking_angle, asymmetry, energy_kev, string_material, list_miller)
    A_cry = A_crystal(energy_kev, string_material, list_miller)
    B_cry = B_crystal(energy_kev, asymmetry, string_material, list_miller)
    return dev_real_o * (1 + 1j * A_cry) + 1j * B_cry


def dev_par_real_h(rocking_angle, asymmetry, energy_kev, string_material, list_miller):
    """ Real part of the deviation parameter for the exit beam """
    br_shift_hc = bragg_shift_hc(energy_kev, asymmetry, string_material, list_miller)
    wd_angle_h = width_angle_h(energy_kev, asymmetry, string_material, list_miller)
    numerator = rocking_angle * np.pi / (180 * 3600) - br_shift_hc.real
    denominator = wd_angle_h.real / 2
    return numerator / denominator


def dev_par_tot_h(rocking_angle, asymmetry, energy_kev, string_material, list_miller):
    """ Deviation parameter for the exit beam """
    dev_real_h = dev_par_real_h(rocking_angle, asymmetry, energy_kev, string_material, list_miller)
    A_cry = A_crystal(energy_kev, string_material, list_miller)
    B_cry = B_crystal(energy_kev, asymmetry, string_material, list_miller)
    return dev_real_h * (1 + 1j * A_cry) + 1j * B_cry


""" REFLECTIVITIES AND WAVEFUNCTION """

def reflectivity_entrance(rocking_angle, asymmetry, energy_kev, string_material, list_miller):
    """ Reflectivity for the entrance beam """
    dev_tot_o = dev_par_tot_o(rocking_angle, asymmetry, energy_kev, string_material, list_miller)
    dev_real_o = dev_par_real_o(rocking_angle, asymmetry, energy_kev, string_material, list_miller)
    chihpl = chihplus(energy_kev, string_material, list_miller)
    chihmi = chihminus(energy_kev, string_material, list_miller)
    chi_correction = abs(chihpl / chihmi)
    wave_entrance = dev_tot_o - np.sign(dev_real_o) * np.sqrt(dev_tot_o ** 2 - 1)
    return chi_correction * (abs(wave_entrance)) ** 2


def reflectivity_exit(rocking_angle, asymmetry, energy_kev, string_material, list_miller):
    """ Reflectivity for the exit beam """
    dev_tot_h = dev_par_tot_h(rocking_angle, asymmetry, energy_kev, string_material, list_miller)
    dev_real_h = dev_par_real_h(rocking_angle, asymmetry, energy_kev, string_material, list_miller)
    chihpl = chihplus(energy_kev, string_material, list_miller)
    chihmi = chihminus(energy_kev, string_material, list_miller)
    chi_correction = abs(chihpl / chihmi)
    wave_exit = dev_tot_h - np.sign(dev_real_h) * np.sqrt((dev_tot_h) ** 2 - 1)
    return chi_correction * (abs(wave_exit)) ** 2


""" In the following functions a new variable is introduced:
    pre_mag is the magnification the beam underwent before hitting the crystal
    
pre_mag = False|number --> False: the beam underwent no magnification before the crystal,
                                   which means pre_mag = 1.0
                           number: the beam underwent a magnification of number before
                                   the crystal
"""
    

def wave_crystal_function_entrance(rocking_angle, asymmetry, energy_kev, string_material, 
                                   list_miller, pre_mag):
    """ Crystalline function for Entrance
    
    returns:
        wave_entrance : entrance crystalline wave
        mag_tot : total magnification, which means magnification produced by the current 
                  crystal by the magnification the beam underwent before the crystal
        mag_now = magnification produced by the crystal
    """
    gam = gamma(energy_kev, asymmetry, string_material, list_miller)
    chihpl = chihplus(energy_kev, string_material, list_miller)
    chihmi = chihminus(energy_kev, string_material, list_miller)
    correction_chi = np.sqrt(chihpl*chihmi)/chihmi
    dev_tot_o = dev_par_tot_o(rocking_angle, asymmetry, energy_kev, string_material, list_miller)
    dev_real_o = dev_par_real_o(rocking_angle, asymmetry, energy_kev, string_material, list_miller)
    wave_entrance_dev = dev_tot_o - np.sign(dev_real_o) * np.sqrt((dev_tot_o) ** 2 - 1)
    wave_entrance = correction_chi * wave_entrance_dev
    mag_now = magnification(energy_kev, asymmetry, string_material, list_miller)
    if pre_mag == False:
        pre_mag = 1.0
    else:
        pre_mag = float(pre_mag)
    mag_tot = mag_now * pre_mag
    return [wave_entrance, mag_tot, mag_now]


def wave_crystal_function_exit(rocking_angle, asymmetry, energy_kev, string_material, 
                               list_miller, pre_mag):
    """ Crystalline function for Exit
    
    returns:
        wave_exit : exit crystalline wave
        mag_tot : total magnification, which means magnification produced by the current 
                  crystal by the magnification the beam underwent before the crystal
        mag_now = magnification produced by the crystal    
    """
    gam = gamma(energy_kev, asymmetry, string_material, list_miller)
    chihpl = chihplus(energy_kev, string_material, list_miller)
    chihmi = chihminus(energy_kev, string_material, list_miller)
    correction_chi = np.sqrt(chihpl*chihmi)/chihmi
    dev_tot_h = dev_par_tot_h(rocking_angle * pre_mag, asymmetry, energy_kev, 
                              string_material, list_miller)
    dev_real_h = dev_par_real_h(rocking_angle * pre_mag, asymmetry, energy_kev, 
                                string_material, list_miller)
    wave_exit_dev = dev_tot_h - np.sign(dev_real_h) * np.sqrt((dev_tot_h) ** 2 - 1)
    wave_exit = correction_chi * wave_exit_dev
    mag_now = magnification(energy_kev, asymmetry, string_material, list_miller)
    if pre_mag == False:
        pre_mag = 1.0
    else:
        pre_mag = float(pre_mag)
    mag_tot = mag_now * pre_mag
    return [wave_exit, mag_tot, mag_now]


def shift_energy(energy_kev, range_kev, string_material, list_miller):
    """ Computes the shift from the Bragg angle caused by a variation of the energy
    and it returns an angle in arcsec.

    energy_kev is the energy at which the reflection should occur,
    whose Bragg angle is taken as zero for the rocking angle.
    range_kev is range of energies explored around the central energy, it should be an array.
    """
    bragg_ang = bragg_angle(energy_kev, string_material, list_miller)
    bragg_ang_range = bragg_angle(energy_kev + range_kev, string_material, list_miller)
    shift_rad = bragg_ang - bragg_ang_range
    shift_arcsec = arcsec(shift_rad)
    return shift_arcsec


def wave_crystal_poly_entrance(rocking_angle, misalignement_angle, asymmetry, energy_kev,
                      range_kev, string_material, list_miller, centred, parallel, pre_mag):
    """ Calculates the crystalline function for a polychromatic beam

    rocking_angle is expressed in arcsec
    misalignement_angle is expressed in arcsec
    asymmetry is expressed in radians
        symmetric Bragg reflection is obtained for asymmetry = - Pi / 2
        asymmetry angle for grazing incidence is negative
    energy_kev is the energy at which the reflection should occur,
        whose Bragg angle is taken as zero for the rocking angle.
    range_kev is range of energies explored around the central energy, it should be an array.
        for a monochromatic beam, range_kev = 0
    string_material is the name of the material 'Name' --> for now can be 'Si' or 'Ge' 
    list_miller is the vector of the Miller index [xyz] --> for now can be [220] or [111] 
    centred = True|False --> True the rocking curve is centred at the zero,
                             False the refraction shift is visible
    parallel = True|False --> respect to the first crystals that undergoes diffraction,
                              True: this crystal is in the parallel, non-dispersive geometry
                              False: this crystal is in the antiparallel, dispersive geometry
    pre_mag = False|number --> magnification the beam underwent before hitting the crystal
                               False: the beam underwent no magnification before the crystal,
                                   which means pre_mag = 1.0
                               number: the beam underwent a magnification of number before
                                   the crystal
    
    returns:
        wave_cry_poly : crystalline wave for polychromatic beam
        mag_tot : total magnification, which means magnification produced by the current 
                  crystal by the magnification the beam underwent before the crystal
        mag_now = magnification produced by the crystal     
    """
    
    if centred == True:
        refraction_shift_complex = bragg_shift_oc(energy_kev, asymmetry, string_material, list_miller)
        refraction_shift = arcsec(refraction_shift_complex.real)
        rocking_origin = refraction_shift
    else:
        rocking_origin = 0

    if parallel == True:
        paral = + 1
    else:
        paral = - 1

    shift_en = shift_energy(energy_kev, range_kev, string_material, list_miller)
    rocking_shift = paral * rocking_angle + rocking_origin + misalignement_angle + shift_en
    wave_cry = wave_crystal_function_entrance(rocking_shift, asymmetry, energy_kev, 
                                              string_material, list_miller, pre_mag)
    wave_cry_poly = wave_cry[0]
    mag_tot = wave_cry[1]
    mag_now = wave_cry[2]
    return [wave_cry_poly, mag_tot, mag_now]


def wave_crystal_poly_exit(rocking_angle, misalignement_angle, asymmetry, energy_kev,
                      range_kev, string_material, list_miller, centred, parallel, pre_mag):
    """ Calculates the crystalline function for a polychromatic beam

    rocking_angle is expressed in arcsec
    misalignement_angle is expressed in arcsec
    asymmetry is expressed in radians
        symmetric Bragg reflection is obtained for asymmetry = - Pi / 2
        asymmetry angle for grazing incidence is negative
    energy_kev is the energy at which the reflection should occur,
        whose Bragg angle is taken as zero for the rocking angle.
    range_kev is range of energies explored around the central energy, it should be an array.
        for a polychromatic beam, range_kev = 0
    string_material is the name of the material 'Name' --> for now can be 'Si' or 'Ge' 
    list_miller is the vector of the Miller index [xyz] --> for now can be [220] or [111] 
    centred = True|False --> True the rocking curve is centred at the zero,
                             False the refraction shift is visible
    parallel = True|False --> respect to the first crystals that undergoes diffraction,
                              True: this crystal is in the parallel, non-dispersive geometry
                              False: this crystal is in the antiparallel, dispersive geometry
    pre_mag = False|number --> magnification the beam underwent before hitting the crystal
                               False: the beam underwent no magnification before the crystal,
                                   which means pre_mag = 1.0
                               number: the beam underwent a magnification of number before
                                   the crystal    
    returns:
        wave_cry_poly : crystalline wave for polychromatic beam
        mag_tot : total magnification, which means magnification produced by the current 
                  crystal by the magnification the beam underwent before the crystal
        mag_now = magnification produced by the crystal         
    """
    
    if centred == True:
        refraction_shift_complex = bragg_shift_hc(energy_kev, asymmetry, string_material, list_miller)
        refraction_shift = arcsec(refraction_shift_complex.real) / pre_mag
        rocking_origin = refraction_shift
    else:
        rocking_origin = 0

    if parallel == True:
        paral = + 1
    else:
        paral = - 1

    shift_en = shift_energy(energy_kev, range_kev, string_material, list_miller)
    rocking_shift = paral * rocking_angle + rocking_origin + misalignement_angle + shift_en
    wave_cry = wave_crystal_function_exit(rocking_shift, asymmetry, energy_kev, 
                                          string_material, list_miller, pre_mag)
    wave_cry_poly = wave_cry[0]
    mag_tot = wave_cry[1]
    mag_now = wave_cry[2]
    return [wave_cry_poly, mag_tot, mag_now]


def frequency_to_arcsec(spatial_frequency, energy_kev):
    """ Transform spatial frequencies into a rocking angle in arcsec """
    rock_angle = arcsec(np.arcsin(wavelength(energy_kev) * spatial_frequency))
    return rock_angle


def arcsec_to_frequency(rocking_angle, energy_kev):
    """ Transform a rocking angle in arcsec to spatial frequencies """
    spat_freq = np.sin(rad(rocking_angle)) / wavelength(energy_kev)
    return spat_freq


def compute_crystal_frequency_entrance(spatial_frequency, misalignement_angle, asymmetry, energy_kev,
                      range_kev, string_material, list_miller, centred, parallel, pre_mag, singularity_shift=-0.01):
    """ Computes the crystalline function in spatial frequencies.

    spatial_frequency is expressed in 1 / meter
    misalignement_angle is expressed in arcsec
    asymmetry is expressed in radians
        symmetric Bragg reflection is obtained for asymmetry = - Pi / 2
        asymmetry angle for grazing incidence is negative
    energy_kev is the energy at which the reflection should occur,
        whose Bragg angle is taken as zero for the rocking angle.
    range_kev is range of energies explored around the central energy, it should be an array.
        for a polychromatic beam, range_kev = 0
    string_material is the name of the material 'Name' --> for now can be 'Si' or 'Ge' 
    list_miller is the vector of the Miller index [xyz] --> for now can be [220] or [111] 
    centred = True|False --> True the rocking curve is centred at the zero,
                             False the refraction shift is visible
    parallel = True|False --> respect to the first crystals that undergoes diffraction,
                              True: this crystal is in the parallel, non-dispersive geometry
                              False: this crystal is in the antiparallel, dispersive geometry
    pre_mag = False|number --> magnification the beam underwent before hitting the crystal
                               False: the beam underwent no magnification before the crystal,
                                   which means pre_mag = 1.0
                               number: the beam underwent a magnification of number before
                                   the crystal  
    singularity_shift = shift in the rocking angle to exclude the singularity in zero. This shift is
                        expressed in units of the Darwin width (half of the FWHM). The standard shift
                        is 0.01 Darwin width.                                
    
    returns:
        wave_cry_freq : crystalline wave for polychromatic beam  in spatial frequencies
        mag_tot : total magnification, which means magnification produced by the current 
                  crystal by the magnification the beam underwent before the crystal
        mag_now = magnification produced by the crystal        
    """
    
    rocking_angle = frequency_to_arcsec(spatial_frequency, energy_kev)
    sing_shift = singularity_shift * arcsec(delta_oc(energy_kev, asymmetry, string_material, list_miller).real)
    wave_cry = wave_crystal_poly_entrance(rocking_angle + sing_shift, misalignement_angle, asymmetry, energy_kev, 
                        range_kev, string_material, list_miller, centred, parallel, pre_mag)
    wave_cry_freq = wave_cry[0]
    mag_tot = wave_cry[1]
    mag_now = wave_cry[2]
    return [wave_cry_freq, mag_tot, mag_now]


def compute_crystal_frequency_exit(spatial_frequency, misalignement_angle, asymmetry, energy_kev,
                      range_kev, string_material, list_miller, centred, parallel, pre_mag, singularity_shift=0.01):
    """ Computes the crystalline function in spatial frequencies.

    spatial_frequency is expressed in 1 / meter
    misalignement_angle is expressed in arcsec
    asymmetry is expressed in radians
        symmetric Bragg reflection is obtained for asymmetry = - Pi / 2
        asymmetry angle for grazing incidence is negative
    energy_kev is the energy at which the reflection should occur,
        whose Bragg angle is taken as zero for the rocking angle.
    range_kev is range of energies explored around the central energy, it should be an array.
        for a polychromatic beam, range_kev = 0
    string_material is the name of the material 'Name' --> for now can be 'Si' or 'Ge' 
    list_miller is the vector of the Miller index [xyz] --> for now can be [220] or [111] 
    centred = True|False --> True the rocking curve is centred at the zero,
                             False the refraction shift is visible
    parallel = True|False --> respect to the first crystals that undergoes diffraction,
                              True: this crystal is in the parallel, non-dispersive geometry
                              False: this crystal is in the antiparallel, dispersive geometry
    pre_mag = False|number --> magnification the beam underwent before hitting the crystal
                               False: the beam underwent no magnification before the crystal,
                                   which means pre_mag = 1.0
                               number: the beam underwent a magnification of number before
                                   the crystal    
    singularity_shift = shift in the rocking angle to exclude the singularity in zero. This shift is
                        expressed in units of the Darwin width (half of the FWHM). The standard shift
                        is 0.01 Darwin width.                               
                                   
    returns:
        wave_cry_freq : crystalline wave for polychromatic beam  in spatial frequencies
        mag_tot : total magnification, which means magnification produced by the current 
                  crystal by the magnification the beam underwent before the crystal
        mag_now = magnification produced by the crystal       
    """
    
    rocking_angle = frequency_to_arcsec(spatial_frequency, energy_kev)
    sing_shift = singularity_shift * arcsec(delta_oc(energy_kev, asymmetry, string_material, list_miller).real)
    wave_cry = wave_crystal_poly_exit(rocking_angle + sing_shift, misalignement_angle, asymmetry, energy_kev,
                      range_kev, string_material, list_miller, centred, parallel, pre_mag)
    wave_cry_freq = wave_cry[0]
    mag_tot = wave_cry[1]
    mag_now = wave_cry[2]
    return [wave_cry_freq, mag_tot, mag_now]


def asymmetry_best(energy_kev, string_material, list_miller):
    """ Calculates the asymmetry angle (miscut angle) that maximizes the resolution,
        obtained maximizing the angular acceptance.
    
    input:
        asymmetry = asymmetry angle in Radians:
            for symmetric Bragg: asymmetry = - np.pi / 2
            for symmetric Laue: asymmetry = 0   
        azimuth = azimuth angle in Radians 
            (rotation in the diffraction plane, rotation angle of the Hirano setup):
            for conventional geometry (maximum magnification): azimuth = 0
            for rotation of 90 Degrees (no magnification): azimuth = np.pi / 2
    
    returns:
        Effective asymmetry angle (miscut angle) in Radians
            for symmetric Bragg: asymmetry = - np.pi / 2
            for symmetric Laue: asymmetry = 0
    """
    
    def angle_width_arcsec(asymmetry):
        """ Calculates the angular acceptance of the crystal.
            Maximum resolution is obtained for maximum angular acceptance
        """
        angle_width_rad = width_angle_o(energy_kev, asymmetry, string_material, list_miller)
        return arcsec(angle_width_rad)
    
    starting_point = - np.pi / 2 - bragg_angle(energy_kev, string_material, list_miller) * 999 / 1000
    # Starting point for the optimization of the asymmetry 
    # It should be as near as possible to the asymmetry for 
    # which diffraction changes from Bragg to Laue
        
    best_asymmetry = opt.fmin(lambda asymmetry: - angle_width_arcsec(asymmetry), starting_point)[0]
    # Optimization of the asymmetry for best magnification 
        
    return best_asymmetry
    

def resolution(energy_kev, asymmetry, string_material, list_miller):
    """ Calculates the resolution of the crystal in um.
    """    
    wave_length = wavelength(energy_kev)
    acceptance_angle_o = 2 * width_angle_o(energy_kev, asymmetry, string_material, list_miller)
    resolution_metres = 2 * wave_length / acceptance_angle_o    
    resolution_um = resolution_metres * 10 ** 6
    return resolution_um
    


pass