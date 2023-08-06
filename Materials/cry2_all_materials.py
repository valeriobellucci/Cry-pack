""" Importing packages """
# Mathematic functions for non-complex numbers.
import math
# Mathematic functions for complex numbers.
import cmath
# NumPy - the fundamental package for scientific computing with Python.
import numpy as np
# Collection of constants
import scipy.constants as sci
# Collection of command style functions that make matplotlib work like MATLAB.
import matplotlib.pyplot as plt


""" Program for calculating crystalline functions.

This program needs to import some .txt tables containing properties of materials
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
    return radians * (180.0 * 3600.0) / math.pi


def rad(arcseconds):
    """ Transforms arcsec in radians """
    return arcseconds * math.pi / (180.0 * 3600.0)


def wavelength(energy_kev):
    """ Retuns the wavelength in metres  """
    return (h_planck * sci.c) / (energy_kev * 1000.0)


def frequency(energy_kev):
    """ Retuns the frequency  """
    return sci.c / wavelength(energy_kev)


def fase(wave_function):
    """ Returns the phase """
    return cmath.phase(wave_function)


def modulo(wave_function):
    """ Returns the modulus """
    return abs(wave_function)


def intensity(wave_function):
    """ Returns the intensity of a wave """
    return (abs(wave_function)) ** 2


def frange(start, stop, step):
    """ range function for float numbers """
    i = start
    while i < stop:
        yield i
        i += step


def printmatrix(matrix):
    """ Prints a matrix in a matrix-like form """
    for i in matrix:
        print i


def txt_to_list(address):
    """ Import a txt file to a list with string and float elements """
    with open(address) as f_txt:
        raw_txt = f_txt.read()

    raw_txt = raw_txt.splitlines()

    to_list = [i.split() for i in raw_txt]

    for i, _ in enumerate(to_list):
        for j, _ in enumerate(to_list[i]):
            if to_list[i][j].translate(None, '.').isdigit():
                to_list[i][j] = float(to_list[i][j])
    return to_list


""" Import txt files containing crystalline constants """


""" Some materials constants """
costantim = txt_to_list("C:\\Python27\\Materials\\materials constants.txt")


""" Position of the atoms in the unit cell """
unitcell = range(2)
unitcell[0] = txt_to_list("C:\\Python27\\Materials\\unit cell Si.txt")
unitcell[1] = txt_to_list("C:\\Python27\\Materials\\unit cell Ge.txt")


""" Form factor of materials - real and imaginary part """
form_factor_real_Si220 = 8.66557
form_factor_real_Si111 = 10.3811
form_factor_real_Ge220 = 23.7972
form_factor_real_Ge111 = 27.3664
form_factor_imag_Si = np.loadtxt("C:\\Python27\\Materials\\f2_Si_50keV.txt")
form_factor_imag_Ge = np.loadtxt("C:\\Python27\\Materials\\f2_Ge_50keV.txt")


""" Polarizability of materials in the forward direction - real and imaginary part """
chio_Si_real = np.loadtxt("C:\\Python27\\Materials\\chio_Si_real.txt")
chio_Si_imag = np.loadtxt("C:\\Python27\\Materials\\chio_Si_imag.txt")
chio_Ge_real = np.loadtxt("C:\\Python27\\Materials\\chio_Ge_real.txt")
chio_Ge_imag = np.loadtxt("C:\\Python27\\Materials\\chio_Ge_imag.txt")


""" Crystalline constants

    Insert the name of the material 'Name' in string_material
    and the vector of the Miller index [xyz] in list_miller
"""


def density(string_material):
    """ Density of a material in g/cm**3 """
    if string_material == 'Si':
        num_mat = 0  # Silicon is the 0th element in the table
    elif string_material == 'Ge':
        num_mat = 1  # Germanium is the 1st element in the table
    return costantim[num_mat][4]


def electron_density(string_material):
    """ Electron density of a material """
    if string_material == 'Si':
        num_mat = 0  # Silicon is the 0th element in the table
    elif string_material == 'Ge':
        num_mat = 1  # Germanium is the 1st element in the table
    dens = costantim[num_mat][4]
    atomic_number = costantim[num_mat][2]
    mass_number = costantim[num_mat][3]
    return dens * sci.N_A * (atomic_number / mass_number) * 10 ** 6


def a1(string_material):
    """ 1st vector of the unit cell (metres) """
    if string_material == 'Si':
        # Silicon is the 0th element in the table
        a1_now = np.array(unitcell[0][1][1:]) * 10 ** -10
    elif string_material == 'Ge':
        # Germanium is the 1st element in the table
        a1_now = np.array(unitcell[1][1][1:]) * 10 ** -10
    return a1_now


def a2(string_material):
    """ 2nd vector of the unit cell (metres) """
    if string_material == 'Si':
        # Silicon is the 0th element in the table
        a2_now = np.array(unitcell[0][2][1:]) * 10 ** -10
    elif string_material == 'Ge':
        # Germanium is the 1st element in the table
        a2_now = np.array(unitcell[1][2][1:]) * 10 ** -10
    return a2_now


def a3(string_material):
    """ 3rd vector of the unit cell (metres) """
    if string_material == 'Si':
        # Silicon is the 0th element in the table
        a3_now = np.array(unitcell[0][3][1:]) * 10 ** -10
    if string_material == 'Ge':
        # Germanium is the 1st element in the table
        a3_now = np.array(unitcell[1][3][1:]) * 10 ** -10
    return a3_now


def volume(string_material):
    """ Volume of the unit cell (metres) """
    cross_product_a1_a2 = np.cross(a1(string_material), a2(string_material))
    dot_product_a3 = np.dot(cross_product_a1_a2, a3(string_material))
    return dot_product_a3


def b1(string_material):
    """ 1st reciprocal vector of the unit cell (metres) """
    cross_product_b1 = np.cross(a2(string_material), a3(string_material))
    return 2 * math.pi * cross_product_b1 / volume(string_material)


def b2(string_material):
    """ 2nd reciprocal vector of the unit cell (metres) """
    cross_product_b2 = np.cross(a3(string_material), a1(string_material))
    return 2 * math.pi * cross_product_b2 / volume(string_material)


def b3(string_material):
    """ 3rd reciprocal vector of the unit cell (metres) """
    cross_product_b3 = np.cross(a1(string_material), a2(string_material))
    return 2 * math.pi * cross_product_b3 / volume(string_material)


def d_space(string_material, list_miller):
    """ Lattice spacing of a diffraction plane """
    d_1 = list_miller[0] * b1(string_material)
    d_2 = list_miller[1] * b2(string_material)
    d_3 = list_miller[2] * b3(string_material)
    return 2 * math.pi / np.linalg.norm(d_1 + d_2 + d_3)


def chio(energy_kev, string_material):
    """ Polarizability in the forward direction
        (interpolation of the tables)
    """
    if string_material == 'Si':
        cho_energy_kev = chio_Si_real[:, 0]
        cho_real = chio_Si_real[:, 1]
        cho_imag = chio_Si_imag[:, 1]
    elif string_material == 'Ge':
        cho_energy_kev = chio_Ge_real[:, 0]
        cho_real = chio_Ge_real[:, 1]
        cho_imag = chio_Ge_imag[:, 1]

    cho_re = np.interp(energy_kev, cho_energy_kev, cho_real)
    cho_im = np.interp(energy_kev, cho_energy_kev, cho_imag)
    cho = cho_re + 1j * cho_im
    return cho


def Fo(energy_kev, string_material):
    """ Structure Factor in the forward direction """
    chio_now = chio(energy_kev, string_material)
    wav_len = wavelength(energy_kev)
    numerator = chio_now * math.pi * volume(string_material)
    denumerator = electron_radius * wav_len ** 2
    return numerator / denumerator


def form_factor_real(string_material, list_miller):
    """ Real part of the form factor.
        (interpolation of the tables)
    """
    miller = list_miller
    mil_if = [abs(miller[0]), abs(miller[1]), abs(miller[2])]
    if (string_material == 'Si') and (mil_if == [2, 2, 0]):
        form_real = form_factor_real_Si220
    elif (string_material == 'Si') and (mil_if == [1, 1, 1]):
        form_real = form_factor_real_Si111
    elif (string_material == 'Ge') and (mil_if == [2, 2, 0]):
        form_real = form_factor_real_Ge220
    elif (string_material == 'Ge') and (mil_if == [1, 1, 1]):
        form_real = form_factor_real_Ge111
    return form_real


def form_factor_imag(energy_kev, string_material):
    """ Imaginary part of the form factor
        (interpolation of the tables)
    """
    if string_material == 'Si':
        form_imag = form_factor_imag_Si
    elif string_material == 'Ge':
        form_imag = form_factor_imag_Ge
    return np.interp(energy_kev, form_imag[:, 0], form_imag[:, 1])


def form_factor(energy_kev, string_material, list_miller):
    """ Form factor
        (interpolation of the tables)
    """
    f_real = form_factor_real(string_material, list_miller)
    f_imag = form_factor_imag(energy_kev, string_material)
    return f_real + 1j * f_imag


def structure_factor_geometric(string_material, list_miller):
    """ Geometric structure factor """
    if (((string_material == 'Si') ^ (string_material == 'Ge'))
            and list_miller == [2, 2, 0]):
        str_factor = 8
    elif (((string_material == 'Si') ^ (string_material == 'Ge'))
          and list_miller == [-2, -2, -0]):
        str_factor = 8
    elif (((string_material == 'Si') ^ (string_material == 'Ge'))
          and list_miller == [1, 1, 1]):
        str_factor = 4 * (1 - 1j)
    elif (((string_material == 'Si') ^ (string_material == 'Ge'))
          and list_miller == [-1, -1, -1]):
        str_factor = 4 * (1 + 1j)
    return str_factor


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
    denominator = math.pi * volume(string_material)
    return - numerator/denominator


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
    return cmath.phase(cmath.sqrt(chihpl * chihmi))


""" HERE START THE CRYSTALLINE FUNCTIONS """


def bragg_angle(energy_kev, string_material, list_miller):
    """  Bragg angle """
    d_now = d_space(string_material, list_miller)
    return math.asin(1 / (2 * d_now) * (h_planck * sci.c) / (energy_kev * 1000))


def gamma_h(energy_kev, asymmetry, string_material, list_miller):
    """ Asymmetry exit cosine """
    return math.cos(asymmetry - bragg_angle(energy_kev, string_material, list_miller))


def gamma_o(energy_kev, asymmetry, string_material, list_miller):
    """ Asymmetry entrance cosine """
    return math.cos(asymmetry + bragg_angle(energy_kev, string_material, list_miller))


def gamma(energy_kev, asymmetry, string_material, list_miller):
    """ Asymmetry ratio """
    # function gamma_h
    g_h = math.cos(asymmetry - bragg_angle(energy_kev, string_material, list_miller))
    # function gamma_o
    g_o = math.cos(asymmetry + bragg_angle(energy_kev, string_material, list_miller))
    return g_h/g_o


""" Darwin width and Bragg shift for symmetric and asymmetric reflections
    (entrance and exit)
"""


def delta_o(energy_kev, string_material, list_miller):
    """ Half of the Darwin width for symmetric reflection """
    bragg_ang = bragg_angle(energy_kev, string_material, list_miller)
    numerator = polariz * (electron_radius * (wavelength(energy_kev)) ** 2)
    denumerator = math.pi * volume(string_material) * math.sin(2 * bragg_ang)
    Fcpl = Fcplus(energy_kev, string_material, list_miller)
    Fcmi = Fcminus(energy_kev, string_material, list_miller)
    Fc_correction = cmath.sqrt(Fcpl * Fcmi)
    return numerator / denumerator * Fc_correction


def bragg_shift_o(energy_kev, string_material, list_miller):
    """ Bragg shift for symmetric reflection """
    bragg_ang = bragg_angle(energy_kev, string_material, list_miller)
    Fo_now = Fo(energy_kev, string_material)
    numerator = (electron_radius * (wavelength(energy_kev)) ** 2 * abs(Fo_now))
    denumerator = 2 * math.pi * volume(string_material) * math.sin(2 * bragg_ang)
    asym_correction = 1 - gamma(energy_kev, -math.pi/2, string_material, list_miller)
    return numerator / denumerator * asym_correction


def darwin_o(energy_kev, string_material, list_miller):
    """ Darwin width for symmetric reflection """
    return 2*delta_o(energy_kev, string_material, list_miller)


def delta_os(energy_kev, asymmetry, string_material, list_miller):
    """ Half of the entrance Darwin width for asymmetric reflection -
        no correction for grazing incidence
    """
    bragg_ang = bragg_angle(energy_kev, string_material, list_miller)
    numerator = polariz * (electron_radius * (wavelength(energy_kev)) ** 2)
    denominator = math.pi * volume(string_material) * math.sin(2 * bragg_ang)
    asym_correction = math.sqrt(abs(gamma(energy_kev, asymmetry, string_material, list_miller)))
    Fcpl = Fcplus(energy_kev, string_material, list_miller)
    Fcmi = Fcminus(energy_kev, string_material, list_miller)
    Fc_correction = cmath.sqrt(Fcpl * Fcmi)
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
    denumerator = 2 * math.pi * volume(string_material) * math.sin(2 * bragg_ang)
    asym_correction = 1 - gamma(energy_kev, asymmetry, string_material, list_miller)
    return numerator / denumerator * asym_correction


def delta_hs(energy_kev, asymmetry, string_material, list_miller):
    """ Half of the exit Darwin width for asymmetric reflection -
        no correction for grazing incidence
    """
    bragg_ang = bragg_angle(energy_kev, string_material, list_miller)
    numerator = polariz * (electron_radius * (wavelength(energy_kev)) ** 2)
    denominator = math.pi * volume(string_material) * math.sin(2 * bragg_ang)
    gam = gamma(energy_kev, asymmetry, string_material, list_miller)
    asym_correction = 1 / (math.sqrt(abs(gam)))
    Fc_correction = cmath.sqrt(Fcplus * Fcminus)
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
    denominator = 2 * math.pi * volume(string_material) * math.sin(2 * bragg_ang)
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
    bragg_corr = chio_now / math.sin(2*bragg_ang)
    chio_corr = (g_o - g_h) * cmath.sqrt(1 - g_o ** 2) * bragg_corr
    radq = cmath.sqrt(g_o ** 2 - chio_corr)
    numerator = -g_o + radq
    denominator = cmath.sqrt(1 - g_o ** 2)
    return numerator / denominator


def bragg_shift_hc(energy_kev, asymmetry, string_material, list_miller):
    """ EXIT Bragg shift for symmetric reflection -
        WITH correction for grazing incidence
    """
    bragg_ang = bragg_angle(energy_kev, string_material, list_miller)
    g_o = gamma_o(energy_kev, asymmetry, string_material, list_miller)
    g_h = gamma_h(energy_kev, asymmetry, string_material, list_miller)
    chio_now = chio(energy_kev, string_material)
    bragg_corr = chio_now / math.sin(2 * bragg_ang)
    chih_corr = (g_o - g_h) * cmath.sqrt(1 - g_h ** 2) * bragg_corr
    radq = cmath.sqrt(g_h**2 - chih_corr)
    numerator = (g_h + radq)
    denominator = cmath.sqrt(1 - g_h ** 2)
    return numerator / denominator


def delta_oc(energy_kev, asymmetry, string_material, list_miller):
    """ Half of the entrance Darwin width for asymmetric reflection -
    WITH correction for grazing incidence.
    """
    bragg_ang = bragg_angle(energy_kev, string_material, list_miller)
    g_h = gamma_h(energy_kev, asymmetry, string_material, list_miller)
    br_shift_oc = bragg_shift_oc(energy_kev, asymmetry, string_material, list_miller)
    delt_o = delta_o(energy_kev, string_material, list_miller)
    bragg_sqrt = (cmath.sqrt(asymmetry + (math.pi / 2) + bragg_ang))
    numerator = delt_o * math.sqrt(abs(g_h)) * bragg_sqrt
    denumerator = (br_shift_oc + (asymmetry + (math.pi / 2) + bragg_ang))
    return numerator / denumerator


def delta_hc(energy_kev, asymmetry, string_material, list_miller):
    """ Half of the exit Darwin width for asymmetric reflection -
    WITH correction for grazing incidence.
    """
    bragg_ang = bragg_angle(energy_kev, string_material, list_miller)
    gam = gamma(energy_kev, asymmetry, string_material, list_miller)
    delt_os = delta_os(energy_kev, asymmetry, string_material, list_miller)
    bragg_abs = abs((asymmetry + (math.pi / 2) - bragg_ang))
    numerator = delt_os * ((abs(gam)) ** -1 * bragg_abs)
    asy_bragg_cor = asymmetry + (math.pi / 2) - bragg_ang
    chio_now = chio(energy_kev, string_material)
    tan_cor = asy_bragg_cor * chio_now * math.tan(bragg_ang)
    denumerator = cmath.sqrt(asy_bragg_cor ** 2 - tan_cor - chio_now)
    return numerator / denumerator


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


def A_crystal_symmetric(energy_kev, string_material, list_miller):
    """ A crystal parameter - only symmetric reflection """
    chihpl = chihplus(energy_kev, string_material, list_miller)
    return -(chihpl).imag / (chihpl).real


def B_crystal_symmetric(energy_kev, string_material, list_miller):
    """ B crystal parameter - only symmetric reflection """
    chio_now = chio(energy_kev, string_material)
    chihpl = chihplus(energy_kev, string_material, list_miller)
    return -(chio_now).imag / (chihpl).real


def A_crystal(energy_kev, string_material, list_miller):
    """ A crystal parameter """
    return -cmath.tan(beta_crystal(energy_kev, string_material, list_miller))


def B_crystal(energy_kev, asymmetry, string_material, list_miller):
    """ B crystal parameter """
    gam = gamma(energy_kev, asymmetry, string_material, list_miller)
    chio_now = chio(energy_kev, string_material)
    numerator = (1 - gam) * chio_now.imag
    chihpl = chihplus(energy_kev, string_material, list_miller)
    chihmi = chihminus(energy_kev, string_material, list_miller)
    chi_corr = cmath.sqrt(chihpl * chihmi)
    beta_now = beta_crystal(energy_kev, string_material, list_miller)
    denumerator = 2 * cmath.cos(beta_now) * math.sqrt(abs(gam)) * polariz * chi_corr
    return numerator / denumerator


""" Deviation parameters """


def dev_par_real_o(rocking_angle, asymmetry, energy_kev, string_material, list_miller):
    """ Real part of the deviation parameter for the entrance beam """
    br_shift_oc = bragg_shift_oc(energy_kev, asymmetry, string_material, list_miller)
    wd_angle_o = width_angle_o(energy_kev, asymmetry, string_material, list_miller)
    numerator = rocking_angle * math.pi / (180 * 3600) - br_shift_oc.real
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
    numerator = rocking_angle * math.pi / (180 * 3600) - br_shift_hc.real
    denumerator = wd_angle_h.real / 2
    return numerator / denumerator


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
    wave_entrance = dev_tot_o - np.sign(dev_real_o) * cmath.sqrt(dev_tot_o ** 2 - 1)
    return chi_correction * (abs(wave_entrance)) ** 2


def reflectivity_exit(rocking_angle, asymmetry, energy_kev, string_material, list_miller):
    """ Reflectivity for the exit beam """
    dev_tot_h = dev_par_tot_h(rocking_angle, asymmetry, energy_kev, string_material, list_miller)
    dev_real_h = dev_par_real_h(rocking_angle, asymmetry, energy_kev, string_material, list_miller)
    chihpl = chihplus(energy_kev, string_material, list_miller)
    chihmi = chihminus(energy_kev, string_material, list_miller)
    chi_correction = abs(chihpl / chihmi)
    wave_exit = dev_tot_h - np.sign(dev_real_h) * cmath.sqrt((dev_tot_h) ** 2 - 1)
    return chi_correction * (abs(wave_exit)) ** 2


def wave_crystal_function(rocking_angle, asymmetry, energy_kev, string_material, list_miller):
    """ Crystalline function -->
        Exit_wave = wave_crystal_function * Entrance_plane_wave
    """
    gam = gamma(energy_kev, asymmetry, string_material, list_miller)
    chihpl = chihplus(energy_kev, string_material, list_miller)
    chihmi = chihminus(energy_kev, string_material, list_miller)
    correction_chi = cmath.sqrt(chihpl*chihmi)/chihmi
    dev_tot_h = dev_par_tot_h(rocking_angle, asymmetry, energy_kev, string_material, list_miller)
    dev_real_h = dev_par_real_h(rocking_angle, asymmetry, energy_kev, string_material, list_miller)
    wave_exit = dev_tot_h - np.sign(dev_real_h) * cmath.sqrt((dev_tot_h) ** 2 - 1)
    return cmath.sqrt(gam) * correction_chi * wave_exit




minrange = -5
maxrange = 15.5
step = 0.025
mnow = 'Si'
asy_now = + 4 * math.pi / 180
mill_ind = [2,2,0]

plt.plot([rocking_angle for rocking_angle in frange(minrange,maxrange,step)]\
     ,[ \
         cmath.phase(wave_crystal_function(rocking_angle, -math.pi/2 + asy_now, 25, mnow, mill_ind)) \
       for rocking_angle in frange(minrange,maxrange,step)])


plt.plot([rocking_angle for rocking_angle in frange(minrange,maxrange,step)]\
     ,[reflectivity_exit(rocking_angle, -math.pi/2 + asy_now, 25, mnow, mill_ind) \
       for rocking_angle in frange(minrange,maxrange,step)])

plt.show()


print(chio(25, mnow))
print(chihplus(25, mnow, [2,2,0]))
print(chihminus(25, mnow, [2,2,0]))
print(beta_crystal(25, mnow, [2,2,0]))

#print(reflectivity_exit(0, -math.pi/2 - 0.06, 25, 'Si',[1,1,1]))

#print(bragg_shift_oc(25, -math.pi/2 - 0.06, 'Si',[1,1,1]))

