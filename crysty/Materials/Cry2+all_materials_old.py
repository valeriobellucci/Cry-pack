# Mathematic functions for non-complex numbers.
import math
# Mathematic functions for complex numbers.
import cmath
# NumPy - the fundamental package for scientific computing with Python.
import numpy as np
# Collection of command style functions that make matplotlib work like MATLAB.
import matplotlib.pyplot as plt
# Collection of constants
import scipy.constants as sci


""" Basic physical constants """


electron_radius = sci.physical_constants['classical electron radius'][0]  # 2.81794*10**-15 m
hplanck = sci.physical_constants['Planck constant in eV s'][0]  # 4.1356692*10**-15 eV s
Cpol = 1.0  # Polarization factor now
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


def wavelength(energy_keV):
    """ Retuns the wavelength in metres  """
    return (hplanck * sci.c) / (energy_keV * 1000.0)


def frequency(energy_keV):
    """ Retuns the frequency  """
    return sci.c / wavelength(energy_keV)


def fase(wave_function):
    """ Returns the phase """
    return cmath.atan(wave_function.imag/wave_function.real)
    #return cmath.phase(wave_function)


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
    for i in range(len(matrix)):
        print(matrix[i])

        
def txt_to_list(address):
    with open(address) as f:
        raw_txt = f.read()

    raw_txt = raw_txt.splitlines()

    to_list = [i.split() for i in raw_txt]

    for i in range(len(to_list)):
        for j in range(len(to_list[i])):
            if to_list[i][j].translate(None,'.').isdigit():
                to_list[i][j] = float(to_list[i][j])
    return to_list


costantim = txt_to_list("C:\\Users\\xx9905\\Documents\\archivio txt\\costanti materiali.txt")
unitcell = range(20)
unitcell[3] = txt_to_list("C:\\Users\\xx9905\\Documents\\archivio txt\\unit cell C.txt")
unitcell[6] = txt_to_list("C:\\Users\\xx9905\\Documents\\archivio txt\\unit cell Si.txt")
unitcell[8] = txt_to_list("C:\\Users\\xx9905\\Documents\\archivio txt\\unit cell Cu.txt")
unitcell[9] = txt_to_list("C:\\Users\\xx9905\\Documents\\archivio txt\\unit cell Ge.txt")
unitcell[17] = txt_to_list("C:\\Users\\xx9905\\Documents\\archivio txt\\unit cell quartz.txt")

form_factor_real_Si220 = np.loadtxt("C:\Python27\Materials\\form_factor_real_Si220.txt")
form_factor_real_Si111 = np.loadtxt("C:\Python27\Materials\\form_factor_real_Si111.txt")
form_factor_imag_Si = np.loadtxt("C:\Python27\Materials\\f2_Si_50keV.txt")


absorption_Si = np.loadtxt("C:\\Python27\\Materials\\absorption_vs_energy\\14.dat")
absorption_Ge = np.loadtxt("C:\\Python27\\Materials\\absorption_vs_energy\\32.dat")


""" Crystalline constants """


def density(string_material):
    if string_material == 'Si':
        num_mat = 6  # Silicon is the 6th element in the table
    elif string_material == 'Ge':
        num_mat = 9  # Germanium is the 9th element in the table
    return costantim[num_mat][4]
def electron_density(string_material):
    if string_material == 'Si':
        num_mat = 6  # Silicon is the 6th element in the table
    elif string_material == 'Ge':
        num_mat = 9  # Germanium is the 9th element in the table
    dens = costantim[num_mat][4]
    atomic_number = costantim[num_mat][2]
    mass_number = costantim[num_mat][3]
    return  dens * sci.N_A * (atomic_number / mass_number) * 10 ** 6;
def a1(string_material):
    if string_material == 'Si':
        a1_now = np.array(unitcell[6][1][1:])*10**-10  # Silicon is the 6th element in the table
    return a1_now
def a2(string_material):
    if string_material == 'Si':
        a2_now = np.array(unitcell[6][2][1:])*10**-10  # Silicon is the 6th element in the table
    return a2_now
def a3(string_material):
    if string_material == 'Si':
        a3_now = np.array(unitcell[6][3][1:])*10**-10  # Silicon is the 6th element in the table
    return a3_now
def volume(string_material):
    return np.dot(np.cross(a1(string_material),a2(string_material)),a3(string_material))
def b1(string_material): return 2*math.pi*np.cross(a2(string_material),a3(string_material))/volume(string_material)
def b2(string_material): return 2*math.pi*np.cross(a3(string_material),a1(string_material))/volume(string_material)
def b3(string_material): return 2*math.pi*np.cross(a1(string_material),a2(string_material))/volume(string_material)
def d_space(string_material, list_miller):
    d_1 = list_miller[0]*b1(string_material)
    d_2 = list_miller[0]*b2(string_material)
    d_3 = list_miller[0]*b3(string_material)
    return 2 * math.pi / np.linalg.norm(d_1 + d_2 + d_3)


def absorption(energy_keV, string_material):
    if string_material == 'Si':
        abso = absorption_Si
    elif string_material == 'Ge':
        abso = absorption_Ge
    abso_energy = abso[:,0] * 1000
    abso_y = abso[:,1] * density(string_material) * 10 ** 2
    return np.interp(energy_keV,abso_energy,abso_y)


def alfa_material(energy_keV, string_material):
    epo = sci.epsilon_0
    freq = frequency(energy_keV)
    den_e = electron_density(string_material)
    wav_len = wavelength(energy_keV)
    absor = absorption(energy_keV, string_material)
    freq_mat = math.sqrt((den_e * sci.e ** 2)/(sci.m_e * epo)) /(2 * math.pi)
    den_e_corr = (2* math.pi * freq *(den_e * sci.e ** 2) / epo)
    abs_cor = (1 / (2 * math.pi)) * wav_len * absor
    freq_cor = (2 * math.pi * freq) ** 2
    frequencies_cor = (4 * math.pi ** 2 * sci.m_e * (freq_mat ** 2 - freq ** 2))
    total_cor = ((freq_cor * abs_cor) * (frequencies_cor ** 2 * abs_cor))
    sqrt_cor = den_e_corr ** 2 - 4 * total_cor
    numerator = - den_e_corr + sqrt_cor
    denumerator = (2 * (freq_cor * 1 / (2 * math.pi) * wav_len * absor))
    return numerator / denumerator


def Chio(energy_keV, string_material):
    epo = sci.epsilon_0
    freq = frequency(energy_keV)
    den_e = electron_density(string_material)
    alfa_mat = alfa_material(energy_keV, string_material)
    freq_mat = math.sqrt((den_e * sci.e ** 2)/(sci.m_e * epo)) /(2 * math.pi)
    frequencies_cor = (4 * math.pi ** 2 * sci.m_e * (freq_mat ** 2 - freq ** 2))
    denumerator_imag = - 2 * math.pi * freq * alfa_mat
    denumerator_real = frequencies_cor
    denumerator = denumerator_real + 1j * denumerator_imag
    numerator = (den_e * sci.e ** 2) / epo
    return numerator / denumerator


def Fo(energy_keV, string_material):
    Chio_now = Chio(energy_keV, string_material)
    wav_len = wavelength(energy_keV)
    numerator = Chio_now * math.pi * volume(string_material)
    denumerator = electron_radius * wav_len ** 2
    return numerator / denumerator


def form_factor_real(energy_keV, string_material, list_miller):
    ''' Real part of the form factor. Insert the name of the material in string_material
        and the name of the Miller index in list_miller
    '''
    miller = list_miller
    mil_if = [abs(miller[0]),abs(miller[1]),abs(miller[2])]
    if (string_material == 'Si') and (mil_if == [2,2,0]):
        form_real = form_factor_real_Si220
    elif (string_material == 'Si') and (mil_if == [1,1,1]):
        form_real = form_factor_real_Si111
    return np.interp(energy_keV,form_real[:,0],form_real[:,1])


def form_factor_imag(energy_keV, string_material):
    ''' Imaginary part of the form factor '''
    if string_material == 'Si':
        form_imag = form_factor_imag_Si
    return np.interp(energy_keV,form_imag[:,0],form_imag[:,1])

def form_factor(energy_keV, string_material, list_miller):
    f_real = form_factor_real(energy_keV, string_material, list_miller)
    f_imag = form_factor_imag(energy_keV, string_material)
    return f_real + 1j * f_imag

def structure_factor_geometric(string_material, list_miller):
    if (string_material == 'Si') and list_miller == [2,2,0]:
        str_factor = 8
    elif (string_material == 'Si') and list_miller == [-2,-2,-0]:
        str_factor = 8
    elif (string_material == 'Si') and list_miller == [1,1,1]:
        str_factor = 4 * (1 - 1j)
    elif (string_material == 'Si') and list_miller == [-1,-1,-1]:
        str_factor = 4 * (1 + 1j)
    return str_factor

def Fc(energy_keV, string_material, list_miller):
    ''' Structure factor (complex) '''
    form_fac = form_factor(energy_keV, string_material, list_miller)
    geom_fac = structure_factor_geometric(string_material, list_miller)
    return form_fac * geom_fac

def Fcplus(energy_keV, string_material, list_miller):
    ''' Structure factor (complex) '''
    return Fc(energy_keV, string_material, list_miller)

def Fcminus(energy_keV, string_material, list_miller):
    ''' Structure factor (complex) '''
    mil_neg = [- list_miller[0], - list_miller[1], - list_miller[2]]
    return Fc(energy_keV, string_material, mil_neg)

def Chih(energy_keV, string_material, list_miller):
    str_factor = Fc(energy_keV, string_material, list_miller)
    numerator = str_factor * electron_radius * (wavelength(energy_keV)) ** 2
    denominator = math.pi * volume(string_material)
    return - numerator/denominator

def Chihplus(energy_keV, string_material, list_miller):
    return Chih(energy_keV, string_material, list_miller)

def Chihminus(energy_keV, string_material, list_miller):
    mil_neg = [- list_miller[0], - list_miller[1], - list_miller[2]]
    return Chih(energy_keV, string_material, mil_neg)

def beta_crystal(energy_keV, string_material, list_miller):
    Chihpl = Chihplus(energy_keV, string_material, list_miller)
    Chihmi = Chihminus(energy_keV, string_material, list_miller)
    return cmath.phase(cmath.sqrt(Chihpl * Chihmi))


""" HERE START THE CRYSTALLINE FUNCTIONS """


def bragg_angle(energy_keV, string_material, list_miller):
    """  Bragg angle"""
    d_now = d_space(string_material, list_miller)
    return math.asin(1 / (2 * d_now) * (hplanck * sci.c) / (energy_keV * 1000))


def gamma_h(energy_keV, asymmetry, string_material, list_miller):
    """ Asymmetry exit cosine """
    return math.cos(asymmetry - bragg_angle(energy_keV, string_material, list_miller))


def gamma_o(energy_keV, asymmetry, string_material, list_miller):
    """ Asymmetry entrance cosine """
    return math.cos(asymmetry + bragg_angle(energy_keV, string_material, list_miller))


def gamma(energy_keV, asymmetry, string_material, list_miller):
    """ Asymmetry ratio """
    # function gamma_h
    g_h = math.cos(asymmetry - bragg_angle(energy_keV, string_material, list_miller))
    # function gamma_o
    g_o = math.cos(asymmetry + bragg_angle(energy_keV, string_material, list_miller))
    return g_h/g_o


""" Darwin width and Bragg shift for symmetric and asymmetric reflections
    (entrance and exit) """


def delta_o(energy_keV, string_material, list_miller):
    """ Half of the Darwin width for symmetric reflection """
    bragg_ang = bragg_angle(energy_keV, string_material, list_miller)
    numerator = Cpol * (electron_radius * (wavelength(energy_keV)) ** 2)
    denumerator = math.pi * volume(string_material) * math.sin(2 * bragg_ang)
    Fcpl = Fcplus(energy_keV, string_material, list_miller)
    Fcmi = Fcminus(energy_keV, string_material, list_miller)
    Fc_correction = cmath.sqrt(Fcpl * Fcmi)
    return numerator / denumerator * Fc_correction


def bragg_shift_o(energy_keV, string_material, list_miller):
    """ Bragg shift for symmetric reflection """
    bragg_ang = bragg_angle(energy_keV, string_material, list_miller)
    Fo_now = Fo(energy_keV, string_material)
    numerator = (electron_radius * (wavelength(energy_keV)) ** 2 * abs(Fo_now))
    denumerator = 2 * math.pi * volume(string_material) * math.sin(2 * bragg_ang)
    asym_correction = 1 - gamma(energy_keV, -math.pi/2, string_material, list_miller)
    return numerator / denumerator * asym_correction


def darwin_o(energy_keV, string_material, list_miller):
    """ Darwin width """
    return 2*delta_o(energy_keV, string_material, list_miller)


def delta_os(energy_keV, asymmetry, string_material, list_miller):
    """ Half of the entrance Darwin width for asymmetric reflection -
        no correction for grazing incidence
    """
    bragg_ang = bragg_angle(energy_keV, string_material, list_miller)
    numerator = Cpol * (electron_radius * (wavelength(energy_keV)) ** 2)
    denominator = math.pi * volume(string_material) * math.sin(2 * bragg_ang)
    asym_correction = math.sqrt(abs(gamma(energy_keV, asymmetry, string_material, list_miller)))
    Fcpl = Fcplus(energy_keV, string_material, list_miller)
    Fcmi = Fcminus(energy_keV, string_material, list_miller)
    Fc_correction = cmath.sqrt(Fcpl * Fcmi)
    return numerator / denominator * asym_correction * Fc_correction


def darwin_os(energy_keV, asymmetry, string_material, list_miller):
    """ Entrance Darwin width for asymmetric reflection -
        no correction for grazing incidence
    """
    return 2 * delta_os(energy_keV, asymmetry, string_material, list_miller)


def bragg_shift_os(energy_keV, asymmetry, string_material, list_miller):
    """ Entrance Bragg shift for symmetric reflection -
        no correction for grazing incidence
    """
    bragg_ang = bragg_angle(energy_keV, string_material, list_miller)
    Fo_now = Fo(energy_keV, string_material)
    numerator = electron_radius * (wavelength(energy_keV)) ** 2 * abs(Fo_now)
    denumerator = 2 * math.pi * volume(string_material) * math.sin(2 * bragg_ang)
    asym_correction = 1 - gamma(energy_keV, asymmetry, string_material, list_miller)
    return numerator / denumerator * asym_correction


def delta_hs(energy_keV, asymmetry, string_material, list_miller):
    """ Half of the exit Darwin width for asymmetric reflection -
        no correction for grazing incidence
    """
    bragg_ang = bragg_angle(energy_keV, string_material, list_miller)
    numerator = Cpol * (electron_radius * (wavelength(energy_keV)) ** 2)
    denominator = math.pi * volume(string_material) * math.sin(2 * bragg_ang)
    gam = gamma(energy_keV, asymmetry, string_material, list_miller)
    asym_correction = 1 / (math.sqrt(abs(gam)))
    Fc_correction = cmath.sqrt(Fcplus * Fcminus)
    return numerator / denominator * asym_correction * Fc_correction


def darwin_hs(energy_keV, asymmetry, string_material, list_miller):
    """ Exit Darwin width for asymmetric reflection -
        no correction for grazing incidence
    """
    return 2 * delta_hs(energy_keV, asymmetry, string_material, list_miller)


def bragg_shift_hs(energy_keV, asymmetry, string_material, list_miller):
    """ Exit Bragg shift for symmetric reflection -
        no correction for grazing incidence
    """
    bragg_ang = bragg_angle(energy_keV, string_material, list_miller)
    Fo_now = Fo(energy_keV, string_material)
    numerator = electron_radius * (wavelength(energy_keV)) ** 2 * abs(Fo_now)
    denominator = 2 * math.pi * volume(string_material) * math.sin(2 * bragg_ang)
    asym_correction = (1 - 1 / gamma(energy_keV, asymmetry, string_material, list_miller))
    return numerator / denominator * asym_correction


"""################ CORRECTIONS FOR GRAZING INCIDENCE #######################"""


def bragg_shift_oc(energy_keV, asymmetry, string_material, list_miller):
    """ Entrance Bragg shift for symmetric reflection -
        WITH correction for grazing incidence
    """
    bragg_ang = bragg_angle(energy_keV, string_material, list_miller)
    g_o = gamma_o(energy_keV, asymmetry, string_material, list_miller)
    g_h = gamma_h(energy_keV, asymmetry, string_material, list_miller)
    Chio_now = Chio(energy_keV, string_material)
    bragg_corr = Chio_now / math.sin(2*bragg_ang)
    Chio_corr = (g_o - g_h) * cmath.sqrt(1 - g_o ** 2) * bragg_corr
    radq = cmath.sqrt(g_o ** 2 - Chio_corr)
    numerator = -g_o + radq
    denominator = cmath.sqrt(1 - g_o ** 2)
    return numerator / denominator


def bragg_shift_hc(energy_keV, asymmetry, string_material, list_miller):
    """ EXIT Bragg shift for symmetric reflection -
        WITH correction for grazing incidence
    """
    bragg_ang = bragg_angle(energy_keV, string_material, list_miller)
    g_o = gamma_o(energy_keV, asymmetry, string_material, list_miller)
    g_h = gamma_h(energy_keV, asymmetry, string_material, list_miller)
    Chio_now = Chio(energy_keV, string_material)
    bragg_corr = Chio_now / math.sin(2 * bragg_ang)
    Chih_corr = (g_o - g_h) * cmath.sqrt(1 - g_h ** 2) * bragg_corr
    radq = cmath.sqrt(g_h**2 - Chih_corr)
    numerator = (g_h + radq)
    denominator = cmath.sqrt(1 - g_h ** 2)
    return numerator / denominator


def delta_oc(energy_keV, asymmetry, string_material, list_miller):
    """ Half of the entrance Darwin width for asymmetric reflection -
    WITH correction for grazing incidence.
    """
    bragg_ang = bragg_angle(energy_keV, string_material, list_miller)
    g_h = gamma_h(energy_keV, asymmetry, string_material, list_miller)
    br_shift_oc = bragg_shift_oc(energy_keV, asymmetry, string_material, list_miller)
    delt_o = delta_o(energy_keV, string_material, list_miller)
    bragg_sqrt = (cmath.sqrt(asymmetry + (math.pi / 2) + bragg_ang))
    numerator = delt_o * math.sqrt(abs(g_h)) * bragg_sqrt
    denumerator = (br_shift_oc + (asymmetry + (math.pi / 2) + bragg_ang))
    return numerator / denumerator


def delta_hc(energy_keV, asymmetry, string_material, list_miller):
    """ Half of the exit Darwin width for asymmetric reflection -
    WITH correction for grazing incidence.
    """
    bragg_ang = bragg_angle(energy_keV, string_material, list_miller)
    gam = gamma(energy_keV, asymmetry, string_material, list_miller)
    delt_os = delta_os(energy_keV, asymmetry, string_material, list_miller)
    bragg_abs = abs((asymmetry + (math.pi / 2) - bragg_ang))
    numerator = delt_os * ((abs(gam)) ** -1 * bragg_abs)
    asy_bragg_cor = asymmetry + (math.pi / 2) - bragg_ang
    Chio_now = Chio(energy_keV, string_material)
    tan_cor = asy_bragg_cor * Chio_now * math.tan(bragg_ang)
    denumerator = cmath.sqrt(asy_bragg_cor ** 2 - tan_cor - Chio_now)
    return numerator / denumerator


def darwin_oc(energy_keV, asymmetry, string_material, list_miller):
    """ Darwin width for the entrance beam -
        WITH correction for grazing incidence
    """
    return 2 * delta_oc(energy_keV, asymmetry, string_material, list_miller)


def darwin_hc(energy_keV, asymmetry, string_material, list_miller):
    """ Darwin width for the exit beam -
        WITH correction for grazing incidence
    """
    return 2 * delta_hc(energy_keV, asymmetry, string_material, list_miller)


def b_delta_oc(energy_keV, asymmetry, string_material, list_miller):
    """ Asymmetry factor for the entrance beam -
        WITH correction for grazing incidence
    """
    delta_oc_real = delta_oc(energy_keV, asymmetry, string_material, list_miller).real
    delta_o_real = delta_o(energy_keV, string_material, list_miller).real
    return (delta_oc_real / delta_o_real) ** 2


def width_angle_o(energy_keV, asymmetry, string_material, list_miller):
    """ Width of the rocking curve for the entrance beam -
        WITH correction for grazing incidence
    """
    return 2 * (delta_oc(energy_keV, asymmetry, string_material, list_miller).real)


def b_delta_hc(energy_keV, asymmetry, string_material, list_miller):
    """ Asymmetry factor for the exit beam -
        WITH correction for grazing incidence
    """
    delta_hc_real = delta_hc(energy_keV, asymmetry, string_material, list_miller).real
    delta_o_real = delta_o(energy_keV, string_material, list_miller).real
    return (delta_hc_real / delta_o_real) ** 2


def width_angle_h(energy_keV, asymmetry, string_material, list_miller):
    """ Width of the rocking curve for the exit beam -
        WITH correction for grazing incidence
    """
    return 2 * (delta_hc(energy_keV, asymmetry, string_material, list_miller).real)


def magnification(energy_keV, asymmetry, string_material, list_miller):
    """ Magnification - WITH correction for grazing incidence """
    delta_oc_real = delta_oc(energy_keV, asymmetry, string_material, list_miller).real
    delta_hc_real = delta_hc(energy_keV, asymmetry, string_material, list_miller).real
    return delta_oc_real / delta_hc_real

########################################################################


""" A & B crystalline parameters """


def A_crystal_symmetric(energy_keV, string_material, list_miller):
    """ A crystal parameter - only symmetric reflection """
    Chihpl = Chihplus(energy_keV, string_material, list_miller)
    return -(Chihpl).imag / (Chihpl).real


def B_crystal_symmetric(energy_keV, string_material, list_miller):
    """ B crystal parameter - only symmetric reflection """
    Chio_now = Chio(energy_keV, string_material)
    Chihpl = Chihplus(energy_keV, string_material, list_miller)
    return -(Chio_now).imag / (Chihpl).real


def A_crystal(energy_keV, string_material, list_miller):
    """ A crystal parameter """
    return -cmath.tan(beta_crystal(energy_keV, string_material, list_miller))


def B_crystal(energy_keV, asymmetry, string_material, list_miller):
    """ B crystal parameter """
    gam = gamma(energy_keV, asymmetry, string_material, list_miller)
    Chio_now = Chio(energy_keV, string_material)
    numerator = (1 - gam) * Chio_now.imag
    Chihpl = Chihplus(energy_keV, string_material, list_miller)
    Chihmi = Chihminus(energy_keV, string_material, list_miller)
    Chi_corr = cmath.sqrt(Chihpl * Chihmi)
    beta_now = beta_crystal(energy_keV, string_material, list_miller)
    denumerator = 2 * cmath.cos(beta_now) * math.sqrt(abs(gam)) * Cpol * Chi_corr
    return numerator / denumerator


""" Deviation parameters """


def dev_par_real_o(rocking_angle, asymmetry, energy_keV, string_material, list_miller):
    """ Real part of the deviation parameter for the entrance beam """
    br_shift_oc = bragg_shift_oc(energy_keV, asymmetry, string_material, list_miller)
    wd_angle_o = width_angle_o(energy_keV, asymmetry, string_material, list_miller)
    numerator = rocking_angle * math.pi / (180 * 3600) - br_shift_oc.real
    denominator = wd_angle_o.real / 2
    return numerator / denominator


def dev_par_tot_o(rocking_angle, asymmetry, energy_keV, string_material, list_miller):
    """ Deviation parameter for the entrance beam """
    dev_real_o = dev_par_real_o(rocking_angle, asymmetry, energy_keV, string_material, list_miller)
    A_cry = A_crystal(energy_keV, string_material, list_miller)
    B_cry = B_crystal(energy_keV, asymmetry, string_material, list_miller)
    return dev_real_o * (1 + 1j * A_cry) + 1j * B_cry


def dev_par_real_h(rocking_angle, asymmetry, energy_keV, string_material, list_miller):
    """ Real part of the deviation parameter for the exit beam """
    br_shift_hc = bragg_shift_hc(energy_keV, asymmetry, string_material, list_miller)
    wd_angle_h = width_angle_h(energy_keV, asymmetry, string_material, list_miller)
    numerator = rocking_angle * math.pi / (180 * 3600) - br_shift_hc.real
    denumerator = wd_angle_h.real / 2
    return numerator / denumerator


def dev_par_tot_h(rocking_angle, asymmetry, energy_keV, string_material, list_miller):
    """ Deviation parameter for the exit beam """
    dev_real_h = dev_par_real_h(rocking_angle, asymmetry, energy_keV, string_material, list_miller)
    A_cry = A_crystal(energy_keV, string_material, list_miller)
    B_cry = B_crystal(energy_keV, asymmetry, string_material, list_miller)
    return dev_real_h * (1 + 1j * A_cry) + 1j * B_cry


""" REFLECTIVITIES AND WAVEFUNCTION """


def reflectivity_entrance(rocking_angle, asymmetry, energy_keV, string_material, list_miller):
    """ Reflectivity for the entrance beam """
    dev_tot_o = dev_par_tot_o(rocking_angle, asymmetry, energy_keV, string_material, list_miller)
    dev_real_o = dev_par_real_o(rocking_angle, asymmetry, energy_keV, string_material, list_miller)
    Chihpl = Chihplus(energy_keV, string_material, list_miller)
    Chihmi = Chihminus(energy_keV, string_material, list_miller)
    Chi_correction = abs(Chihpl / Chihmi)
    wave_entrance = dev_tot_o - np.sign(dev_real_o) * cmath.sqrt(dev_tot_o ** 2 - 1)
    return Chi_correction * (abs(wave_entrance)) ** 2


def reflectivity_exit(rocking_angle, asymmetry, energy_keV, string_material, list_miller):
    """ Reflectivity for the exit beam """
    dev_tot_h = dev_par_tot_h(rocking_angle, asymmetry, energy_keV, string_material, list_miller)
    dev_real_h = dev_par_real_h(rocking_angle, asymmetry, energy_keV, string_material, list_miller)
    Chihpl = Chihplus(energy_keV, string_material, list_miller)
    Chihmi = Chihminus(energy_keV, string_material, list_miller)
    Chi_correction = abs(Chihpl / Chihmi)
    wave_exit = dev_tot_h - np.sign(dev_real_h) * cmath.sqrt((dev_tot_h) ** 2 - 1)
    return Chi_correction * (abs(wave_exit)) ** 2


def wave_crystal_function(rocking_angle, asymmetry, energy_keV, string_material, list_miller):
    """ Crystalline function -->
        Exit_wave = wave_crystal_function * Entrance_plane_wave
    """
    gam = gamma(energy_keV, asymmetry, string_material, list_miller)
    Chihpl = Chihplus(energy_keV, string_material, list_miller)
    Chihmi = Chihminus(energy_keV, string_material, list_miller)
    correction_Chi = cmath.sqrt(Chihpl*Chihmi)/Chihmi
    dev_tot_h = dev_par_tot_h(rocking_angle, asymmetry, energy_keV, string_material, list_miller)
    dev_real_h = dev_par_real_h(rocking_angle, asymmetry, energy_keV, string_material, list_miller)
    wave_exit = dev_tot_h - np.sign(dev_real_h) * cmath.sqrt((dev_tot_h) ** 2 - 1)
    return cmath.sqrt(gam) * correction_Chi * wave_exit



'''
minrange = -5
maxrange = 7.5
step = 0.025

plt.plot([rocking_angle for rocking_angle in frange(minrange,maxrange,step)]\
     ,[ \
         fase(wave_crystal_function(rocking_angle, -math.pi/2, 25, 'Si', [2,2,0])) \
       for rocking_angle in frange(minrange,maxrange,step)])

plt.plot([rocking_angle for rocking_angle in frange(minrange,maxrange,step)]\
     ,[ \
         cmath.phase(wave_crystal_function(rocking_angle, -math.pi/2, 25, 'Si', [2,2,0])) \
       for rocking_angle in frange(minrange,maxrange,step)])


plt.plot([rocking_angle for rocking_angle in frange(minrange,maxrange,step)]\
     ,[reflectivity_entrance(rocking_angle, -math.pi/2, 25, 'Si', [2,2,0]) \
       for rocking_angle in frange(minrange,maxrange,step)])

plt.show()
'''

print(Chio(25, 'Si'))
print(Chihplus(25, 'Si', [2,2,0]))
print(Chihminus(25, 'Si', [2,2,0]))
print(beta_crystal(25, 'Si', [2,2,0]))

#print(reflectivity_exit(0, -math.pi/2 - 0.06, 25, 'Si',[1,1,1]))

#print(bragg_shift_oc(25, -math.pi/2 - 0.06, 'Si',[1,1,1]))

