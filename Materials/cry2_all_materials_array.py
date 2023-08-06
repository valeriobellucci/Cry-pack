""" Importing packages """
# NumPy - the fundamental package for scientific computing with Python.
import numpy as np
# Collection of constants
import scipy.constants as sci


""" Program for calculating crystalline functions.

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
    return (h_planck * sci.c) / (energy_kev * 1000.0)


def frequency(energy_kev):
    """ Retuns the frequency  """
    return sci.c / wavelength(energy_kev)


def fase(wave_function):
    """ Returns the phase """
    return np.angle(wave_function)


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
costantim = txt_to_list("/home/ws/xx9905/dev/Materials/materials constants.txt")


""" Position of the atoms in the unit cell """
unitcell = range(2)
unitcell[0] = txt_to_list("/home/ws/xx9905/dev/Materials/unit cell Si.txt")
unitcell[1] = txt_to_list("/home/ws/xx9905/dev/Materials/unit cell Ge.txt")


""" Form factor of materials - real and imaginary part """
form_factor_real_Si220 = 8.66557
form_factor_real_Si111 = 10.3811
form_factor_real_Ge220 = 23.7972
form_factor_real_Ge111 = 27.3664
form_factor_imag_Si = np.loadtxt("/home/ws/xx9905/dev/Materials/f2_Si_50keV.txt")
form_factor_imag_Ge = np.loadtxt("/home/ws/xx9905/dev/Materials/f2_Ge_50keV.txt")


""" Polarizability of materials in the forward direction - real and imaginary part """
chio_Si_real = np.loadtxt("/home/ws/xx9905/dev/Materials/chio_Si_real.txt")
chio_Si_imag = np.loadtxt("/home/ws/xx9905/dev/Materials/chio_Si_imag.txt")
chio_Ge_real = np.loadtxt("/home/ws/xx9905/dev/Materials/chio_Ge_real.txt")
chio_Ge_imag = np.loadtxt("/home/ws/xx9905/dev/Materials/chio_Ge_imag.txt")


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
    return 2 * np.pi * cross_product_b1 / volume(string_material)


def b2(string_material):
    """ 2nd reciprocal vector of the unit cell (metres) """
    cross_product_b2 = np.cross(a3(string_material), a1(string_material))
    return 2 * np.pi * cross_product_b2 / volume(string_material)


def b3(string_material):
    """ 3rd reciprocal vector of the unit cell (metres) """
    cross_product_b3 = np.cross(a1(string_material), a2(string_material))
    return 2 * np.pi * cross_product_b3 / volume(string_material)


def d_space(string_material, list_miller):
    """ Lattice spacing of a diffraction plane """
    d_1 = list_miller[0] * b1(string_material)
    d_2 = list_miller[1] * b2(string_material)
    d_3 = list_miller[2] * b3(string_material)
    return 2 * np.pi / np.linalg.norm(d_1 + d_2 + d_3)


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
    '''cho = -0.511*0.0001 + 1j * 0.26386*0.000001'''
    return cho


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
    denominator = np.pi * volume(string_material)
    chih_now = - numerator/denominator
    '''chih_now = 0.105*0.0001 + 1j*0.255*0.000001'''
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
    return np.arcsin(1 / (2 * d_now) * (h_planck * sci.c) / (energy_kev * 1000))


def gamma_h(energy_kev, asymmetry, string_material, list_miller):
    """ Asymmetry exit cosine """
    return np.cos(asymmetry - bragg_angle(energy_kev, string_material, list_miller))


def gamma_o(energy_kev, asymmetry, string_material, list_miller):
    """ Asymmetry entrance cosine """
    return np.cos(asymmetry + bragg_angle(energy_kev, string_material, list_miller))


def gamma(energy_kev, asymmetry, string_material, list_miller):
    """ Asymmetry ratio """
    # function gamma_h
    g_h = np.cos(asymmetry - bragg_angle(energy_kev, string_material, list_miller))
    # function gamma_o
    g_o = np.cos(asymmetry + bragg_angle(energy_kev, string_material, list_miller))
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


def wave_crystal_function(rocking_angle, asymmetry, energy_kev, string_material, list_miller):
    """ Crystalline function -->
        Exit_wave = wave_crystal_function * Entrance_plane_wave
    """
    gam = gamma(energy_kev, asymmetry, string_material, list_miller)
    chihpl = chihplus(energy_kev, string_material, list_miller)
    chihmi = chihminus(energy_kev, string_material, list_miller)
    correction_chi = np.sqrt(chihpl*chihmi)/chihmi
    dev_tot_h = dev_par_tot_h(rocking_angle, asymmetry, energy_kev, string_material, list_miller)
    dev_real_h = dev_par_real_h(rocking_angle, asymmetry, energy_kev, string_material, list_miller)
    wave_exit = dev_tot_h - np.sign(dev_real_h) * np.sqrt((dev_tot_h) ** 2 - 1)
    return correction_chi * wave_exit


def shift_energy(energy_kev, range_kev, string_material, list_miller):
    """ Computes the shift from the Bragg angle provoked by a variation of the energy
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


def wave_crystal_poly(rocking_angle, misalignement_angle, asymmetry, energy_kev,
                      range_kev, string_material, list_miller, centred, parallel):
    """ Calculates the crystalline function for a polychromatic beam

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
    """
    
    if centred == True:
        refraction_shift_complex = bragg_shift_hc(energy_kev, asymmetry, string_material, list_miller)
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
    wave_cry_poly = wave_crystal_function(rocking_shift, asymmetry,
                                          energy_kev, string_material, list_miller)
    return wave_cry_poly


def frequency_to_arcsec(spatial_frequency, energy_kev):
    """ Transform spatial frequencies into a rocking angle in arcsec """
    rock_angle = arcsec(np.arcsin(wavelength(energy_kev) * spatial_frequency))
    return rock_angle


def arcsec_to_frequency(rocking_angle, energy_kev):
    """ Transform a rocking angle in arcsec to spatial frequencies """
    spat_freq = np.sin(rad(rocking_angle)) / wavelength(energy_kev)
    return spat_freq


def compute_crystal_frequency(spatial_frequency, misalignement_angle, asymmetry, energy_kev,
                      range_kev, string_material, list_miller, centred, parallel):
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
    """
    rocking_angle = frequency_to_arcsec(spatial_frequency, energy_kev)
    wave_cry_freq = wave_crystal_poly(rocking_angle, misalignement_angle, asymmetry, energy_kev,
                      range_kev, string_material, list_miller, centred, parallel)
    return wave_cry_freq



'''
# Load this package to initialize the plots

# Collection of command style functions that make matplotlib work like MATLAB.
import matplotlib.pyplot as plt

mnow = 'Si'
mill_ind = [2,2,0]
energy_keV_now = 25
asy_now = - (np.pi / 2) - 4.19 * np.pi / 180
step = 0.025
maxrange = 2 * arcsec(width_angle_h(energy_keV_now, asy_now, mnow, mill_ind))
minrange = - maxrange - step/2
maxrange_freq = arcsec_to_frequency(maxrange, energy_keV_now)
minrange_freq = arcsec_to_frequency(minrange, energy_keV_now)
step_freq = maxrange_freq / 100
maxrange_keV = 0.002
step_keV = maxrange_keV / 200
minrange_keV = - maxrange_keV - step_keV/2

mis_angle_now = 0
range_kev_now = 0
rock_angl_now = 0.1
cent_now = False
paral_now = True
'''


'''
# Plots intensity and phase vs rocking angle

vector_rocking_angle = np.arange(minrange,maxrange,step)
vector_cry = wave_crystal_function(vector_rocking_angle, asy_now,
                                   energy_keV_now, mnow, mill_ind)
vector_cry_poly = wave_crystal_poly(vector_rocking_angle, mis_angle_now, asy_now, energy_keV_now,
                               range_kev_now, mnow, mill_ind, cent_now, paral_now)

plt.plot(vector_rocking_angle, intensity(vector_cry_poly))
plt.plot(vector_rocking_angle, np.angle(vector_cry_poly))
plt.show()
'''

'''
# Plots intensity and phase vs energy 

vector_range_kev = np.arange(-0.002,0.002,0.002/200)
vector_cry_poly = wave_crystal_poly(rock_angl_now, mis_angle_now, asy_now, energy_keV_now,
                               vector_range_kev, mnow, mill_ind, cent_now, paral_now)

plt.plot(vector_range_kev, intensity(vector_cry_poly))
plt.plot(vector_range_kev, np.angle(vector_cry_poly))
plt.show()
'''

'''
# Plots intensity and phase vs spatial frequency

vector_freq = np.arange(minrange_freq,maxrange_freq,step_freq)
vector_cry_freq = compute_crystal_frequency(vector_freq, mis_angle_now, asy_now, energy_keV_now,
                               range_kev_now, mnow, mill_ind, cent_now, paral_now)

plt.plot(vector_freq, intensity(vector_cry_freq))
plt.plot(vector_freq, np.angle(vector_cry_freq))
plt.show()

print maxrange_freq
print minrange_freq
print maxrange
print minrange
'''


# Load this package to initialize the plots

# Collection of command style functions that make matplotlib work like MATLAB.
import matplotlib.pyplot as plt

mnow = 'Ge'
mill_ind = [2,2,0]
energy_keV_now = 10.76
asy_now = - (np.pi / 2) # 0.12287  0.28623
maxrange = 1 * arcsec(width_angle_h(energy_keV_now, asy_now, mnow, mill_ind))
step = maxrange / 200
minrange = - maxrange - step/2
maxrange_freq = arcsec_to_frequency(maxrange, energy_keV_now)
minrange_freq = arcsec_to_frequency(minrange, energy_keV_now)
step_freq = maxrange_freq / 100
maxrange_keV = 0.002
step_keV = maxrange_keV / 200
minrange_keV = - maxrange_keV - step_keV/2

mis_angle_now = 0
range_kev_now = 0
rock_angl_now = 0.1
cent_now = True
paral_now = True

'''
# Plots intensity and phase vs rocking angle

vector_freq = np.arange(-1681818,1678533,3285)
vector_rocking_angle = np.arange(minrange,maxrange,step)

vector_cry_freq = compute_crystal_frequency(vector_freq, mis_angle_now, asy_now, energy_keV_now,
                               range_kev_now, mnow, mill_ind, cent_now, paral_now)
#vector_cry_rock = wave_crystal_poly(vector_rocking_angle, mis_angle_now, asy_now, energy_keV_now,
#                               range_kev_now, mnow, mill_ind, cent_now, paral_now)

plt.plot(vector_freq, vector_cry_freq.real)
#plt.plot(vector_rocking_angle, vector_cry_poly.imag)
plt.show()
'''
print arcsec(width_angle_o(energy_keV_now, asy_now, mnow, mill_ind))
print arcsec(width_angle_o(energy_keV_now, - (np.pi / 2), mnow, mill_ind))


# Plots intensity and phase vs rocking angle

vector_rocking_angle = np.arange(minrange,maxrange,step)
vector_cry_poly = wave_crystal_poly(vector_rocking_angle, mis_angle_now, asy_now, energy_keV_now,
                               range_kev_now, mnow, mill_ind, cent_now, paral_now)
#vector_cry_poly = [vector_rocking_angle, intensity(wave_crystal_poly(vector_rocking_angle, mis_angle_now, asy_now, energy_keV_now,
  #                             range_kev_now, mnow, mill_ind, cent_now, paral_now))]

print [np.array((vector_rocking_angle,intensity(vector_cry_poly))).transpose()]
plt.plot(vector_rocking_angle, intensity(vector_cry_poly))
plt.show()