import math  # Mathematic functions for non-complex numbers.
import cmath # Mathematic functions for complex numbers.
import numpy as np # NumPy - the fundamental package for scientific computing with Python.
import matplotlib.pyplot as plt #  Collection of command style functions that make matplotlib work like MATLAB.
# import scipy.constants as sci

re = 2.81794*10**-15
hplanck = 4.1356692*10**-15
c = 2.99792458*10**8
Na = 6.02*10**23
#Changed: \[Epsilon)o = 8.854*10^-12 --> Epo
Epo = 8.854*10**-12
e = 1.602*10**-19
me = 9.109*10**-31

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


''' string_material '''
def density(string_material):
    if string_material == 'Si':
        den = costantim[6][5]  # Silicon is the 6th element in the table
    return den
'''def electron_density(material): return  costantim[material][5]*Na* \
        (costantim[[material, 3]]/costantim[[material, 4]])*10**6;'''
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
def b1(material): return 2*math.pi*np.cross(a2(material),a3(material))/volume(material)
def b2(material): return 2*math.pi*np.cross(a3(material),a1(material))/volume(material)
def b3(material): return 2*math.pi*np.cross(a1(material),a2(material))/volume(material)
def d_space(list_material, list_miller):
    d_1 = list_miller[0]*b1(list_material)
    d_2 = list_miller[0]*b2(list_material)
    d_3 = list_miller[0]*b3(list_material)
    return 2 * math.pi / np.linalg.norm(d_1 + d_2 + d_3)



form_factor_real_Si220 = np.loadtxt("C:\Python27\Materials\\form_factor_real_Si220.txt")
form_factor_real_Si111 = np.loadtxt("C:\Python27\Materials\\form_factor_real_Si111.txt")
form_factor_imag_Si = np.loadtxt("C:\Python27\Materials\\f2_Si_50keV.txt")




def form_factor_real(energy_keV, string_material, list_miller):
    ''' Real part of the form factor. Insert the name of the material in string_material
        and the name of the Miller index in list_miller
    '''
    if (string_material == 'Si') and (list_miller == [2,2,0]):
        form_real = form_factor_real_Si220
    elif (string_material == 'Si') and (list_miller == [1,1,1]):
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
    if (string_material == 'Si') and (list_miller == [2,2,0]):
        str_factor = 8
    elif (string_material == 'Si') and (list_miller == [-2,-2,-0]):
        str_factor = 8
    elif (string_material == 'Si') and (list_miller == [1,1,1]):
        str_factor = 4 * (1 - 1j)
    elif (string_material == 'Si') and (list_miller == [-1,-1,-1]):
        str_factor = 4 * (1 + 1j)
    return str_factor
    

def Fc(energy_keV, string_material, list_miller):
    ''' Structure factor (complex) '''
    form_fac = form_factor(energy_keV, string_material, list_miller)
    geom_fac = structure_factor_geometric(string_material, list_miller)
    return form_fac * geom_fac
    
def Chih(energy_keV, string_material, list_miller):
    str_factor = Fc(energy_keV, string_material, list_miller)
    numerator = str_factor * electron_radius * (wavelength(energy_keV)) ** 2
    denominator = math.pi * volume(material)
    return - numerator/denominator

def Chihplus(energy_keV, string_material, list_miller):
    return Chih(energy_keV, string_material, list_miller)

def Chihminus(energy_keV, string_material, list_miller):
    return Chih(energy_keV, string_material, - list_miller)

def beta_crystal(energy_keV, string_material, list_miller):
    Chihpl = Chihplus(energy_keV, string_material, list_miller)
    Chihmi = Chihminus(energy_keV, string_material, list_miller)
    cmath.phase(cmath.sqrt(Chihpl * Chihmi))

print(
d_space('Si', [2,2,0])
)


'''
print(fimag[6])

print(fimag[6][:,0])
#plt.plot(fimag[6])

plt.plot(fimag[6][:,0],fimag[6][:,1])
axes = plt.gca()
axes.set_xlim([20,30])
axes.set_ylim([0,0.1])

plt.show()
'''
'''
def ffactor(energy_keV, h, k, l, material):
    fsum1 = ffa1 * math.exp(-ffb1 / (4 * (d_space(h,k,l,material) * 10 ** 10) ** 2))
    fsum2 = ffa2 * math.exp(-ffb2 / (4 * (d_space(h,k,l,material) * 10 ** 10) ** 2))
    fsum3 = ffa3 * math.exp(-ffb3 / (4 * (d_space(h,k,l,material) * 10 ** 10) ** 2))
    fsum4 = ffa4 * math.exp(-ffb1 / (4 * (d_space(h,k,l,material) * 10 ** 10) ** 2))
    return ffc(material) + (fsum1 + fsum2 + fsum3 + fsum4) + 1j * fdue(energy_keV, material)
'''

