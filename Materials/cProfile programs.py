import cProfile
import math
import cry2_all_materials_vector as cry
mnow = 'Si'
cProfile.run('cry.wave_crystal_function(4, -math.pi/2 + 4 * math.pi / 180, 25, mnow, [2,2,0])')
