"""
Benchmark M32b

 Calcultates root system benchmark M32b (static root system, age dependent conductivities)

 D. Leitner, 2019

"""

import os
import matplotlib.pyplot as plt
from vtk_tools import *
from math import *
import van_genuchten as vg

# go to the right place
path = os.path.dirname(os.path.realpath(__file__))
os.chdir(path)
os.chdir("../../../build-cmake/rosi_benchmarking/roots_1p")

# delete old result
#os.system("rm xylem-depth-00001.vtp")

# run dumux
#os.system("./rootsystem input/xylem-depth.input")

# plot
p_, z_ = read3D_vtp_data("mri_default_coupled_velocity2-00101.vtp")
h_ = vg.pa2head(p_)
plt.plot(h_, z_[:, 2], "r+")  # cell data
xmin = min(h_)
xmax = max(h_)
print("from ", xmin, "to", xmax, " cm pressure head")

plt.ylabel("Depth (m)")
plt.xlabel("Xylem pressure (cm)")
np.savetxt("xylem-depth", np.vstack((100 * z_[:, 2], h_)), delimiter = ',')
plt.show()
