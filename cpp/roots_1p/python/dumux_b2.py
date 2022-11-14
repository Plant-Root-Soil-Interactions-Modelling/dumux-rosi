"""
Benchmark M32a

 Calcultates root system benchmark M32a (static root system, constant conductivities)

 D. Leitner, 2019

"""
import sys; sys.path.append("../../../../CPlantBox/src/python_modules/")

import os
import matplotlib.pyplot as plt
from vtk_tools import *
import van_genuchten as vg

# go to the right place
path = os.path.dirname(os.path.realpath(__file__))
os.chdir(path)
os.chdir("../../../build-cmake/cpp/roots_1p")

# delete old result
os.system("rm benchmark2-00001.vtp")

# run dumux
os.system("./rootsystem input/b2.input")

# plot
p_, z_ = read3D_vtp_data("benchmark2-00001.vtp")
h_ = vg.pa2head(p_)
plt.plot(h_, z_[:, 2], "r+")  # cell data
plt.ylabel("Depth (m)")
plt.xlabel("Xylem pressure (cm)")
np.savetxt("dumux_m32a", np.vstack((100 * z_[:, 2], h_)), delimiter = ',')
# plt.show()

