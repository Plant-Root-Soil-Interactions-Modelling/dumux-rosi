# Benchmark 2 (in 1D)
#
# compares the dumux solution to the analytical solution (Figure 3, Vanderborght et al 2005)
#
# D. Leitner, 2018
#
import sys; sys.path.append("../../../../CPlantBox/src/python_modules/")
import sys; sys.path.append("../../../python/soil/")  # for the analytical solutions

import os
import matplotlib.pyplot as plt
from vtk_tools import *
import van_genuchten as vg
import analytic_b2

# go to the right place
path = os.path.dirname(os.path.realpath(__file__))
os.chdir(path)
os.chdir("../../../build-cmake/cpp/soil_richards")

# run dumux
os.system("./richards1d input/b2_1d.input")
# os.system("./richards1d input/b2_1d.input -Soil.Grid.Cells 1000")  # high res looks nice

# result dumux jan1 (Figure 2a)
p_, z_ = read3D_vtp_data("benchmark1d_2-00001.vtp", 2, None)  # pressure [PA] is at index 2
h_ = vg.pa2head(p_)
plt.plot(h_, z_[:, 0] * 100, "r+")  # z coordinate is at position 0 in 1D models

np.savetxt("dumux1d_b2", np.vstack((z_[:, 0], h_)), delimiter = ",")

plt.show()
