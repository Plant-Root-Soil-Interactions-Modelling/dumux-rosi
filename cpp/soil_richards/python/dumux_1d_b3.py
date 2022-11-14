# Benchmark 3 (in 1D)
#
# compares the dumux solution to the analytical solution (Figure 4abc Vanderborght et al 2005)
#
# D. Leitner, 2018
#
import sys; sys.path.append("../../../../CPlantBox/src/python_modules/")
import sys; sys.path.append("../../../python/soil/")  # for the analytical solutions

import os
import matplotlib.pyplot as plt
from analytic_b3 import *
from vtk_tools import *
import van_genuchten as vg

# go to the right place
path = os.path.dirname(os.path.realpath(__file__))
os.chdir(path)
os.chdir("../../../build-cmake/cpp/soil_richards")

# run dumux
os.system("./richards1d input/b3a_1d.input")
os.system("./richards1d input/b3b_1d.input")
os.system("./richards1d input/b3c_1d.input")

ex = []  # list for data for export

# result (Figure 4a)
for i in range(0, 3):
    p_, z_ = read3D_vtp_data("benchmark1d_3a-0000" + str(i + 1) + ".vtp", 2)  # pressure [PA] is at index 2
    h_ = vg.pa2head(p_)
    theta_ = vg.water_content(h_, sand)
    ax1.plot(theta_, z_[:, 0] * 100, "r+")
    ex.append(z_[:, 0])
    ex.append(theta_)

# result (Figure 4b)
for i in range(0, 3):
    p_, z_ = read3D_vtp_data("benchmark1d_3b-0000" + str(i + 1) + ".vtp", 2)
    h_ = vg.pa2head(p_)
    theta_ = vg.water_content(h_, loam)
    ax2.plot(theta_, z_[:, 0] * 100, "r+")
    ex.append(z_[:, 0])
    ex.append(theta_)

# result (Figure 4c)
for i in range(0, 3):
    p_, z_ = read3D_vtp_data("benchmark1d_3c-0000" + str(i + 1) + ".vtp", 2)
    h_ = vg.pa2head(p_)
    theta_ = vg.water_content(h_, clay)
    ax3.plot(theta_, z_[:, 0] * 100, "r+")
    ex.append(z_[:, 0])
    ex.append(theta_)

np.savetxt("dumux1d_b3", np.vstack(ex), delimiter = ",")

# plt.show()

