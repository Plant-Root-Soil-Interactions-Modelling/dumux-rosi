# Benchmark 3 (in 3D)
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
np_ = 1  # number of processors
if np_ == 1:
    os.system("./richards3d input/b3a_3d.input")
#     os.system("./richards3d input/b3b_3d.input")
#     os.system("./richards3d input/b3c_3d.input")
else:
    pass
#    os.system("mpirun -n " + str(np_) + " ./richards3d input/b3a_3d.input -Grid.Overlap 1")
#     os.system("mpirun -n " + str(np_) + " ./richards3d input/b3b_3d.input -Grid.Overlap 1")
#     os.system("mpirun -n " + str(np_) + " ./richards3d input/b3c_3d.input -Grid.Overlap 1")

ex = []  # list for data for export

# result (Figure 4a)
for i in range(0, 3):
    p_, z_ = read3D_data("benchmark3d_3a-0000" + str(i + 1), np_, 2)
    h_ = vg.pa2head(p_)
    theta_ = vg.water_content(h_, sand)
    ax1.plot(theta_, z_[:, 2] * 100, "r+")
    ex.append(z_[:, 2])
    ex.append(theta_)

# result (Figure 4b)
for i in range(0, 3):
    p_, z_ = read3D_data("benchmark3d_3b-0000" + str(i + 1), np_, 2)
    h_ = vg.pa2head(p_)
    theta_ = vg.water_content(h_, loam)
    ax2.plot(theta_, z_[:, 2] * 100, "r+")
    ex.append(z_[:, 2])
    ex.append(theta_)

# result (Figure 4c)
for i in range(0, 3):
    p_, z_ = read3D_data("benchmark3d_3c-0000" + str(i + 1), np_, 2)
    h_ = vg.pa2head(p_)
    theta_ = vg.water_content(h_, clay)
    ax3.plot(theta_, z_[:, 2] * 100, "r+")
    ex.append(z_[:, 2])
    ex.append(theta_)

np.savetxt("dumux3d_b3", np.vstack(ex), delimiter=",")

plt.show()

