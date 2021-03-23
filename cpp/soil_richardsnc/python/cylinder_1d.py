# radially symmetric cylinder (in 1D)
#
# compares the dumux solution to the comsol solution
#
# Just a quick check that water movement including a tracer is plausible.
#
# D. Leitner, 2020
#
""" TODO not working: Dune reported error: NumericalProblem check parameters """
import sys; sys.path.append("../../../../CPlantBox/src/python_modules/")

import os
import matplotlib.pyplot as plt
from vtk_tools import *
import van_genuchten as vg

# go to the right place
path = os.path.dirname(os.path.realpath(__file__))
os.chdir(path)
os.chdir("../../../build-cmake/cpp/soil_richardsnc")

# run dumux
os.system("./richardsnc1d_cyl input/cylinder_1d.input")
 
# plot dumux results
# c_, z_ = read3D_data("cylinder_1d-00001.vtp", data_index = 13)
# plt.plot((z_[:,0] - 0.0002) * 100, c_, "b:",)
# c_, z_ = read3D_data("cylinder_1d-00002.vtp", data_index = 13)
# plt.plot((z_[:,0] - 0.0002) * 100, c_, "b",)

# read and plot comsol data
os.chdir("../../../build-cmake/rosi_benchmarking/soil_richardsnc/python")
data = np.loadtxt("cylinder_1d_Comsol_P.txt", skiprows=8)
z_comsol = data[:, 0]
h_comsol = data[:, 25]
h_comsol2 = data[:, -1]
plt.plot(z_comsol, h_comsol, "r:", z_comsol, h_comsol2, "r")

plt.xlabel('distance from the root surface (cm)')
plt.ylabel('concentration (mol)')
plt.legend(["dumux, 10d", "dumux, 20d", "comsol, 10d", "comsol, 20d"], loc='lower right')

plt.show()

