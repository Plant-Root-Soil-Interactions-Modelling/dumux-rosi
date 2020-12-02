# radially symmetric cylinder (in 1D)
#
# compares the dumux solution to the comsol solution
#
# Just a quick check that water movement including a tracer is plausible.
#
# D. Leitner, 2020
#
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
os.system("./richardsnc1d_cyl input/cylinder_1d_water.input")

# plot dumux results
p_, z_ = read3D_vtp_data("cylinder_1d_water-00001.vtp",2)
h_ = vg.pa2head(p_)
plt.plot((z_[:,0] - 0.0002) * 100, h_, "b",)

# read and plot comsol data
os.chdir("../../../build-cmake/cpp/soil_richardsnc/python")
data = np.loadtxt("cylinder_1d_Comsol_water.txt", skiprows=8)
z_comsol = data[:, 0]
h_comsol = data[:, 25]
plt.plot(z_comsol, h_comsol, "r")

plt.xlabel('distance from the root surface (cm)')
plt.ylabel('pressure head (cm)')
plt.legend(["dumux", "comsol"], loc='lower right')

plt.show()

