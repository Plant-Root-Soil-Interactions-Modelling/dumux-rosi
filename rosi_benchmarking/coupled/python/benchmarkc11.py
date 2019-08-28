# Run Benchmakrk C11

import os
import matplotlib.pyplot as plt
from vtk_tools import *
import van_genuchten as vg

# go to the right place
path = os.path.dirname(os.path.realpath(__file__))
os.chdir(path)
os.chdir("../../../build-cmake/rosi_benchmarking/coupled")

soiltype = 1  # sand, loam, clay
trans = 0 * 0.1 * 0.02 * 0.02 * 1.e-6 * 1000.  # cm/day -> kg/day
print(trans)
os.system("./coupled_seq input/benchmarkC11.input -Soil.Layer.Number {} -RootSystem.Collar.Transpiration {}"
          .format(soiltype, trans))

# Figure
s_, p_, y_ = read3D_vtp("benchmarkC11-00006", 1, 1)
h1_ = vg.pa2head(p_)
plt.plot(y_ * 100, h1_, "r+")
plt.ylabel('$\psi$ (cm)')
plt.xlabel('y axis (cm)')
plt.show()

# p2_, z2_ = read3D_vtp_data("benchmarkC11R-00001.vtp", False)
