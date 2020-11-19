""" Benchmark C12, soil part """
import sys; sys.path.append("../../../../CPlantBox/src/python_modules/")

import os
import matplotlib.pyplot as plt
from vtk_tools import *
import van_genuchten as vg

# go to the right place
path = os.path.dirname(os.path.realpath(__file__))
os.chdir(path)
os.chdir("../../../build-cmake/cpp/soil_richards")

np_ = 1  # number of processors
# run dumux
if np_ == 1:
    os.system("./richards_ug input/benchmarkC12ug_3d.input")
else:
    os.system("mpirun -n " + str(np_) + " ./richardsUG input/benchmarkC12ug_3d.input -Grid.Overlap 0")

# Figure
p_, z1_ = read3D_data("benchmarkC12ug-00001", np_)
h1_ = vg.pa2head(p_)
plt.plot(h1_, z1_ * 100, "r+")
plt.xlabel('$\psi$ (cm)')
plt.ylabel('z axis (cm)')
plt.show()
