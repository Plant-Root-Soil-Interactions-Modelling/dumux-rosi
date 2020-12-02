""" Benchmark C11, single root in thin soil layer, soil part """
import sys; sys.path.append("../../../../CPlantBox/src/python_modules/")

import os
import matplotlib.pyplot as plt
from vtk_tools import *
import van_genuchten as vg

# go to the right place
path = os.path.dirname(os.path.realpath(__file__))
os.chdir(path)
os.chdir("../../../build-cmake/cpp/soil_richards")

soiltype = 1  # sand, loam, clay

np_ = 1  # number of processors
# run dumux
if np_ == 1:
    os.system("./richards3d input/benchmarkC11_3d.input -Soil.Layer.Number {}".format(soiltype))
else:
    os.system("mpirun -n " + str(np_) + " ./richards3d input/benchmarkC11_3d.input -Soil.Layer.Number {} -Grid.Overlap 0".format(soiltype))

# Figure
p_, z_ = read3D_data("benchmarkC11-00001", np_, 1)
h1_ = vg.pa2head(p_)
plt.plot(z_[:, 1] * 100, h1_, "r+")
plt.ylabel('$\psi$ (cm)')
plt.xlabel('y axis (cm)')
plt.show()
