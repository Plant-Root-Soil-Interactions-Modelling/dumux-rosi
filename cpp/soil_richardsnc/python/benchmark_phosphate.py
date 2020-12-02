""" 
Phophate benchmark (to compare to Comsol solution), single root in thin soil layer, soil part 
"""

import sys; sys.path.append("../../../../CPlantBox/src/python_modules/")

import os
import matplotlib.pyplot as plt
from vtk_tools import *
import van_genuchten as vg

# go to the right place
path = os.path.dirname(os.path.realpath(__file__))
os.chdir(path)
os.chdir("../../../build-cmake/cpp/soil_richardsnc")

soiltype = 1  # sand, loam, clay

np_ = 1  # number of processors
# run dumux
if np_ == 1:
    os.system("./richardsnc3d input/benchmark_phosphate.input -Soil.Layer.Number {}".format(soiltype))
else:
    os.system("mpirun -n " + str(np_) + " ./richards3d input/benchmarkC11_3d.input -Soil.Layer.Number {} -Grid.Overlap 0".format(soiltype))

# Figure
p_, y_ = read3D_data("benchmark_phosphate-00001", np_, 2)
h1_ = vg.pa2head(p_)
plt.plot(h1_, y_[:,2] * 100, "r+")
plt.xlabel('$\psi$ (cm)')
plt.ylabel('z axis (cm)')
plt.show()
