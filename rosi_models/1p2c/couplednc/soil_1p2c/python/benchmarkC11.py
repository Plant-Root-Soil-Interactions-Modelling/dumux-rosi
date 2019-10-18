""" Benchmark C11, single root in thin soil layer, soil part """

import os
import matplotlib.pyplot as plt
from vtk_tools import *
import van_genuchten as vg

# go to the right place
path = os.path.dirname(os.path.realpath(__file__))
os.chdir(path)
os.chdir("../../../build-cmake/rosi_benchmarking/soil")

soiltype = 1  # sand, loam, clay

np_ = 1  # number of processors
# run dumux
if np_ == 1:
    os.system("./richards3d benchmarks_3d/benchmarkC11.input -Soil.Layer.Number {}".format(soiltype))
else:
    os.system("mpirun -n " + str(np_) + " ./richards3d benchmarks_3d/benchmarkC11.input -Soil.Layer.Number {} -Grid.Overlap 0".format(soiltype))

# Figure
s_, p_, y_ = read3D_vtp("benchmarkC11-00001", np_, 1)
h1_ = vg.pa2head(p_)
plt.plot(y_ * 100, h1_, "r+")
plt.ylabel('$\psi$ (cm)')
plt.xlabel('y axis (cm)')
plt.show()
