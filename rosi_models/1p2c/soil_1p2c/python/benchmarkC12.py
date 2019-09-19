""" Benchmark C12, water uptake in roots """

import os
import matplotlib.pyplot as plt
from vtk_tools import *
import van_genuchten as vg

# go to the right place
path = os.path.dirname(os.path.realpath(__file__))
os.chdir(path)
os.chdir("../../../build-cmake/rosi_benchmarking/soil")

np_ = 4  # number of processors
# run dumux
if np_ == 1:
    os.system("./richards3d benchmarks_3d/benchmarkC12.input")
else:
    os.system("mpirun -n " + str(np_) + " ./richards3d benchmarks_3d/benchmarkC12.input -Grid.Overlap 0")

# Figure
s_, p_, y_ = read3D_vtp("benchmarkC12-00001", np_, 1)
h1_ = vg.pa2head(p_)
plt.plot(y_ * 100, h1_, "r+")
plt.ylabel('$\psi$ (cm)')
plt.xlabel('y axis (cm)')
plt.show()
