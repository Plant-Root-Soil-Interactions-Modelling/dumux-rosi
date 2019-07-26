""" Wet Top Scenario Felicien, 1 cm^3 resolution """

import os
import matplotlib.pyplot as plt
from vtk_tools import *
import van_genuchten as vg

# go to the right place
path = os.path.dirname(os.path.realpath(__file__))
os.chdir(path)
os.chdir("../../../build-cmake/rosi_benchmarking/soil")

np_ = 8  # number of processors

# run dumux
if np_ == 1:
    os.system("./richards3d benchmarks_3d/swtop.input")
else:
    os.system("mpirun -n " + str(np_) + " ./richards3d benchmarks_3d/swtop.input -Grid.Overlap 0")

# Figure
s_, p_, z1_ = read3D_vtp("swtop-00001", np_)
h1_ = vg.pa2head(p_)
plt.plot(h1_, (z1_ - 1.26) * 100, "r+")
plt.xlabel('$\psi$ (cm)')
plt.ylabel('Depth (cm)')

plt.show()

