""" sunflower_HLCT, soil part """
import sys; sys.path.append("../../../python/modules/")

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
    os.system("./richards3d input/sunflower_HLCT.input")
else:
    os.system("mpirun -n " + str(np_) + " ./richards3d input/sunflower_HLCT.input -Grid.Overlap 0")

# Figure
s_, p_, z1_ = read3D_vtp("sunflower_HLCT-00001", np_)
h1_ = vg.pa2head(p_)
plt.plot(h1_, z1_ * 100, "r+")
plt.xlabel('$\psi$ (cm)')
plt.ylabel('z axis (cm)')
plt.show()
