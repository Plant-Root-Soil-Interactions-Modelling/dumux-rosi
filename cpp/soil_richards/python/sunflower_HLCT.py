""" sunflower_HLCT, soil part NOT WORKING? """
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
    os.system("./richards3d input/sunflower_HLCT.input")
else:
    os.system("mpirun -n " + str(np_) + " ./richards3d input/sunflower_HLCT.input -Grid.Overlap 0")

# Figure
p_, z1_ = read3D_data("sunflower_HLCT-00001", np_, 2)
h1_ = vg.pa2head(p_)
plt.plot(h1_, z1_[:,2] * 100, "r+")
plt.xlabel('$\psi$ (cm)')
plt.ylabel('z axis (cm)')
plt.show()
