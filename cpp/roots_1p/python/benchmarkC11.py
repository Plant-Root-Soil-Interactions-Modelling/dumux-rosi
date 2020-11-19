""" Benchmark C11, single root in thin soil layer, root part (TODO?)"""
import sys; sys.path.append("../../../../CPlantBox/src/python_modules/")

import os
import matplotlib.pyplot as plt
from vtk_tools import *
import van_genuchten as vg

# go to the right place
path = os.path.dirname(os.path.realpath(__file__))
os.chdir(path)
os.chdir("../../../build-cmake/cpp/roots_1p")

# run dumux
os.system("./rootsystem input/benchmarkC11.input")

# plot
p_, z_ = read3D_vtp_data("benchmarkC11-00001.vtp")
h_ = vg.pa2head(p_)
print(h_.shape, p_.shape)  # ???

plt.plot(h_, z_[:, 2], "r+")  # cell data
plt.ylabel("Depth (m)")
plt.xlabel("Xylem pressure (cm)")
plt.show()
