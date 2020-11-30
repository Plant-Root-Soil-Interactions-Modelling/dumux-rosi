# Run small coupled test senarios
import sys; sys.path.append("../../../../CPlantBox/src/python_modules/")

import os
import matplotlib.pyplot as plt
from vtk_tools import *
import van_genuchten as vg

# go to the right place
path = os.path.dirname(os.path.realpath(__file__))
os.chdir(path)
os.chdir("../../../build-cmake/cpp/coupled_1p_richards")

os.system("./coupled input/benchmark_phosphate.input")

# Figure
p_, z1_ = read3D_data("benchmark_phosphate2-00011.vtp")
h1_ = vg.pa2head(p_)
plt.plot(h1_, z1_[:,0] * 100, "r+")
plt.xlabel('$\psi$ (cm)')
plt.ylabel('Depth (cm)')

plt.show()
