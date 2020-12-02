""" Phospahte benchmark for comparison with Comsol, single root in thin soil layer, root part """
import sys; sys.path.append("../../../../CPlantBox/src/python_modules/")

import os
import matplotlib.pyplot as plt
from vtk_tools import *
import van_genuchten as vg

# go to the right place
path = os.path.dirname(os.path.realpath(__file__))
os.chdir(path)
os.chdir("../../../build-cmake/cpp/roots_1pnc")

# run dumux
os.system("./rootsystem_1p2c input/benchmark_phosphate.input")

# plot
p_, z_ = read3D_vtp_data("benchmark_phosphate-00001.vtp")
h_ = vg.pa2head(p_)
plt.plot(h_, z_[:, 2], "r+")  # cell data
plt.ylabel("Depth (m)")
plt.xlabel("Xylem pressure (cm)")
plt.show()
