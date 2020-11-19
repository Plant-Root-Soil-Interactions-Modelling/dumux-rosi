"""
TODO
"""
import sys; sys.path.append("../../../../CPlantBox/src/python_modules/")

import os
import matplotlib.pyplot as plt
from vtk_tools import *
from math import *
import van_genuchten as vg

name = "singleroot"  # this name should be unique

# go to the right place
path = os.path.dirname(os.path.realpath(__file__))
os.chdir(path)
os.chdir("../../../build-cmake/cpp/roots_1p/")

# run dumux
os.system("./rootsystem_periodic input/" + name + ".input")

# move results to folder 'name'
if not os.path.exists("results_" + name):
    os.mkdir("results_" + name)
os.system("mv " + name + "* " + "results_" + name + "/")
os.system("cp input/" + name + ".input " + "results_" + name + "/")

# plot
p_, z_ = read3D_vtp_data("results_" + name + "/" + name + "-00001.vtp")
h_ = vg.pa2head(p_)
plt.plot(h_, z_[:, 2], "r+")  # cell data
plt.ylabel("Depth (m)")
plt.xlabel("Xylem pressure (cm)")
plt.savefig("results_" + name + ".pdf")
plt.show()
