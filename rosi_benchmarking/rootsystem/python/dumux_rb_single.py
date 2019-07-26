#
# compares the dumux solution of 1d geometry created by rootbox to its analytical solution
#
# D. Leitner, 2019
#

import os
import matplotlib.pyplot as plt
from vtk_tools import *
from math import *
import van_genuchten as vg

import dumux_b1

# go to the right place
path = os.path.dirname(os.path.realpath(__file__))
os.chdir(path)
os.chdir("../../../build-cmake/rosi_benchmarking/rootsystem")
os.system("rm rb_single-00001.vtp")  # delete old results

# run dumux
os.system("./rootsystem_rb input/rb_single.input")

# plot
p_, z_ = read3D_vtp_data("rb_single-00001.vtp", False)
h_ = vg.pa2head(p_)

# plt.plot(h_, z_[:, 2] + 0.03, "r+")  # cell data
# plt.show()

dumux_b1.ax1.plot(h_, z_[:, 2] + 0.03, "r+")  # cell data
dumux_b1.plt.show()
