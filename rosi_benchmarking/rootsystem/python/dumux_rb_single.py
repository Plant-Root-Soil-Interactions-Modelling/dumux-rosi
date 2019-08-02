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
os.system("rm rb_singleT-00001.vtp")  # delete old results

# parameter
trans = 0.0017;  # 0.0017  # kg /day
print(trans)

# run dumux
os.system("./rootsystem_rb input/rb_single.input")
p_, z_ = read3D_vtp_data("rb_single-00001.vtp", False)
h_ = vg.pa2head(p_)

os.system("./rootsystem_rb input/rb_single.input -RootSystem.Collar.Transpiration {} -Problem.Name rb_singleT".format(trans))
p2_, z2_ = read3D_vtp_data("rb_singleT-00001.vtp", False)  #
h2_ = vg.pa2head(p2_)

with open("rb_singleT_actual_transpiration.txt", 'r') as f:
    d = np.loadtxt(f, delimiter = ',')

dumux_b1.ax1.plot(h_, z_[:, 2] + 0.03, "r+")  # cell data
dumux_b1.ax3.plot(h2_, z2_[:, 2] + 0.03, "r+")  # cell data
dumux_b1.plt.show()

