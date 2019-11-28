#
# Root box root sytem simulation (see rb_rootsytem.input)
#
# D. Leitner, 2018
#

import os
import matplotlib.pyplot as plt
from vtk_tools import *
from math import *
import van_genuchten as vg

# go to the right place
path = os.path.dirname(os.path.realpath(__file__))
os.chdir(path)
os.chdir("../../../build-cmake/rosi_benchmarking/roots_1p")

# run dumux
os.system("./rootsystem_rb input/rb_rootsystem_trans.input")

# plot
p_, z_ = read3D_vtp_data("rb_rootsystem_trans-00001.vtp", False)
h_ = vg.pa2head(p_)
plt.plot(h_, z_[:, 2], "r+")  # cell data
plt.ylabel("Depth (m)")
plt.xlabel("Xylem pressure (cm)")
plt.show()
