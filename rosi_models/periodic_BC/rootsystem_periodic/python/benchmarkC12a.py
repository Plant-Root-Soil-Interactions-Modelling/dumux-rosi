#
# Calcultates root system benchmark C1.2a (water uptake by a root system, constant conductivities)
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
os.chdir("../../../build-cmake/rosi_benchmarking/rootsystem")

# delete old result
os.system("rm benchmarkC12a-00001.vtp")

# run dumux
os.system("./rootsystem input/benchmarkC12a.input")

# plot
p_, z_ = read3D_vtp_data("benchmarkC12a-00001.vtp", False)
h_ = vg.pa2head(p_)
plt.plot(h_, z_[:, 2], "r+")  # cell data
plt.ylabel("Depth (m)")
plt.xlabel("Xylem pressure (cm)")
np.savetxt("dumux_b2", np.vstack((z_[:, 2], h_)), delimiter = ',')
plt.show()

