# Benchmark 1 (in 1D)
#
# compares the dumux solution to the analytical solution (Figure 2abc Vanderborght et al 2005)
#
# D. Leitner, 2018
#
import sys; sys.path.append("../../../../CPlantBox/src/python_modules/")
import sys; sys.path.append("../../../python/soil/")  # for the analytical solutions

import os
import matplotlib.pyplot as plt
from analytic_b1 import *
from vtk_tools import *
import van_genuchten as vg

# fig, (ax1, ax2, ax3) = plt.subplots(1, 3)

# go to the right place
path = os.path.dirname(os.path.realpath(__file__))
os.chdir(path)
os.chdir("../../../build-cmake/cpp/soil_richards")

# run dumux
os.system("./richards1d input/b1a_1d.input")
os.system("./richards1d input/b1b_1d.input")
os.system("./richards1d input/b1c_1d.input")

# Figure 2a
p_, z1_ = read3D_data("benchmark1d_1a-00001.vtp",1, 2, None)  # pressure [PA] is at index 2
h1_ = vg.pa2head(p_)
ax1.plot(h1_, z1_[:, 0] * 100, "r+")  # z coordinate is at position 0 in 1D models

# Figure 2b
p_, z2_ = read3D_data("benchmark1d_1b-00001.vtp", 1, 2, None)
h2_ = vg.pa2head(p_)
ax2.plot(h2_, z2_[:, 0] * 100, "r+")

# Figure 2c
p_, z3_ = read3D_data("benchmark1d_1c-00001.vtp",1, 2, None)
h3_ = vg.pa2head(p_)
ax3.plot(h3_, z3_[:, 0] * 100, "r+")

# print(z1_.shape, z2_.shape, z3_.shape, h1_.shape, h2_.shape, h3_.shape)
np.savetxt("dumux1d_b1", np.vstack((z1_[:, 0], h1_, z2_[:, 0], h2_, z3_[:, 0], h3_)), delimiter = ",")

plt.show()
