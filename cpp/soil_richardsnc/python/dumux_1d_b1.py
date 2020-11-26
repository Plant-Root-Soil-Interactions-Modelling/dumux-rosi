# Benchmark 1 (in 1D)
#
# compares the dumux solution to the analytical solution (Figure 2abc Vanderborght et al 2005)
#
# Just a quick check that water movement including a tracer is plausible.
#
# D. Leitner, 2020
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
os.chdir("../../../build-cmake/cpp/soil_richardsnc")

# run dumux
os.system("./richardsnc1d input/b1a_1d.input")
os.system("./richardsnc1d input/b1b_1d.input")
os.system("./richardsnc1d input/b1c_1d.input")

# Figure 2a
p_, z1_ = read3D_vtp_data("benchmark1d_1a-00001.vtp", 2)  # matric potential
h1_ = vg.pa2head(p_)
ax1.plot(h1_, z1_[:,0] * 100, "r+")
ax1b = ax1.twiny()
c_, z_ = read3D_vtp_data("benchmark1d_1a-00000.vtp", 13)  # Solute concentration
ax1b.plot(c_, z_[:,0] * 100, "g:")
print(max(c_))
c_, z_ = read3D_vtp_data("benchmark1d_1a-00001.vtp", 13)
ax1b.plot(c_, z_[:,0] * 100, "r:")

# Figure 2b
p_, z2_ = read3D_vtp_data("benchmark1d_1b-00001.vtp",2)
h2_ = vg.pa2head(p_)
ax2.plot(h2_, z2_[:,0] * 100, "r+")
ax2b = ax2.twiny()
c_, z_ = read3D_vtp_data("benchmark1d_1b-00000.vtp", 13)  # Solute concentration
ax2b.plot(c_, z_[:,0] * 100, "g:")
c_, z_ = read3D_vtp_data("benchmark1d_1b-00001.vtp", 13)
ax2b.plot(c_, z_[:,0] * 100, "r:")

# Figure 2c
p_, z3_ = read3D_vtp_data("benchmark1d_1c-00001.vtp",2)
h3_ = vg.pa2head(p_)
ax3.plot(h3_, z3_[:,0] * 100, "r+")
ax3b = ax3.twiny()
c_, z_ = read3D_vtp_data("benchmark1d_1c-00000.vtp", 13)  # Solute concentration
ax3b.plot(c_, z_[:,0] * 100, "g:")
c_, z_ = read3D_vtp_data("benchmark1d_1c-00001.vtp", 13)
ax3b.plot(c_, z_[:,0] * 100, "r:")

plt.show()
