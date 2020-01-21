# Benchmark 2 (in 3D)
#
# compares the dumux solution to the analytical solution (Figure 3, Vanderborght et al 2005)
#
# For RichardsNC checks the advective and diffusive flow for plausibility
#
# D. Leitner, 2020
#
import os
import matplotlib.pyplot as plt
from vtk_tools import *
import van_genuchten as vg

fig, ax1 = plt.subplots()

import analytic_b2

# go to the right place
path = os.path.dirname(os.path.realpath(__file__))
os.chdir(path)
os.chdir("../../../build-cmake/rosi_benchmarking/soil_richardsnc")

# run dumux
os.system("./richardsnc3d input/b2_3d.input")

ax2 = ax1.twiny()

# Soil matric potential
s_, p_, z_ = read3D_vtp("benchmark3d_2-00000")
h_ = vg.pa2head(p_)
ax1.plot(h_, z_ * 100, "g+")
s_, p_, z_ = read3D_vtp("benchmark3d_2-00001")
h_ = vg.pa2head(p_)
ax1.plot(h_, z_ * 100, "r+")

ax1.legend(["analytic", "initial", "numeric"])

# Solute concentration
s_, c_, z_ = read3D_vtp("benchmark3d_2-00000", 1, 13)
ax2.plot(c_, z_ * 100, "g:")
s_, c_, z_ = read3D_vtp("benchmark3d_2-00001", 1, 13)
ax2.plot(c_, z_ * 100, "r:")

ax2.legend(["initial", "final"])

# np.savetxt("dumux1d_b2", np.vstack((z_, h_)), delimiter = ",")

plt.show()
