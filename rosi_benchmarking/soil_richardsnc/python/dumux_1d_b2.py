# Benchmark 2 (in 1D)
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
# os.system("./richardsnc1d input/b2_1d.input")
os.system("./richardsnc1d input/b2_1d.input -Soil.Grid.Cells 1000")  # high res looks nice

ax2 = ax1.twiny()

# Soil matric potential
s_, p_, z_ = read1D_vtp_data("benchmark1d_2-00000.vtp")
h_ = vg.pa2head(p_)
ax1.plot(h_, z_ * 100, "g+")
s_, p_, z_ = read1D_vtp_data("benchmark1d_2-00001.vtp")
h_ = vg.pa2head(p_)
ax1.plot(h_, z_ * 100, "r+")

ax1.legend(["analytic", "initial", "numeric"])

# Solute concentration
s_, c_, z_ = read1D_vtp_data("benchmark1d_2-00000.vtp", 13)
ax2.plot(c_, z_ * 100, "g:")
s_, c_, z_ = read1D_vtp_data("benchmark1d_2-00001.vtp", 13)
ax2.plot(c_, z_ * 100, "r:")

ax2.legend(["initial", "final"])

# np.savetxt("dumux1d_b2", np.vstack((z_, h_)), delimiter = ",")

plt.show()
