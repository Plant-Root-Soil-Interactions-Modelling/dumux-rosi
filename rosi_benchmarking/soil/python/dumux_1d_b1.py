# Benchmark 1 (in 1D)
#
# compares the dumux solution to the analytical solution (Figure 2abc Vanderborght et al 2005)
#
# D. Leitner, 2018
#

import os
import matplotlib.pyplot as plt
from analytic_b1 import *
from vtk_tools import *
import van_genuchten as vg

# go to the right place
path = os.path.dirname(os.path.realpath(__file__))
os.chdir(path)
os.chdir("../../../build-cmake/rosi_benchmarking/soil")

# run dumux
os.system("./richards1d benchmarks_1d/b1a.input")
os.system("./richards1d benchmarks_1d/b1b.input")
os.system("./richards1d benchmarks_1d/b1c.input")

# Figure 2a
s_, p_ = read1D_vtp_data("benchmark1d_1a-00001.vtp", False)
z_ = np.linspace(0, -200, len(s_))
h_ = vg.pa2head(p_)
ax1.plot(h_, z_, "r+")

# Figure 2b
s_, p_ = read1D_vtp_data("benchmark1d_1b-00001.vtp", False)
z_ = np.linspace(0, -200, len(s_))
h_ = vg.pa2head(p_)
ax2.plot(h_, z_, "r+")

# Figure 2c
s_, p_ = read1D_vtp_data("benchmark1d_1c-00001.vtp", False)
z_ = np.linspace(0, -200, len(s_))
h_ = vg.pa2head(p_)
ax3.plot(h_, z_, "r+")

plt.show()

