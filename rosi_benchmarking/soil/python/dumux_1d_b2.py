# Benchmark 2 (in 1D)
#
# compares the dumux solution to the analytical solution (Figure 3, Vanderborght et al 2005)
#
# D. Leitner, 2018
#
import os
import matplotlib.pyplot as plt
from vtk_tools import *
import van_genuchten as vg
import analytic_b2

# go to the right place
path = os.path.dirname(os.path.realpath(__file__))
os.chdir(path)
os.chdir("../../../build-cmake/rosi_benchmarking/soil")

# run dumux
os.system("./richards1d benchmarks_1d/b2.input")

# result dumux jan1 (Figure 2a)
s_, p_ = read1D_vtp_data("benchmark1d_2-00001.vtp", False)
z_ = np.linspace(0, -54, len(s_))
h_ = vg.pa2head(p_)
plt.plot(h_, z_, "r+")

plt.show()
