# Benchmark 1 (3D unstructured grid)
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
np_ = 4  # number of processors
if np_ == 1:
    pass
    os.system("./richardsUG benchmarks_ug/b1a.input")
    os.system("./richardsUG benchmarks_ug/b1b.input")
    os.system("./richardsUG benchmarks_ug/b1c.input")
else:
    pass
    os.system("mpirun -n " + str(np_) + " ./richardsUG benchmarks_ug/b1a.input -Grid.Overlap 0")
    os.system("mpirun -n " + str(np_) + " ./richardsUG benchmarks_ug/b1b.input -Grid.Overlap 0")
    os.system("mpirun -n " + str(np_) + " ./richardsUG benchmarks_ug/b1c.input -Grid.Overlap 0")

# Figure 2a
s_, p_, z1_ = read3D_vtp("benchmarkUG_1a-00001", np_)
h1_ = vg.pa2head(p_)
ax1.plot(h1_, (z1_ - 2) * 100, "r+")

# Figure 2b
s_, p_, z2_ = read3D_vtp("benchmarkUG_1b-00001", np_)
h2_ = vg.pa2head(p_)
ax2.plot(h2_, (z2_ - 2) * 100, "r+")

# Figure 2c
s_, p_, z3_ = read3D_vtp("benchmarkUG_1c-00001", np_)
h3_ = vg.pa2head(p_)
ax3.plot(h3_, (z3_ - 2) * 100, "r+")

np.savetxt("dumuxUG_b1", np.vstack((z1_ - 2, h1_, z2_ - 2, h2_, z3_ - 2, h3_)), delimiter = ",")

plt.show()

