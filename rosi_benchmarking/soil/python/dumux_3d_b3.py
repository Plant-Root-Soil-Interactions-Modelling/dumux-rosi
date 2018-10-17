# Benchmark 3 (in 3D)
#
# compares the dumux solution to the analytical solution (Figure 4abc Vanderborght et al 2005)
#
# D. Leitner, 2018
#

import os
import matplotlib.pyplot as plt
from analytic_b3 import *
from vtk_tools import *
import van_genuchten as vg

# go to the right place
path = os.path.dirname(os.path.realpath(__file__))
os.chdir(path)
os.chdir("../../../build-cmake/rosi_benchmarking/soil")

# run dumux
np_ = 16  # number of processors
if np_ == 1:
    os.system("./richards3d benchmarks_3d/b3a.input")
    os.system("./richards3d benchmarks_3d/b3b.input")
    os.system("./richards3d benchmarks_3d/b3c.input")
else:
    os.system("mpirun -n " + str(np_) + " ./richards3d benchmarks_3d/b3a.input -Grid.Overlap 0")
    os.system("mpirun -n " + str(np_) + " ./richards3d benchmarks_3d/b3b.input -Grid.Overlap 0")
    os.system("mpirun -n " + str(np_) + " ./richards3d benchmarks_3d/b3c.input -Grid.Overlap 0")

ex = []  # list for data for export

# result (Figure 4a)
for i in range(0, 3):
    s_, p_, z_ = read3D_vtp("benchmark3d_3a-0000" + str(i + 1) + ".vtp", np_)
    h_ = vg.pa2head(p_)
    theta_ = vg.water_content(h_, sand)
    ax1.plot(theta_, (z_ - 2) * 100, "r+")
    ex.append(z_ - 2)
    ex.append(theta_)

# result (Figure 4b)
for i in range(0, 3):
    s_, p_, z_ = read3D_vtp("benchmark3d_3b-0000" + str(i + 1) + ".vtp", np_)
    h_ = vg.pa2head(p_)
    theta_ = vg.water_content(h_, loam)
    ax2.plot(theta_, (z_ - 2) * 100, "r+")
    ex.append(z_ - 2)
    ex.append(theta_)

# result (Figure 4c)
for i in range(0, 3):
    s_, p_, z_ = read3D_vtp("benchmark3d_3c-0000" + str(i + 1) + ".vtp", np_)
    h_ = vg.pa2head(p_)
    theta_ = vg.water_content(h_, clay)
    ax3.plot(theta_, (z_ - 2) * 100, "r+")
    ex.append(z_ - 2)
    ex.append(theta_)

np.savetxt("dumux3d_b3", np.vstack(ex), delimiter = ",")

plt.show()

