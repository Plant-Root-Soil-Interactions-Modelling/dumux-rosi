# Benchmark 1 (3D unstructured grid)
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
np_ = 1  # number of processors
if np_ == 1:
    os.system("./richards_ug input/b1a_ug.input")
    os.system("./richards_ug input/b1b_ug.input")
    os.system("./richards_ug input/b1c_ug.input")
else:
    os.system("mpirun -n " + str(np_) + " ./richards_ug input/b1a_ug.input -Grid.Overlap 1")
    os.system("mpirun -n " + str(np_) + " ./richards_ug input/b1b_ug.input -Grid.Overlap 1")
    os.system("mpirun -n " + str(np_) + " ./richards_ug input/b1c_ug.input -Grid.Overlap 1")

# Figure 2a
p_, z1_ = read3D_data("benchmarkUG_1a-00001", np_, 2)
h1_ = vg.pa2head(p_)
ax1.plot(h1_, (z1_[:, 2] - 2) * 100, "r+")

# Figure 2b
p_, z2_ = read3D_data("benchmarkUG_1b-00001", np_, 2)
h2_ = vg.pa2head(p_)
ax2.plot(h2_, (z2_[:, 2] - 2) * 100, "r+")

# Figure 2c
p_, z3_ = read3D_data("benchmarkUG_1c-00001", np_, 2)
h3_ = vg.pa2head(p_)
ax3.plot(h3_, (z3_[:, 2] - 2) * 100, "r+")

np.savetxt("dumuxUG_b1", np.vstack((z1_[:, 2] - 2, h1_, z2_[:, 2] - 2, h2_, z3_[:, 2] - 2, h3_)), delimiter = ",")

plt.show()

