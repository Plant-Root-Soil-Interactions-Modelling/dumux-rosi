# Benchmark 1 (in 3D)
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
import time

# go to the right place
path = os.path.dirname(os.path.realpath(__file__))
os.chdir(path)
os.chdir("../../../build-cmake/cpp/soil_richards")

# # run dumux
t = time.time()
np_ = 1  # number of processors
if np_ == 1:
    os.system("./richards3d input/b1a_3d.input")
#     os.system("./richards3d input/b1b_3d.input")
#     os.system("./richards3d input/b1c_3d.input")
else:
    os.system("mpiexec -n " + str(np_) + " ./richards3d input/b1a_3d.input")
#     os.system("mpirun -n " + str(np_) + " ./richards3d input/b1b_3d.input")
#     os.system("mpirun -n " + str(np_) + " ./richards3d input/b1c_3d.input")

elapsed = time.time() - t
print("Time elapsed", elapsed)

# Figure 2a
p_, z1_ = read3D_data("benchmark3d_1a-00001", np_, 2)
h1_ = vg.pa2head(p_)
print(h1_)
ax1.plot(h1_, z1_[:, 2] * 100, "r+")

# Figure 2b
p_, z2_ = read3D_data("benchmark3d_1b-00001", np_, 2)
h2_ = vg.pa2head(p_)
ax2.plot(h2_, z2_[:, 2] * 100, "r+")

# Figure 2c
p_, z3_ = read3D_data("benchmark3d_1c-00001", np_, 2)
h3_ = vg.pa2head(p_)
ax3.plot(h3_, z3_[:, 2] * 100, "r+")

np.savetxt("dumux3d_b1", np.vstack((z1_[:, 2], h1_, z2_[:, 2], h2_, z3_[:, 2], h3_)), delimiter=",")

plt.show()

