# Benchmark 1 (in 3D)
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
import time

# go to the right place
path = os.path.dirname(os.path.realpath(__file__))
os.chdir(path)
os.chdir("../../../build-cmake/rosi_benchmarking/soil_richards")

# run dumux
t = time.time()
np_ = 1  # number of processors
if np_ == 1:
    pass
    os.system("./richards3d input/b1a_3d.input")
    os.system("./richards3d input/b1b_3d.input")
    os.system("./richards3d input/b1c_3d.input")
else:
    os.system("mpirun -n " + str(np_) + " ./richards3d input/b1a_3d.input -Grid.Overlap 0")
    os.system("mpirun -n " + str(np_) + " ./richards3d input/b1b_3d.input -Grid.Overlap 0")
    os.system("mpirun -n " + str(np_) + " ./richards3d input/b1c_3d.input -Grid.Overlap 0")

elapsed = time.time() - t
print("Time elapsed", elapsed)

# Figure 2a
s_, p_, z1_ = read3D_vtp("benchmark3d_1a-00001", np_)
h1_ = vg.pa2head(p_)
ax1.plot(h1_, z1_ * 100, "r+")

# Figure 2b
s_, p_, z2_ = read3D_vtp("benchmark3d_1b-00001", np_)
h2_ = vg.pa2head(p_)
ax2.plot(h2_, z2_ * 100, "r+")

# Figure 2c
s_, p_, z3_ = read3D_vtp("benchmark3d_1c-00001", np_)
h3_ = vg.pa2head(p_)
ax3.plot(h3_, z3_ * 100, "r+")

# np.savetxt("dumux3d_b1", np.vstack((z1_, h1_, z2_, h2_, z3_, h3_)), delimiter = ",")

plt.show()

