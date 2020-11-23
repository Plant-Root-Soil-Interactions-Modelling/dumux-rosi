# Benchmark 2 (in 3D)
#
# compares the dumux solution to the analytical solution (Figure 3, Vanderborght et al 2005)
#
# D. Leitner, 2018
#
import sys; sys.path.append("../../../../CPlantBox/src/python_modules/")
import sys; sys.path.append("../../../python/soil/")  # for the analytical solutions

import os
import matplotlib.pyplot as plt
from vtk_tools import *
import van_genuchten as vg
import analytic_b2

# go to the right place
path = os.path.dirname(os.path.realpath(__file__))
os.chdir(path)
os.chdir("../../../build-cmake/cpp/soil_richards")

# run dumux
np_ = 1  # number of processors
if np_ == 1:
    os.system("./richards3d input/b2_3d.input")
else:
    os.system("mpirun -n " + str(np_) + " ./richards3d input/b2_3d.input -Grid.Overlap 0")

# result dumux jan1 (Figure 2a)
p_, z_ = read3D_data("benchmark3d_2-00001", np_, 2)
h_ = vg.pa2head(p_)
plt.plot(h_, z_[:, 2] * 100, "r+")

np.savetxt("dumux3d_b2", np.vstack((z_[:, 2], h_)), delimiter = ",")

plt.show()
