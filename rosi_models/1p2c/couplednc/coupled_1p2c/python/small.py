# Run decoupled senarios

import os
import matplotlib.pyplot as plt
from vtk_tools import *
import van_genuchten as vg

# go to the right place
path = os.path.dirname(os.path.realpath(__file__))
os.chdir(path)
os.chdir("../../../build-cmake/rosi_benchmarking/coupled")

np_ = 1  # number of processors

# run dumux
if np_ == 1:
    os.system("./coupled_seq input/small.input")
else:
    os.system("mpirun -n " + str(np_) + " ./coupled_seq input/small.input -Grid.Overlap 0")

# Figure
s_, p_, z1_ = read3D_vtp("smallS-00001", np_)
# s2_, p2_, z2_ = read3D_vtp("smallR-00001", 1)
p2_, z2_ = read3D_vtp_data("smallR-00001.vtp", False)

h1_ = vg.pa2head(p2_)
plt.plot(h1_, z2_ * 100, "r+")
plt.xlabel('$\psi$ (cm)')
plt.ylabel('Depth (cm)')

plt.show()

