# Bot Scenario Felicien

import os
import matplotlib.pyplot as plt
from vtk_tools import *
import van_genuchten as vg

# go to the right place
path = os.path.dirname(os.path.realpath(__file__))
os.chdir(path)
os.chdir("../../../build-cmake/rosi_benchmarking/soil")

# run dumux
np_ = 1  # number of processors
if np_ == 1:
    os.system("./richards1d benchmarks_3d/swbot.input")
else:
    os.system("mpirun -n " + str(np_) + " ./richards1d benchmarks_3d/swbot.input -Grid.Overlap 0")

# Figure
s_, p_, z1_ = read1D_vtp_data("swbot-00000.vtp", False)
h1_ = vg.pa2head(p_)
plt.plot(h1_, (z1_ - 2) * 100, "r+")

plt.show()

