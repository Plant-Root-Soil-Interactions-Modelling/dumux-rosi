import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src")

import functional.van_genuchten as vg
from fv.fv_grid import *
import fv.fv_richards as richards  # Python solver

import numpy as np
from scipy import sparse
import scipy.sparse.linalg as LA
import matplotlib.pyplot as plt
import os
import time

""" 
Cylindrical 1D model (Python solver) Water movement only

higher accuracy results start to differ more (?)
"""

ndof = 50  # ndof = 1000, dann passts
nodes = np.logspace(np.log10(0.02), np.log10(0.6), ndof + 1)  # small spatial resolution requires smaller time step
# nodes = np.linspace(0.02, 0.6, ndof + 1)
grid = FVGrid1Dcyl(nodes)

loam = [0.045, 0.43, 0.04, 1.6, 50]

sim_times = [0., 0.001, 0.1, 0.2, 1, 2, 10, 20]  # [10, 20]
# rich = richards.FV_Richards(grid, loam)  # general solver using UMFPack
rich = richards.FVRichards1D(grid, loam)  # specialized 1d solver (direct banded is sufficient)
rich.x0 = np.ones((ndof,)) * (-100)
rich.bc[(0, 0)] = ["flux_out", [-0.1, -15000, grid.center(0)]]  # dx = grid.centers[0], distance between root mid 0, and first cell mid
# rich.bc[(0, 0)] = ["rootsystem", [-1500, 100000]]

t = time.time()
h = rich.solve(sim_times, 0.001, True)
print("elapsed time", time.time() - t)

# u = rich.darcy_velocity()
# plt.title("Darcy velocities (cm/day)")
# plt.plot(rich.grid.centers(), u, "b*")
# plt.show()

col = ["r*", "g*", "b", "m", "c", "y"] * 5
for i in range(0, len(sim_times)):
    plt.plot(rich.grid.centers(), h[i,:], col[i], label = "Time {:g} days".format(sim_times[i]))
plt.xlabel("cm")
plt.ylabel("matric potential (cm)")
plt.legend()

data = np.loadtxt("cylinder_1d_Comsol.txt", skiprows = 8)
z_comsol = data[:, 0]
plt.plot(z_comsol + 0.02, data[:, 25], "k", label = "comsol 10 days")
plt.plot(z_comsol + 0.02, data[:, -1], "k:", label = "comsol 20 days")

plt.xlabel("distance from root axis (cm)")
plt.xlabel("soil matric potential (cm)")
plt.legend()
plt.show()

