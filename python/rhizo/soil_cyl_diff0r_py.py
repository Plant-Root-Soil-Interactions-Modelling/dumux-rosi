import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src")

import functional.van_genuchten as vg
from fv.fv_grid import *
import fv.fv_advectiondiffusion as ad  # Python solver
import fv.fv_richards as richards  # Python solver

import numpy as np
from scipy import sparse
import scipy.sparse.linalg as LA
import matplotlib.pyplot as plt
import os
import time

""" 
Cylindrical 1D model (Pyhton) Diffusion only, zero sink, no water movement

richards model is not solved (since it does nothing), otherwise use fv_system to solve both
"""
ndof = 200
# nodes = np.logspace(np.log10(0.02), np.log10(0.6), ndof + 1)
nodes = np.linspace(0.02, 0.6, ndof + 1)
grid = FVGrid1Dcyl(nodes)

loam = [0.045, 0.43, 0.04, 1.6, 50]
soil = vg.Parameters(loam)
theta = vg.water_content(-100., soil)
rich = richards.FVRichards1D(grid, loam)  # specialized 1d solver (direct banded is sufficient)
rich.x0 = np.ones((ndof,)) * (-100)  # [cm] initial soil matric potential

ad = ad.FVAdvectionDiffusion_richards(grid, rich)
ad.x0 = np.ones((ndof,)) * 0.01  # [g/cm] initial concentration
ad.b0 = np.ones((ndof,)) * 140  # [1] buffer power
ad.D0 = np.ones((ndof,)) * 1.e-5 * 24.* 3600.  # [cm2/day]
ad.bc[(0, 0)] = ["concentration", [0., grid.center(0), np.array([-1])]]

sim_times = [ 10., 20.]  # days 25, 30
maxDt = 0.01

t = time.time()
c = ad.solve(sim_times, maxDt)
print("elapsed time", time.time() - t)

col = ["r*", "g*", "b*", "c*", "m*", "y*", ]
for i in range(0, len(sim_times)):
    plt.plot(ad.grid.centers(), c[i,:], col[i], label = "Time {:g} days".format(sim_times[i]))
plt.xlabel("cm")
plt.ylabel("solute concentration (g/cm3)")
data = np.loadtxt("c_diff_results.txt", skiprows = 8)
z_comsol = data[:, 0]
plt.plot(z_comsol + 0.02, data[:, 25], "k", label = "comsol 10 days")
plt.plot(z_comsol + 0.02, data[:, -1], "k:", label = "comsol 20 days")
plt.xlabel("distance from root axis (cm)")
plt.legend()
plt.show()

