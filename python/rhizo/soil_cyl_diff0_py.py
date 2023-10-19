import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src")

import functional.van_genuchten as vg
from fv.fv_grid import *
import fv.fv_advectiondiffusion as ad  # Python solver

import numpy as np
from scipy import sparse
import scipy.sparse.linalg as LA
import matplotlib.pyplot as plt
import os
import time

""" 
Cylindrical 1D model (Pyhton) Diffusion only, zero sink, no water movement (no modell)
"""

loam = [0.045, 0.43, 0.04, 1.6, 50]
soil = vg.Parameters(loam)
theta = vg.water_content(-100., soil)
phi = loam[1]
sw = theta / phi

ndof = 1000
# nodes = np.logspace(np.log10(0.02), np.log10(0.6), ndof + 1)
nodes = np.linspace(0.02, 0.6, ndof + 1)
# grid = FVGrid1D(nodes)
grid = FVGrid1Dcyl(nodes)

ad = ad.FVAdvectionDiffusion1D(grid)
ad.bc[(0, 0)] = ["concentration", [0., grid.center(0), np.array([-1])]]

D = np.ones((ndof,)) * 1.e-5 * 24.* 3600.  # [cm2/day]
b = np.ones((ndof,)) * (140 + theta)  # [1] buffer power
c0 = np.ones((ndof,)) * 0.01  # [g/cm] initial concentration
u = np.zeros((ndof,))  # [cm/day] velocity field
ad.D = phi * (sw ** 3) * np.cbrt(theta) * D  # D * theta * 0.5 * 0.25  # phi * (sw ** 3) * np.cbrt(theta) * D
ad.b = b
ad.x0 = c0

sim_times = [ 10., 20.]  # days 25, 30
maxDt = 0.001
print("CFL", ad.cfl(maxDt))
print("Suggested maximal time step (for c=0.1)", ad.cfl_dt(0.1))

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

