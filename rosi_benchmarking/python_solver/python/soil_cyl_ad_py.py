import sys
sys.path.append("../../../build-cmake/rosi_benchmarking/python_solver/")

import numpy as np
from scipy import sparse
import scipy.sparse.linalg as LA
import matplotlib.pyplot as plt

import solver.van_genuchten as vg
from solver.fv_grid import *
import solver.fv_advectiondiffusion as ad  # Python solver
import solver.fv_richards as richards  # Python solver
import solver.fv_system as system  # Python solver

import matplotlib.pyplot as plt
import os
import time

""" 
Cylindrical 1D model (Pyhton) Advection and Diffusion 
"""

ndof = 1000
# nodes = np.logspace(np.log10(0.02), np.log10(0.6), ndof + 1)
nodes = np.linspace(0.02, 0.6, ndof + 1)
grid = FVGrid1Dcyl(nodes)

loam = [0.045, 0.43, 0.04, 1.6, 50]
soil = vg.Parameters(loam)
theta = vg.water_content(-100., soil)
rich = richards.FVRichards1D(grid, loam)  # specialized 1d solver (direct banded is sufficient)
dx = grid.nodes[1] - grid.center(0)
rich.x0 = np.ones((ndof,)) * (-100)  # [cm] initial soil matric potential
rich.bc[(0, 0)] = ["flux_out", [-0.1, -15000, 2 * dx]]

ad = ad.FVAdvectionDiffusion1D_richards(grid, rich)
ad.x0 = np.ones((ndof,)) * 0.01  # [g/cm] initial concentration
ad.b0 = np.ones((ndof,)) * 140.  # [1] buffer power
ad.D0 = np.ones((ndof,)) * 1.e-5 * 24.* 3600.  # [cm2/day]
dx = grid.nodes[1] - grid.center(0)
ad.bc[(0, 0)] = ["concentration", [0., 2 * dx, np.array([-1])]]

sim_times = [ 1, 10., 20.]  # days 25, 30
maxDt = 0.001

eqns = system.FVSystem(grid)
eqns.add_eqn(rich)  # eqn 0
eqns.add_eqn(ad)  # eqn 1

t = time.time()
c = eqns.solve(sim_times, maxDt)
print("elapsed time", time.time() - t)

fig, (ax1, ax2) = plt.subplots(1, 2)
col = ["r*", "g*", "b*", "c*", "m*", "y*", ]

for i in range(0, len(sim_times)):
    ax1.plot(ad.grid.centers(), c[i][0, :], col[i], label = "Time {:g} days".format(sim_times[i]))
for i in range(0, len(sim_times)):
    ax2.plot(ad.grid.centers(), c[i][1, :], col[i], label = "Time {:g} days".format(sim_times[i]))

os.chdir("../../../build-cmake/rosi_benchmarking/soil_richards/python")
data = np.loadtxt("cylinder_1d_Comsol.txt", skiprows = 8)
z_comsol = data[:, 0]
ax1.plot(z_comsol + 0.02, data[:, 25], "k", label = "comsol 10 days")
ax1.plot(z_comsol + 0.02, data[:, -1], "k:", label = "comsol 20 days")
ax1.set_xlabel("distance from root axis (cm)")
ax1.set_ylabel("soil matric potential (cm)")
ax1.legend()

os.chdir("../../../build-cmake/rosi_benchmarking/soil_richardsnc/python")
data = np.loadtxt("c_results.txt", skiprows = 8)  # buffer power = 100
# data = np.loadtxt("c_results_nob.txt", skiprows=8) # buffer power = 0

z_comsol = data[:, 0]
ax2.plot(z_comsol + 0.02, data[:, 2], "r", label = "1 hour")
ax2.plot(z_comsol + 0.02, data[:, 25], "g", label = "comsol 10 days")
ax2.plot(z_comsol + 0.02, data[:, -1], "b:", label = "comsol 20 days")
ax2.set_xlabel('distance from the root axis (cm)')
ax2.set_ylabel('solute concentration (g/cm3)')
ax2.legend()

plt.xlabel("distance from root axis (cm)")
plt.legend()
plt.show()

