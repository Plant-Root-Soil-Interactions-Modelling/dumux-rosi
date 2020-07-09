import sys
sys.path.append("../../../build-cmake/rosi_benchmarking/python_solver/")

import numpy as np
from scipy import sparse
import scipy.sparse.linalg as LA
import matplotlib.pyplot as plt

import solver.van_genuchten as vg
from solver.fv_grid import *
import solver.fv_advectiondiffusion as ad  # Python solver
import solver.fv_richards as rich  # Python solver

import matplotlib.pyplot as plt
import os
import time

""" 
Cylindrical 1D model (Pyhton) Advection and Diffusion 
"""

ndof = 100
nodes = np.linspace(0.02, 0.6, ndof + 1)
grid = FVGrid1Dcyl(nodes)

loam = [0.045, 0.43, 0.04, 1.6, 50]
soil = vg.Parameters(loam)
theta = vg.water_content(-100., soil)
print(theta)

# richards = rich.FV_Richards(grid, loam)  # general solver using UMFPack
richards = rich.FVRichards1D(grid, loam)  # specialized 1d solver (direct banded is sufficient)
richards.x0 = np.ones((ndof,)) * (-100)
dx = 0.5 * (nodes[1] - 0.01)
richards.bc[(0, 0)] = ["flux_out", [-0.1, -15000, dx]]

ad = ad.FVAdvectionDiffusion_richards(grid, richards)
dx = grid.center(0) - nodes[0]
ad.bc[(0, 0)] = ["concentration", [0., 2 * dx]]
D = np.ones((ndof,)) * 1.e-5 * 24.* 3600.  # [cm2/day]
b = np.ones((ndof,)) * (140 + theta)  # [1] buffer power # should be + theta !?
c0 = np.ones((ndof,)) * 0.01  # [g/cm] initial concentration
u = np.zeros((ndof,))  # [cm/day] velocity field
ad.D = D * theta * 0.5 * (0.5 * 0.5)  # DiffusivityConstantTortuosity: porosity * saturation * tau * diffCoeff; todo (???????)!!!!!!!!!!!!!!!!!1
ad.b = b
ad.x0 = c0

sim_times = [ 10., 20.]  # days 25, 30
maxDt = 0.01

t = time.time()
c = ad.solve(sim_times, maxDt)
print("elapsed time", time.time() - t)

fig, (ax1, ax2) = plt.subplots(1, 2)
col = ["r*", "g*", "b*", "c*", "m*", "y*", ]

for i in range(0, len(sim_times)):
    plt.plot(ad.grid.centers(), c[i, :], col[i], label = "Time {:g} days".format(sim_times[i]))
plt.xlabel("cm")
plt.ylabel("solute concentration (g/cm3)")

os.chdir("../../../build-cmake/rosi_benchmarking/soil_richardsnc/python")
data = np.loadtxt("c_diff_results.txt", skiprows = 8)
z_comsol = data[:, 0]
plt.plot(z_comsol + 0.02, data[:, 25], "k", label = "comsol 10 days")
plt.plot(z_comsol + 0.02, data[:, -1], "k:", label = "comsol 20 days")

plt.xlabel("distance from root axis (cm)")
plt.legend()
plt.show()

