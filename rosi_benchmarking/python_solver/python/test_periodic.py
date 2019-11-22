import sys
sys.path.append("../../../build-cmake/rosi_benchmarking/python_solver/")
from richardsyaspsolver import *

import os
import time

import numpy as np
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
import matplotlib.pyplot as plt

from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

s = RichardsYaspSolver()  # the one and only

#
# Set up the problem
#
s.initialize([""])
s.setParameter("Problem.Name", "periodicity")
s.createGrid([-0.25, -0.25, -0.5], [0.25, 0.25, 0.], [19, 19, 19], "true true false")  # 125000
s.setVanGenuchtenParameter(0.08, 0.43, 0.04, 1.6, 50.)  # Loam
s.setHomogeneousInitialConditions(-100, True)  # cm pressure head, hydraulic equilibrium
s.setBCTopBot("constantFlux", 0, "freeDrainage", 0.)
s.initializeProblem()

if rank == 0:
    print("\nGrid bounds", s.getGridBounds(), "\n")

points = s.getDofCoordinates()  # gathered in rank = 0
cells = s.getCellCenters()  # gathered in rank = 0

dof = 0;
if rank == 0:
    dof = points.shape[0]
dof = comm.bcast(dof, root = 0)

# Test picking
p = [0., 0.02, -0.01]
id = s.pickCell(p)
print("Total dof of rank", rank, "=", dof, "picked id", id)
comm.barrier()
if rank == 0:
    print("Picked cell ", cells[id])
    print("Distance to element center", np.linalg.norm(cells[id] - p), "cm")

p = [0., 0.25, -0.4]
id2 = s.pickCell(p)
print("Total dof of rank", rank, "=", dof, "picked id", id)
comm.barrier()
if rank == 0:
    print("Picked cell ", cells[id2])
    print("Distance to element center", np.linalg.norm(cells[id2] - p), "cm")

sources = { id: 1.e-3, id2: 1.e-3}  # gIdx: value
s.setSource(sources)

# # Show inital condition
# x = np.array(s.getSolution())
# if rank == 0:
#     print("Water volume", s.getWaterVolume(), "cm3")
# #     plt.plot(s.toHead(x), points[:, 2] * 100, "r*")
# #     plt.show()

dt = 3600  # ten days in seconds
s.ddt = 1  # s, initial internal time step

for i in range(0, 24):
    if rank == 0:
        print(i, "*** External time step ", dt, "***", "simulation time ", s.simTime / 3600 / 24)
    s.simulate(dt)
    print("Water volume", s.getWaterVolume(), "cm3")
#

s.writeDumuxDGF("periodicTest1")
# x = np.array(s.getSolution())
# # x = np.array(s.getSaturation())
#
# if rank == 0:
#     plt.plot(s.toHead(x), points[:, 2] * 100, "r*")
#     plt.show()
