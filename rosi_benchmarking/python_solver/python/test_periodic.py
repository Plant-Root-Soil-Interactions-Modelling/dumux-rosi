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

s.initialize([""])
s.setParameter("Problem.Name", "periodicity")
s.createGrid([-0.25, -0.25, -0.5], [0.25, 0.25, 0.], [49, 49, 49], "false true true")  # 125000
s.setVanGenuchtenParameter(0.08, 0.43, 0.04, 1.6, 50.)  # Loam
if rank == 0:
    print("trying", s.getGridBounds())
s.setHomogeneousInitialConditions(-100, True)  # cm pressure head, hydraulic equilibrium
s.setBCTopBot("constantFlux", 0, "freeDrainage", 0.)
s.initializeProblem()

points = s.getDofCoordinates()
id = s.pickCell([0., 0.2, -0.1])
if rank == 0:
    dof = points.shape[0]
    source = np.zeros((dof, 1))
    print("DOF", dof, "Picked id", id)
source[dof] = 5  # g / day
s.setSource(source)

dt = 3600  # ten days in seconds
s.ddt = 1  # s, initial internal time step

for i in range(0, 10 * 24):
    if rank == 0:
        print(i, "*** External time step ", dt, "***", "Simulation time ", s.simTime / 3600 / 24)
    s.simulate(dt)

x = np.array(s.getSolution())
# x = np.array(s.getSaturation())

if rank == 0:
    plt.plot(toHead(x), points[:, 2] * 100, "r*")
    plt.show()
