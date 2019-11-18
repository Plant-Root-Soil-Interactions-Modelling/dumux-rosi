import sys
sys.path.append("../../../build-cmake/rosi_benchmarking/python_solver/")
from richardsyaspsolver import *

import os
import time

from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

s = RichardsYaspSolver()

# s.initialize([""])
  # wrong param group (todo)
# s.createGrid([-1, -1, -1], [1,1,1], [2, 2, 2], "false, false, false")
# s.createGrid("../grids/b1.dgf") # does only take intervals

s.initialize(["-input", "../input/b1a_3d.input"])  # , "-Grid.Overlap", "0"
s.createGrid()
s.initializeProblem()

if rank == 0:

    print("\nNumber of points", len(s.getPoints()))
    print("Number of cells", len(s.getCellCenters()))

    print()
    print(s, "\n")

    print("Bounding box ", s.getGridBounds())  # it is not optimal to pass parameters via strings

#     print()
#     print("Initial total water ", s.getgetWaterVolume())

    print()
    coord = s.getDofCorrdinates()
    print("DOF ", coord.shape[0] - 1)
    sol = s.getSolution()

t = time.time()

# dt = 3600 * 24  # a day in seconds
# s.ddt = 360  # s, initial internal time step
# for i in range(0, 10):
#     if rank == 0:
#         print("*************** External time step ", i, dt, "****************", "Simulation time ", s.simTime / 3600 / 24, "days, internal time step", s.ddt / 3600 / 24, "days")
#     s.simulate(dt, -1)
    #
    # do wild stuff
    #

if rank == 0:
    print("\nOverall time ", time.time() - t)
    print("done")

# print(len(s.initialValues[0]))
# print(s.initialValues[0])
# print(len(s.solution[0]))
# print(s.solution)

