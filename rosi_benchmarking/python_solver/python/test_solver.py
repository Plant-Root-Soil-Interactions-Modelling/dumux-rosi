import sys
sys.path.append("../../../build-cmake/rosi_benchmarking/python_solver/")
from richardsyaspsolver import *

import os
import time

from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
import matplotlib.pyplot as plt

from mpi4py import MPI

""" in the far future this will be a unit test """


def toHead(pa):  # Pascal (kg/ (m s^2)) to cm pressure head
    g = 9.81  # gravitational acceleration (m/s^2)
    rho = 1.e3  # density of water, (kg/m^3)
    ref = 1.e5  # Pa
    return (pa - ref) * 100 / rho / g


comm = MPI.COMM_WORLD
rank = comm.Get_rank()

s = RichardsYaspSolver()

# s.initialize([""])
# wrong param group (todo)
# s.createGrid([-1, -1, -1], [1,1,1], [2, 2, 2], "false, false, false")
# s.createGrid("../grids/b1.dgf") # YASP does only take intervals

s.initialize(["-input", "../input/b1a_3d.input"])  # , "-Grid.Overlap", "0"

s.createGrid()
s.initializeProblem()

# if rank == 0:
#     print("\nNumber of points", len(s.getPoints()))
#     print("Number of cells", len(s.getCellCenters()))
#     print()
#     print(s, "\n")
#
# bounds = s.getGridBounds()
# if rank == 0:
#     print("Bounding box ", bounds)  # it is not optimal to pass parameters via strings

# TODO
#     print()
#     print("Initial total water ", s.getgetWaterVolume())

# print()
# ind = s.getDofIndices()  # collect global ids
# if rank == 0:
#     print()
#     print("Indices: ", len(ind))
#     print("Indices unique: ", len(list(set(ind))))
#     print()

# points = s.getDofCorrdinates()  # in case of box, these are grid vertices
# if rank == 0:
#     print("Coordinates ", points.shape[0])
# #     fig = plt.figure()
# #     ax = fig.add_subplot(111, projection = '3d')
# #     ax.set_xlabel('X Label')
# #     ax.set_ylabel('Y Label')
# #     ax.set_zlabel('Z Label')
# #     ax.scatter(points[:, 0], points[:, 1], points[:, 2])
# #     plt.show()

# sol = s.getSolution()
# if rank == 0:
#     print("Initial solution ", sol.shape[0])
# s.plotZ(0, toHead)

t = time.time()

dt = 10 * 3600 * 24  # a day in seconds
s.ddt = 360  # s, initial internal time step
for i in range(0, 10):
    if rank == 0:
        print("*************** External time step ", i, dt, "****************", "Simulation time ", s.simTime / 3600 / 24, "days, internal time step", s.ddt / 3600 / 24, "days")
    s.simulate(dt, -1)
    #
    # do wild stuff
    #

if rank == 0:
    print("\nOverall time ", time.time() - t)
    print("done")

s.plotZ(0, toHead, "Pressure head (cm)")

# print(len(s.initialValues[0]))
# print(s.initialValues[0])
# print(len(s.solution[0]))
# print(s.solution)

