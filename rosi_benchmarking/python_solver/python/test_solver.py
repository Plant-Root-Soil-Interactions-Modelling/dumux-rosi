import sys
sys.path.append("../../../build-cmake/rosi_benchmarking/python_solver/")
from richardsyaspsolver import *

import os
import time

import numpy as np
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
import matplotlib.pyplot as plt

from mpi4py import MPI

""" in the far future this will be a unit test """


def toHead(pa):  # Pascal (kg/ (m s^2)) to cm pressure head
    g = 9.81  # gravitational acceleration (m/s^2)
    rho = 1.e3  # density of water, (kg/m^3)
    ref = 1.e5  # Pa
    return (pa - ref) * 100 / rho / g


def plot_interpolated_Z(eq = 0, conv = lambda x: x, xlabel = "Pressure", xy = (0, 0)):
    """ plots the current solution along the z axis (not working, todo) """
    self.checkInitialized()
    bounds = self.getGridBounds()
    z_ = np.linspace(bounds[2], bounds[5], 200)
    y_ = np.ones(z_.shape) * xy[1]
    x_ = np.ones(z_.shape) * xy[0]
    xi = np.hstack((x_, y_, z_))
    sol = self.interpolate(xi, eq)  # for this all processes are needed
    if rank == 0:
        plt.plot(conv(sol), z_ * 100, "*")  #
        plt.xlabel(xlabel)
        plt.ylabel("Depth (cm)")
        plt.show()


def plotZ(points, v, xlabel = "Pressure"):
    """ plots x along the z axis """
    if rank == 0:
        plt.plot(v, points[:, 2] * 100, "*")
        plt.xlabel(xlabel)
        plt.ylabel("Depth (cm)")
        plt.show()


comm = MPI.COMM_WORLD
rank = comm.Get_rank()

s = RichardsYaspSolver()  # the one and only

# s.initialize([""])
# s.createGrid([-1, -1, -1], [1,1,1], [2, 2, 2], "false, false, false")
# s.createGrid("../grids/b1.dgf") # YASP does only take intervals
s.initialize(["-input", "../input/b1a_3d.input"])  # , "-Grid.Overlap", "0"
# s.createGrid("Soil")
s.createGrid([-0.05, -0.05, -2], [0.05, 0.05, 0], [9, 9, 199], "false false false")
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

ind = s.getDofIndices()  # collect global ids
if rank == 0:
    print()
    print("Indices: ", len(ind))
    l = list(set(ind))
    print("Indices unique: ", len(l))
    l = np.sort(l)
    print("from ", l[0], "to", l[-1])
    print()

# points = s.getDofCorrdinates()  # in case of box, these are grid vertices
# if rank == 0:
#     print("Coordinates ", points.shape[0])
#     #print(points)
#     fig = plt.figure()

#     ax = fig.add_subplot(111, projection = '3d')
#     ax.set_xlabel('X Label')
#     ax.set_ylabel('Y Label')
#     ax.set_zlabel('Z Label')
#     ax.scatter(points[:, 0], points[:, 1], points[:, 2])
#     plt.show()

# sol = s.getSolution()
# if rank == 0:
#     print("Initial solution ", sol.shape[0])
# s.plotZ(0, toHead)

# bounds = s.getGridBounds()
# if rank == 0:
#     print("Bounding box ", bounds)  # it is not optimal to pass parameters via strings
#
# x = s.getSolution()
# s.plotZ(toHead(x), "Pressure head (cm)")

# input("Press Enter to continue...")

print()
t = time.time()

dt = 10 * 3600 * 24  # a day in seconds
s.ddt = 1  # s, initial internal time step
for i in range(0, 10):
    if rank == 0:
        print("*************** External time step ", i, dt, "****************", "Simulation time ", s.simTime / 3600 / 24, "days, internal time step", s.ddt / 3600 / 24, "days")
    s.simulate(dt, -1)

    #
    # optionally, do wild stuff
    #

if rank == 0:
    print("\nOverall time ", time.time() - t)

ind = s.getDofIndices()  # collect global ids
if rank == 0:
    print("Dof indices")
    print("Indices: ", len(ind))
    print("Indices unique: ", len(list(set(ind))))
    print()

bounds = s.getGridBounds()
if rank == 0:
    print("Bounding box:")
    print("Bounding box ", bounds)  # it is not optimal to pass parameters via strings

x = np.array(s.getSolution())
if rank == 0:
    print("Solution:")
    print(x.shape)
    print("first", toHead(x[0]))
    print("last", toHead(x[-1]))
    print("last", toHead(x[-2]))

points = s.getDofCoordinates()
if rank == 0:
    print("Dof coordinates")
    # print(points)
    print(len(points))
    print(type(points))
    print(points.shape)

plotZ(points, toHead(x), "Pressure head (cm)")

