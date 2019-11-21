import sys
from _socket import AF_AX25
sys.path.append("../../../build-cmake/rosi_benchmarking/python_solver/")
from richardsyaspsolver import *

import os
import time

import numpy as np
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
import matplotlib.pyplot as plt

from analytic_b1 import *  # plots the analytical solutions to ax1, ax2, ax3

from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

""" We solve infiltration Benchmark M2.1 with a 3D YaspGrid but no resolution in x and y (for speed), 
    also works parallel with mpiexec
"""


def solve(number, s_):
    """ solve benchmark 1abc , number = 'a', 'b', or 'c'  """

    s = RichardsYaspSolver()
    s_.append(s)  # remember (for plot)
    s.initialize(["-input", "../input/b1" + number + "_3d.input"])
    s.createGrid([-0.05, -0.05, -2], [0.05, 0.05, 0], [1, 1, 199], "false false false")
    s.initializeProblem()

    t = time.time()

    dt = 10 * 3600 * 24  # a day in seconds
    s.ddt = 1  # s, initial internal time step
    for i in range(0, 10):
        if rank == 0:
            print(i, "*** External time step ", dt / 3600 / 24, " days ****", "Simulation time",
                  s.simTime / 3600 / 24, "days, internal time step", s.ddt / 3600 / 24, "days")
        s.simulate(dt)
        # optionally, do wild stuff

    points = s.getDofCoordinates()
    x = np.array(s.getSolution())
    if rank == 0:
        if number == 'a':
            ax1.plot(s.toHead(x), points[:, 2] * 100, "r*")
        elif number == 'b':
            ax2.plot(s.toHead(x), points[:, 2] * 100, "r*")
        elif number == 'c':
            ax3.plot(s.toHead(x), points[:, 2] * 100, "r*")


if __name__ == "__main__":
    s = []
    # solve('a', s)
    # solve('b', s)
    solve('c', s)
    if rank == 0:
        plt.show()  # TODO plot is messed up, results are correct

