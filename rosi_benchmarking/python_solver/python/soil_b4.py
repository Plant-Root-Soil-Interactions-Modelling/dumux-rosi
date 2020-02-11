import sys
sys.path.append("../../../build-cmake/rosi_benchmarking/python_solver/")

from richardssp import RichardsSP  # C++ part (Dumux binding)
from solver.richards import RichardsWrapper  # Python part

import matplotlib.pyplot as plt
import numpy as np

from analytic_b4 import *  # plots the analytical solutions to ax1, ax2, ax3

from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

""" 
Dynamic evaporation (Benchmark M2.2) with a 3D SPGrid but low resolution in x and y (for speed), 

everything scripted, no input file needed

also works parallel with mpiexec
"""


def solve(soil, simtime, evap, NZ, ic = -200):
    cpp_base = RichardsSP()
    s = RichardsWrapper(cpp_base)
    s.initialize()

    s.setHomogeneousIC(ic)  # cm pressure head
    s.setTopBC("atmospheric", 0.5, [[0., 1.e10], [evap, evap]])  #  [cm/day] atmospheric is with surface run-off
    s.setBotBC("freeDrainage")  # BC
    s.createGrid([-5., -5., -100.], [5., 5., 0.], [1, 1, NZ])  # [cm]
    s.setVGParameters([soil])
    s.initializeProblem()
    idx_top = s.pickCell([0.0, 0.0, 0.0])  # index to watch flux

    N = 200
    dt = simtime / N
    s.ddt = 1.e-5  # initial Dumux time step [days]
    maxDt = 1  # maximal Dumux time step [days]

    x_, y_ = [], []

    for i in range(0, N):

        if rank == 0:
            print("***** external time step", dt, " d, simulation time", s.simTime, "d, internal time step", s.ddt, "d")

        s.solve(dt, maxDt)

        f = s.getNeumann(idx_top)

        if rank == 0:
            x_.append(s.simTime)
            y_.append(f)

    return x_, y_


if __name__ == "__main__":

    sand = [0.045, 0.43, 0.15, 3, 1000]
    loam = [0.08, 0.43, 0.04, 1.6, 50]
    clay = [0.1, 0.4, 0.01, 1.1, 10]

    # low res
    xa, za = solve(sand, 1, -0.1, 99, -40)
    xb, zb = solve(loam, 10, -0.1, 99)
    xc, zc = solve(loam, 2, -0.3, 99)
    xd, zd = solve(clay, 6, -0.3, 99)

    # high res
    xa2, za2 = solve(sand, 1, -0.1, 399, -40)
    xb2, zb2 = solve(loam, 10, -0.1, 399)
    xc2, zc2 = solve(loam, 2, -0.3, 399)
    xd2, zd2 = solve(clay, 6, -0.3, 399)

    if rank == 0:
        ax1.plot(xa, za, "r")
        ax2.plot(xb, zb, "r")
        ax3.plot(xc, zc, "r")
        ax4.plot(xd, zd, "r")
        ax1.plot(xa2, za2, "g--")
        ax2.plot(xb2, zb2, "g--")
        ax3.plot(xc2, zc2, "g--")
        ax4.plot(xd2, zd2, "g--")
        plt.show()
