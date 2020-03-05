import sys
sys.path.append("../../../build-cmake/rosi_benchmarking/python_solver/")

from dumux_rosi import RichardsSP  # C++ part (Dumux binding)
from solver.richards import RichardsWrapper  # Python part

import matplotlib.pyplot as plt
import numpy as np

from analytic_b4 import *  # plots the analytical solutions to ax1, ax2, ax3

from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

""" 
Benchmark M2.2
"""
xb2, zb2 = solve(loam, 10, -0.1, 399)


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

