import sys
sys.path.append("../../../build-cmake/rosi_benchmarking/python_solver/")

from richardssp import RichardsSP  # C++ part (Dumux binding)
from solver.richards import RichardsWrapper  # Python part

import matplotlib.pyplot as plt

from analytic_b2 import *  # plots the analytical solutions to ax1, ax2, ax3

from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

""" 
Steady state evaporation with a 3D SPGrid but no resolution in x and y (for speed)

everything scripted, no input file needed

also works parallel with mpiexec
"""

cpp_base = RichardsSP()
s = RichardsWrapper(cpp_base)
s.initialize()

N = 52  # resolution
s.createGrid([-0.05, -0.05, -.53], [0.05, 0.05, 0], [1, 1, N])  # [m]
s.setICZ([0, -53] , [-.53, 0])  # inital guess
s.setTopBC("constantFlux", -0.5)  #  [cm/day]
s.setBotBC("constantPressure", 0.)
loam = [0.08, 0.43, 0.04, 1.6, 50]
s.setVGParameters([loam])
s.initializeProblem()

dt = 3600 * 24  # a day in seconds
simtime = 356  # days
s.ddt = 1  # s, initial internal time step

for i in range(0, simtime):

    if rank == 0:
        print("*****", "#", i, "external time step", dt / 3600 / 24, " days, simulation time",
              s.simTime / 3600 / 24, "days, internal time step", s.ddt / 3600 / 24, "days")

    s.simulate(dt)

z = s.getDofCoordinates()
x = s.getSolution()

if rank == 0:
    print("\n", s)  # some info
    plt.plot(RichardsWrapper.toHead(x), z[:, 2] * 100, "r*")
    plt.show()

