import sys
sys.path.append("../../../build-cmake/rosi_benchmarking/python_solver/")

from richardssp import RichardsSP  # C++ part (Dumux binding)
from solver.richards import RichardsWrapper  # Python part

import matplotlib.pyplot as plt

from analytic_b1 import *  # plots the analytical solutions to ax1, ax2, ax3

from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

""" 
Steady state infiltration with a 3D SPGrid but no resolution in x and y (for speed), 
approximated with simulation time of one year

everything scripted, no input file needed

also works parallel with mpiexec
"""


def solve(soils):
    cpp_base = RichardsSP()
    s = RichardsWrapper(cpp_base)
    s.initialize()

    s.setHomogeneousIC(-50, False)
    s.setTopBC("constantFlux", 0.5)  #  [cm/day]
    s.setBotBC("freeDrainage")
    s.createGrid([-0.05, -0.05, -2], [0.05, 0.05, 0], [1, 1, 199])
    s.setLayersZ([2, 2, 1, 1], [-2., -0.5, -0.5, 0.])
    s.setVGParameters(soils)
    s.initializeProblem()

    if rank == 0:
        print(s)

    dt = 3600 * 24  # a day in seconds
    simtime = 356  # days
    s.ddt = 1  # s, initial internal time step

    for i in range(0, simtime):

        if rank == 0:
            print("*****", "#", i, "external time step", dt / 3600 / 24, " days, simulation time",
                  s.simTime / 3600 / 24, "days, internal time step", s.ddt / 3600 / 24, "days")

        s.simulate(dt)

    points = s.getDofCoordinates()
    x = s.getSolution()
    return x, points


if __name__ == "__main__":

    sand = [0.045, 0.43, 0.15, 3, 1000 ]
    loam = [0.08, 0.43, 0.04, 1.6, 50]
    clay = [0.1, 0.4, 0.01, 1.1, 10 ]

    xa, za = solve([loam, sand])
    xb, zb = solve([sand, loam])
    xc, zc = solve([clay, sand])

    if rank == 0:
        ax1.plot(RichardsWrapper.toHead(xa), za[:, 2] * 100, "r*")
        ax2.plot(RichardsWrapper.toHead(xb), zb[:, 2] * 100, "r*")
        ax3.plot(RichardsWrapper.toHead(xc), zc[:, 2] * 100, "r*")
        plt.show()

