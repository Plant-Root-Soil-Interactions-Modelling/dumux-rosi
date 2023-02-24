import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src")

from analytic_b1 import *  # plots the analytical solutions to ax1, ax2, ax3
from rosi_richards import RichardsSP  # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part

import matplotlib.pyplot as plt
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()

""" 
Steady state infiltration with a 3D SPGrid but no resolution in x and y (for speed), 
approximated with simulation time of one year

everything scripted, no input file needed

also works parallel with mpiexec
"""


def solve(soils):

    s = RichardsWrapper(RichardsSP())
    s.initialize()
    s.createGrid([-5., -5., -200.], [5., 5., 0.], [1, 1, 199])  # [cm]
    s.setHomogeneousIC(-50.)  # cm pressure head
    s.setTopBC("constantFlux", 0.5)  #  [cm/day]
    s.setBotBC("freeDrainage")
    s.setLayersZ([2, 2, 1, 1], [-200., -50., -50., 0.])  # sample points ([1], [cm])
    s.setVGParameters(soils)
    s.initializeProblem()

    if rank == 0:
        print(s)

    dt = 1  # a days
    simtime = 356  # days
    s.ddt = 1.e-3  # days

    for i in range(0, simtime):

        if rank == 0:
            print("*****", "#", i, "external time step", dt, " d, simulation time", s.simTime, "d, internal time step", s.ddt, "d")

        s.solve(dt)

    points = s.getDofCoordinates()
    x = s.getSolution()

    return x, points


if __name__ == "__main__":

    sand = [0.045, 0.43, 0.15, 3, 1000]
    loam = [0.08, 0.43, 0.04, 1.6, 50]
    clay = [0.1, 0.4, 0.01, 1.1, 10]

    xa, za = solve([loam, sand])
    xb, zb = solve([sand, loam])
    xc, zc = solve([clay, sand])

    if rank == 0:
        ax1.plot(RichardsWrapper.to_head(xa), za[:, 2], "r*")
        ax2.plot(RichardsWrapper.to_head(xb), zb[:, 2], "r*")
        ax3.plot(RichardsWrapper.to_head(xc), zc[:, 2], "r*")
        plt.show()

