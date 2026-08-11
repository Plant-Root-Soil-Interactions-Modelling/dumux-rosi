"""
Steady state infiltration with a 3D SPGrid but no resolution in x and y (for speed),
approximated with simulation time of one year (s.solveSteadyState() should also work, but is not tested here)

everything scripted, no input file needed

also works parallel with mpiexec
"""

import matplotlib.pyplot as plt
from analytic_b1 import *  # plots the analytical solutions to ax1, ax2, ax3
from mpi4py import MPI
from rosi.richards import RichardsWrapper  # Python part
from rosi.rosi_richards import RichardsSP  # C++ part (Dumux binding)

comm = MPI.COMM_WORLD
rank = comm.Get_rank()


def solve(soils):

    s = RichardsWrapper(RichardsSP())
    s.initialize()
    s.createGrid([-5.0, -5.0, -200.0], [5.0, 5.0, 0.0], [1, 1, 199], periodic=False)  # [cm] perid
    s.setHomogeneousIC(-50.0)  # cm pressure head
    s.setTopBC("constantFlux", 0.5)  #  [cm/day]
    s.setBotBC("freeDrainage")
    s.setLayersZ([2, 2, 1, 1], [-200.0, -50.0, -50.0, 0.0])  # sample points ([1], [cm])
    s.setVGParameters(soils)
    s.initializeProblem()

    idx = s.pick([0.0, 0.0, -0.01])
    print("Top cell index is", idx)

    if rank == 0:
        print(s)

    dt = 1  # a days
    simtime = 356  # days
    s.ddt = 1.0e-3  # days

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
