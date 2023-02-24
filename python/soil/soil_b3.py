import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src")

from rosi_richards import RichardsSP  # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part
from analytic_b3 import *  # plots the analytical solutions to ax1, ax2, ax3

import matplotlib.pyplot as plt
import numpy as np
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()

""" 
Dynamic infiltration (Benchmark M2.1) with a 3D SPGrid but no resolution in x and y (for speed), 

everything scripted, no input file needed

also works parallel with mpiexec
"""


def solve(soil, simtimes):
    s = RichardsWrapper(RichardsSP())
    s.initialize()
    s.setTopBC("atmospheric", 0.5, [[0., 1.e10], [100., 100.]])  #  [cm/day] atmospheric is with surface run-off
    s.setBotBC("freeDrainage")
    N = 199
    s.createGrid([-5., -5., -200.], [5., 5., 0.], [1, 1, N])  # [cm] N
    s.setHomogeneousIC(-400.)  # cm pressure head
    s.setVGParameters([soil])
    s.initializeProblem()

    dt_ = np.diff(simtimes)
    s.ddt = 1.e-5  # initial dumux time step [days]

    z_, x_ = [], []
    for dt in dt_:

        if rank == 0:
            print("***** external time step", dt, " d, simulation time", s.simTime, "d, internal time step", s.ddt, "d")

        s.solve(dt)

        points = s.getDofCoordinates()
        theta = s.getWaterContent()

        if rank == 0:
            z_.append(points[:, 2])
            x_.append(theta)

    return x_, z_


if __name__ == "__main__":

    sand = [0.045, 0.43, 0.15, 3, 1000]
    loam = [0.08, 0.43, 0.04, 1.6, 50]
    clay = [0.1, 0.4, 0.01, 1.1, 10]

    xa, za = solve(sand, [0., 0.1, 0.2, 0.3])
    xb, zb = solve(loam, [0., 0.2, 0.5, 1.])
    xc, zc = solve(clay, [0., 0.1, 0.2, 0.5])

    if rank == 0:
        for i in range(0, 3):
            ax1.plot(xa[i], za[i], "r*")
            ax2.plot(xb[i], zb[i], "r*")
            ax3.plot(xc[i], zc[i], "r*")
        plt.show()

