import sys
sys.path.append("../../../build-cmake/rosi_benchmarking/python_solver/")

from richardssp import RichardsSP  # C++ part (Dumux binding)
from solver.richards import RichardsWrapper  # Python part

import matplotlib.pyplot as plt
import numpy as np

from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

""" 
Root uptake with a sink term (Benchmark C11) with a 3D SPGrid but low resolution z (for speed), 

everything scripted, no input file needed

also works parallel with mpiexec
"""

cpp_base = RichardsSP()
s = RichardsWrapper(cpp_base)
s.initialize()


def solve(soil, simtimes, N, q_r):
    """ 
    soil            soil type
    simtimes        simulation output times
    N               sub steps per day
    q_r             root flux [cm/s]       
    """

    q_r = q_r / 100 * 1000  # kg/m2/s
    q_r = q_r * 1000 / 1.e4  # g/cm2/s = 0.5
    q_r = q_r * 2 * 0.2 * np.pi * 1  # g/day
    print("Qr as sink", q_r)

    s.setHomogeneousIC(-100)  # cm pressure head
    s.setTopBC("noflux")
    s.setBotBC("noflux")
    s.createGrid([-.971, -.971, -1.], [.971, .971, 0.], [20, 20, 3])  # [cm]
    s.setVGParameters([soil])
    s.initializeProblem()
    idx_top = s.pickCell([0.0, 0.0, -.5])  # index to watch flux

    sources = { idx_top:-q_r }  # gIdx: value [ kg/s ]
    s.setSource(sources)

    dt_ = np.diff(simtimes)
    s.ddt = 1.e-5  # initial Dumux time step [days]

    points = s.getDofCoordinates()
    x = None

    for dt in dt_:

        print("Water volume", s.getWaterVolume(), "cm3")
        if rank == 0:
            print("***** external time step", dt, " d, simulation time", s.simTime, "d, internal time step", s.ddt, "d")

        s.simulate(dt)

        x = s.toHead(s.getSolution())
        x0 = x[idx_top]  # MID

        if x0 < -15000:
            print("Simualtion time at -15000 cm", x0, "after", s.simTime, "days")
            return x, points

    return x, points


if __name__ == "__main__":

    sand = [0.045, 0.43, 0.15, 3, 1000]
    loam = [0.08, 0.43, 0.04, 1.6, 50]
    clay = [0.1, 0.4, 0.01, 1.1, 10]  # UNITS ?!!!!!qr, qs, alpha, n, ks

    clay_times = np.linspace(0, 3, 300)

    x_final, z = solve(clay, clay_times, 10, 0.05)

    if rank == 0:

        plt.plot(z[:, 1], x_final, "r*")
        plt.xlabel("cm")
        plt.ylabel("cm pressure head")
        plt.show()

