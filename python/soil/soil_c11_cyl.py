import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src")

from rosi_richards_cyl import RichardsCylFoam  # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part

import matplotlib.pyplot as plt
import numpy as np
import time
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()

""" 
Root uptake with a sink term (Benchmark C11) with a 1D cylindrical model (Dumux Solver)

everything scripted, no input file needed
"""

r_root = 0.02  # cm
r_out = 0.6


def solve(soil, simtimes, q_r, N):
    """ 
    soil            soil type
    simtimes        simulation output times
    q_r             root flux [cm/day]   
    N               spatial resolution
    """

    cpp_base = RichardsCylFoam()
    s = RichardsWrapper(cpp_base)

    s.initialize()
    s.setHomogeneousIC(-100)  # cm pressure head
    s.setOuterBC("fluxCyl")
    s.setInnerBC("fluxCyl", -qj)
    s.setVGParameters([soil[0:5]])
    s.createGrid1d(np.linspace(r_root, r_out, N))  # [cm]
    # s.createGrid([0.02], [0.6], [N])  # [cm]

    s.initializeProblem()
    s.setCriticalPressure(-np.inf)  # disable switch to dirichlet at critical pressure, otherwise we do not notice, when the switch happens
    idx_top = s.pickCell([r_root])  # index for sink
    x_ = s.getDofCoordinates()
    # s.setRegularisation(1.e-6, 1.e-6)

    dt_ = np.diff(simtimes)
    s.ddt = 1.e-5  # initial Dumux time step [days]

    for dt in dt_:

        if rank == 0:
            print("***** external time step {:.3f} d, simulation time {:g} d, internal time step {:g} d Water volume cm3".format(dt, s.simTime, s.ddt))

        s.solve(dt)

        x0 = s.getSolutionHeadAt(idx_top)
        if x0 < -15000:
            y = s.getSolutionHead()
            return y, x_, s.simTime

    y = s.getSolutionHead()
    return y, x_, s.simTime


if __name__ == "__main__":

    sand = [0.045, 0.43, 0.15, 3, 1000, "Sand"]
    loam = [0.08, 0.43, 0.04, 1.6, 50, "Loam"]
    clay = [0.1, 0.4, 0.01, 1.1, 10, "Clay"]

    sim_times = np.linspace(0, 25, 250)  # temporal resolution of 0.1 d
    fig, ax = plt.subplots(2, 3, figsize = (14, 14))

    if rank == 0:  # measure simulation wall time
        t0 = time.time()

    jobs = ([sand, 0.1, 0, 0], [loam, 0.1, 0, 1], [clay, 0.1, 0, 2], [sand, 0.05, 1, 0], [loam, 0.05, 1, 1], [clay, 0.05, 1, 2])
    for soil, qj, i, j in jobs:
        y, x, t = solve(soil, sim_times, qj, 41)
        if rank == 0:
            ax[i, j].plot(x, y, "r*")
            ax[i, j].set_xlabel("r (cm)")
            ax[i, j].set_ylabel("water potential (cm)")
            ax[i, j].legend(["final: {:.3f} d".format(t)])
            ax[i, j].title.set_text(soil[5] + ", q = {:.2f} cm/d".format(qj))

    if rank == 0:  # measure simulation wall time
        print("Elapsed time: ", time.time() - t0, "s")
        plt.show()

