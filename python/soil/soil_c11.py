import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src")

from rosi_richards import RichardsSP  # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part

import matplotlib.pyplot as plt
import numpy as np
import time

from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()

""" 
Root uptake with the classic sink term (Benchmark C11) with a 3D SPGrid but low resolution z (for speed), 

everything scripted, no input file needed

also works parallel with mpiexec (only slightly faster, due to overhead)
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

    q_r = q_r * 1  # * 1 g/cm3 = g/cm2/day
    q_r = q_r * 2 * r_root * np.pi * 1.  # g/day
    print("Qr as sink", q_r)

    cpp_base = RichardsSP()
    s = RichardsWrapper(cpp_base)
    s.initialize()
    s.setHomogeneousIC(-100)  # cm pressure head
    s.setTopBC("noflux")
    s.setBotBC("noflux")
    l = np.sqrt((r_out * r_out - r_root * r_root) * np.pi) / 2  # same area as cylindrical
    s.createGrid([-l, -l, -1.], [l, l, 0.], [N, N, 1])  # [cm]
    s.setVGParameters([soil[0:5]])
    s.initializeProblem()
    s.setCriticalPressure(-15000)
    idx_top = s.pickCell([0., 0., -.5])  # index for sink

    sources = { idx_top:-q_r }  # gIdx: value [ g/day ]
    s.setSource(sources)

    nsp = N  # number of sample points for output
    y_ = np.linspace(0, l, nsp)
    y_ = np.expand_dims(y_, 1)
    x_ = np.hstack((np.zeros((nsp, 1)), y_, np.zeros((nsp, 1))))

    dt_ = np.diff(simtimes)
    s.ddt = 1.e-5  # initial Dumux time step [days]

    for dt in dt_:

        vol = s.getWaterVolume()  # don't move to output, needs all ranks

        if rank == 0:
            print("***** external time step {:.3f} d, simulation time {:.3f} d, internal time step {:.3f} d Water volume {:.3f} cm3".
                  format(dt, s.simTime, s.ddt, vol))

        s.solve(dt)

        x0 = s.to_head(s.getSolutionAt(idx_top))
        if x0 < -15000:
            if rank == 0:
                print("Simulation time at -15000 cm > {:.3f} cm after {:.3f} days".format(float(x0), s.simTime))
            y = s.to_head(s.interpolateNN(x_))
            return y, x_[:, 1], s.simTime

    y = s.to_head(s.interpolateNN(x_))
    return y, x_[:, 1], s.simTime


if __name__ == "__main__":

    sand = [0.045, 0.43, 0.15, 3, 1000, "Sand"]
    loam = [0.08, 0.43, 0.04, 1.6, 50, "Loam"]
    clay = [0.1, 0.4, 0.01, 1.1, 10, "Clay"]

    sim_times = np.linspace(0, 25, 250)  # temporal resolution of 0.01 d
    fig, ax = plt.subplots(2, 3, figsize = (14, 14))

    if rank == 0:
        t0 = time.time()

    t_ = []
    jobs = ([sand, 0.1, 0, 0], [loam, 0.1, 0, 1], [clay, 0.1, 0, 2], [sand, 0.05, 1, 0], [loam, 0.05, 1, 1], [clay, 0.05, 1, 2])
    for soil, qj, i, j in jobs:
        y, x, t = solve(soil, sim_times, qj, 41)
        t_.append(t)
        if rank == 0:
            ax[i, j].plot(x, y, "r*")
            ax[i, j].set_xlabel("r (cm)")
            ax[i, j].set_ylabel("water potential (cm)")
            ax[i, j].legend(["final: {:.3f} d".format(t)])
            ax[i, j].title.set_text(soil[5] + ", q = {:.2f} cm/d".format(qj))

    if rank == 0:
        print("Elapsed time: ", time.time() - t0, "s")
        print(t_)
        plt.show()

