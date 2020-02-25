import sys
sys.path.append("../../../build-cmake/rosi_benchmarking/python_solver/")

from solver.xylem_flux import XylemFluxPython  # Python hybrid solver
import solver.plantbox as pb
import solver.rsml_reader as rsml
from dumux_rosi import RichardsSP  # C++ part (Dumux binding)
from solver.richards import RichardsWrapper  # Python part

import van_genuchten as vg

from math import *
import numpy as np
import matplotlib.pyplot as plt
import timeit
import cProfile

from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

""" 
Benchmark M1.1 single root in in soil

also works parallel with mpiexec (only slightly faster, due to overhead)
"""

cpp_base = RichardsSP()
s = RichardsWrapper(cpp_base)
s.initialize()


def solve(soil, simtimes, q_r, N):
    """ 
    soil            soil type
    simtimes        simulation output times
    q_r             root flux [cm/day]   
    N               spatial resolution
    """

    r_root = 0.02  # cm
    r_out = 0.6

    q_r = q_r * 1  # * 1 g/cm3 = g/cm2/day
    q_r = q_r * 2 * r_root * np.pi * 1.  # g/day
    print("Qr as sink", q_r, "g day-1")
    trans = q_r

    """ Root problem"""
    n1 = pb.Vector3d(0, 0, 0)
    n2 = pb.Vector3d(0, 0, -1)
    seg = pb.Vector2i(0, 1)
    rs = pb.MappedSegments([n1, n2], [seg], [r_root])  # a single root
    r = XylemFluxPython(rs)
    r.setKr([1.e-7])
    r.setKx([1.e-7])

    """ Soil problem """
    s.setHomogeneousIC(-100)  # cm pressure head
    s.setTopBC("noflux")
    s.setBotBC("noflux")
    l = np.sqrt((r_out * r_out - r_root * r_root) * np.pi) / 2  # same area as cylindrical
    s.createGrid([-l, -l, -1.], [l, l, 0.], [N, N, 1])  # [cm]
    s.setVGParameters([soil[0:5]])
    s.initializeProblem()

    """ Coupling (map indices) """
    picker = lambda x, y, z : s.pick(x, y, z)
    r.rs.setSoilGrid(picker)
    cci = picker(0, 0, 0)  # collar cell index

    """ Numerical solution """
    start_time = timeit.default_timer()
    nsp = N  # number of sample points for output
    y_ = np.linspace(0, l, nsp)
    y_ = np.expand_dims(y_, 1)
    x_ = np.hstack((np.zeros((nsp, 1)), y_, np.zeros((nsp, 1))))
    sx = s.getSolutionHead()  # inital condition, solverbase.py

    dt_ = np.diff(simtimes)
    s.ddt = 1.e-5  # initial Dumux time step [days]

    rx_hom = r.solve_neumann(0., -trans)  # xylem_flux.py
    # homogeneous solution is constant and indepent of soil

    for dt in dt_:

        if rank == 0:  # Root part is not parallel
            rx = r.getSolution(rx_hom, sx)  # class XylemFlux is defined in MappedOrganism.h
            fluxes = r.soilFluxes(0., rx_hom)  # class XylemFlux is defined in MappedOrganism.h  !CARE only approx (1 seg)
            cflux = r.collar_flux(0., rx, sx)
            # print("fluxes", fluxes, " = ", trans, "=", cflux, "g day-1")
        else:
            fluxes = None
        fluxes = comm.bcast(fluxes, root = 0)  # Soil part runs parallel

        s.setSource(fluxes)  # g day-1, richards.py
        s.solve(dt)

        x0 = s.toHead(s.getSolutionAt(cci))
        if x0 < -15000:
            if rank == 0:
                print("Simulation time at -15000 cm > {:.3f} cm after {:.3f} days".format(float(x0), s.simTime))
            y = s.interpolateNN(x_)
            return y, x_[:, 1], s.simTime

    y = s.toHead(s.interpolateNN(x_))
    return y, x_[:, 1], s.simTime


if __name__ == "__main__":

    sand = [0.045, 0.43, 0.15, 3, 1000, "Sand"]
    loam = [0.08, 0.43, 0.04, 1.6, 50, "Loam"]
    clay = [0.1, 0.4, 0.01, 1.1, 10, "Clay"]

    sim_times = np.linspace(0, 25, 250)  # temporal resolution of 0.1 d
    fig, ax = plt.subplots(2, 3, figsize = (14, 14))

    if rank == 0:
        t0 = timeit.default_timer()

    jobs = ([sand, 0.1, 0, 0], [loam, 0.1, 0, 1], [clay, 0.1, 0, 2], [sand, 0.05, 1, 0], [loam, 0.05, 1, 1], [clay, 0.05, 1, 2])
    for soil, qj, i, j in jobs:
        y, x, t = solve(soil, sim_times, qj, 40)
        if rank == 0:
            ax[i, j].plot(x, s.toHead(y), "r*")
            ax[i, j].set_xlabel("r (cm)")
            ax[i, j].set_ylabel("water potential (cm)")
            ax[i, j].legend(["final: {:.3f} d".format(t)])
            ax[i, j].title.set_text(soil[5] + ", q = {:.2f} cm/d".format(qj))

    if rank == 0:
        print("Elapsed time: ", timeit.default_timer() - t0, "s")
        plt.show()

