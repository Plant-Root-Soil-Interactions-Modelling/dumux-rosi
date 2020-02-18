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

from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()


def sinusoidal(t):
    return np.sin(2. * pi * np.array(t) - 0.5 * pi) + 1.

""" 
Benchmark M1.2 static root system in soil
(does not work in parallel)
"""

""" Parameters """
kz = 4.32e-2  # [cm^3/day]
kr = 1.728e-4  # [1/day]
trans = 6.4  # cm3 /day
p0 = -1000  # dircichlet bc at top
sim_time = 7  # [day] for task b
wilting_point = -10000  # cm
loam = [0.08, 0.43, 0.04, 1.6, 50]

""" Root problem """
r = XylemFluxPython("../grids/RootSystem.rsml")
r.setKr([kr])
r.setKx([kz])
nodes = r.get_nodes()

""" Soil problem """
cpp_base = RichardsSP()
s = RichardsWrapper(cpp_base)
s.initialize()
s.createGrid([-4., -4., -20.], [4., 4., 0.], [16, 16, 40])  # [cm]
# s.createGrid([-4., -4., -20.], [4., 4., 0.], [8, 8, 20])  # [cm]
s.setHomogeneousIC(-669.8 - 10, True)  # cm pressure head, equilibrium
s.setTopBC("noFlux")
s.setBotBC("noFlux")
s.setVGParameters([loam])
s.initializeProblem()

""" Coupling (map indices) """
picker = lambda x, y, z : s.pick(x, y, z)
r.rs.setSoilGrid(picker)
cci = picker(nodes[0, 0], nodes[0, 1], nodes[0, 2])  # collar cell index

""" Numerical solution (a) """
start_time = timeit.default_timer()
x_, y_, w_, cpx, cps = [], [], [], [], []
sx = s.getSolutionHead()  # inital condition, solverbase.py

dt = 120. / (24 * 3600)  # [days] Time step must be very small
N = sim_time * round(1. / dt)
t = 0.

for i in range(0, N):

    if rank == 0:  # Root part is not parallel
        rx_hom = r.solve(t, -trans * sinusoidal(t), sx[cci], wilting_point)  # xylem_flux.py
        rx = r.getSolution(rx_hom, sx)  # class XylemFlux is defined in MappedOrganism.h
        fluxes = r.soilFluxes(t, rx_hom)  # class XylemFlux is defined in MappedOrganism.h
    else:
        fluxes = None

    fluxes = comm.bcast(fluxes, root = 0)  # Soil part runs parallel
    s.setSource(fluxes)  # richards.py
    s.solve(dt)

    sx = s.getSolutionHead()  # richards.py
    water = s.getWaterVolume()

    if rank == 0:
        n = round(float(i) / float(N) * 100.)
        print("[" + ''.join(["*"]) * n + ''.join([" "]) * (100 - n) + "]")
        f = float(r.collar_flux(t, rx, sx))  # exact root collar flux
        x_.append(t)
        y_.append(f)
        w_.append(water)
        cpx.append(rx[0])
        cps.append(float(sx[cci]))

        # print("Time:", t, ", collar flux", f, "cm^3/day at", rx[0], "cm xylem ", float(sx_old[cci]), "cm soil", "; domain water", s.getWaterVolume(), "cm3")

    t += dt

""" Plot """
if rank == 0:
    print ("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s")
    fig, ax1 = plt.subplots()
    ax1.plot(x_, trans * sinusoidal(x_), 'k')  # potential transpiration
    ax1.plot(x_, -np.array(y_), 'g')  # actual transpiration (neumann)
    ax2 = ax1.twinx()
    ax2.plot(x_, np.cumsum(-np.array(y_) * dt), 'c--')  # cumulative transpiration (neumann)
    ax1.set_xlabel("Time [d]")
    ax1.set_ylabel("Transpiration $[cm^3 d^{-1}]$")
    ax1.legend(['Potential', 'Actual', 'Cumulative'], loc = 'upper left')
    plt.show()

