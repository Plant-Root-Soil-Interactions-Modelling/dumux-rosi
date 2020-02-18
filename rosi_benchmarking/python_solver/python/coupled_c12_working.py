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

""" 
Benchmark M1.2 static root system in soil
(does not work in parallel)
"""


def sinusoidal(t):
    return np.sin(2. * pi * np.array(t) - 0.5 * pi) + 1.


""" Parameters """
kz = 4.32e-2  # [cm^3/day]
kr = 1.728e-4  # [1/day]
trans = 6.4  # cm3 /day
p0 = -1000  # dircichlet bc at top
sim_time = 7  # [day] for task b
wilting_point = -10000  # cm
loam = [0.08, 0.43, 0.04, 1.6, 50]

# sp = vg.Parameters(loam)
# wp_wc = vg.water_content(wilting_point, sp)
# print("Water content at wilting point is ", wp_wc)
# cell_volume = 1  # cm3

""" root problem """
r = XylemFluxPython("../grids/RootSystem.rsml")
r.setKr([kr])
r.setKx([kz])
nodes = r.get_nodes()

""" soil problem """
cpp_base = RichardsSP()
s = RichardsWrapper(cpp_base)
s.initialize()
s.createGrid([-4., -4., -20.], [4., 4., 0.], [8, 8, 20])  # [cm]
s.setHomogeneousIC(-659.8 - 7.5, True)  # cm pressure head, equilibrium
s.setTopBC("noFlux")
s.setBotBC("noFlux")
s.setVGParameters([loam])
s.initializeProblem()

""" coupling """
picker = lambda x, y, z : s.pick(x, y, z)  # mpi
r.rs.setSoilGrid(picker)
cci = picker(nodes[0, 0], nodes[0, 1], nodes[0, 2])  # collar cell index

""" Numerical solution (a) """
fig, (ax1, ax2) = plt.subplots(1, 2)

sx = s.toHead(s.getSolution())  # inital condition
rx = r.solve(0., 0., True)
# print("collar cell index", cci)
# print("collar soil pressure", sx[cci], "cm")

dt = 60. / (24 * 3600)  # [days] Time step must be very small
t = 0.

start_time = timeit.default_timer()
x_, y_, w_, cpx, cps = [], [], [], [], []
for i in range(0, sim_time * round(1. / dt)):

    rx_hom = r.solve(t, -trans * sinusoidal(t), sx[cci], wilting_point)  # xylem_flux.py
    # rx_hom = r.solve_dirichlet(t, -10000, sx[cci]) # xylem_flux.py
    # rx_hom = r.solve_neumann(t, -trans * sinusoidal(t)) # xylem_flux.py

    rx = r.getSolution(rx_hom, sx)  # class XylemFlux is defined in MappedOrganism.h

    fluxes = r.soilFluxes(t, rx_hom)
#     print("Summed fluxes of root to soil ", np.sum(np.array(list(fluxes.values()))), "cm3 day-1")
#     for cellIdx, flux in fluxes.items():
#         fluxes[cellIdx] = min(fluxes[cellIdx], 0.)

    s.setSource(fluxes)  # richards.py
    s.solve(dt)
    sx = s.toHead(s.getSolution())

    f = float(r.collar_flux(t, rx, sx))  # exact root collar flux
    x_.append(t)
    y_.append(f)
    w_.append(s.getWaterVolume())
    cpx.append(rx[0])
    cps.append(float(sx[cci]))
    # print("Time:", t, ", collar flux", f, "cm^3/day at", rx[0], "cm xylem ", float(sx_old[cci]), "cm soil", "; domain water", s.getWaterVolume(), "cm3")

    t += dt

# print ("Coupled example solved in", timeit.default_timer() - start, " s")

ax1.plot(x_, trans * sinusoidal(x_), "k")
ax1.plot(x_, -np.array(y_), "r")
ax1.set_xlabel("Time (days)")
ax1.set_ylabel("Transpiration (cm$^3$ day$^{-1})$")

# ax2.plot(x_, cpx, "g")
# ax2.plot(x_, cps, "b")
# ax2.set_xlabel("Time (days)")
# ax2.set_ylabel("Water potential (cm)")
# ax2.legend(["xylem", "soil"])

ax2.plot(x_, w_, "b")
ax2.set_xlabel("Time (days)")
ax2.set_ylabel("Water volume cm$^3$$")

# ax1.plot(rx, nodes[:, 2] , "r*")
# ax1.set_xlabel("Xylem pressure (cm)")
# ax1.set_ylabel("Depth (m)")
# ax1.set_title("Constant conductivities")

# print()
#
# """ Numerical solution (b) """
#
# kr0 = np.array([[0, 1.14e-03], [2, 1.09e-03], [4, 1.03e-03], [6, 9.83e-04], [8, 9.35e-04], [10, 8.90e-04], [12, 8.47e-04], [14, 8.06e-04], [16, 7.67e-04], [18, 7.30e-04], [20, 6.95e-04], [22, 6.62e-04], [24, 6.30e-04], [26, 5.99e-04], [28, 5.70e-04], [30, 5.43e-04], [32, 5.17e-04]])
# kr1 = np.array([[0, 4.11e-03], [1, 3.89e-03], [2, 3.67e-03], [3, 3.47e-03], [4, 3.28e-03], [5, 3.10e-03], [6, 2.93e-03], [7, 2.77e-03], [8, 2.62e-03], [9, 2.48e-03], [10, 2.34e-03], [11, 2.21e-03], [12, 2.09e-03], [13, 1.98e-03], [14, 1.87e-03], [15, 1.77e-03], [16, 1.67e-03], [17, 1.58e-03]])
# r.setKrTables([kr0[:, 1], kr1[:, 1], kr1[:, 1], kr1[:, 1]], [kr0[:, 0], kr1[:, 0], kr1[:, 0], kr1[:, 0]])
#
# kx0 = np.array([[0, 6.74e-02], [2, 7.48e-02], [4, 8.30e-02], [6, 9.21e-02], [8, 1.02e-01], [10, 1.13e-01], [12, 1.26e-01], [14, 1.40e-01], [16, 1.55e-01], [18, 1.72e-01], [20, 1.91e-01], [22, 2.12e-01], [24, 2.35e-01], [26, 2.61e-01], [28, 2.90e-01], [30, 3.21e-01], [32, 3.57e-01]])
# kx1 = np.array([[0, 4.07e-04], [1, 5.00e-04], [2, 6.15e-04], [3, 7.56e-04], [4, 9.30e-04], [5, 1.14e-03], [6, 1.41e-03], [7, 1.73e-03], [8, 2.12e-03], [9, 2.61e-03], [10, 3.21e-03], [11, 3.95e-03], [12, 4.86e-03], [13, 5.97e-03], [14, 7.34e-03], [15, 9.03e-03], [16, 1.11e-02], [17, 1.36e-02]])
# r.setKxTables([kx0[:, 1], kx1[:, 1], kx1[:, 1], kx1[:, 1]], [kx0[:, 0], kx1[:, 0], kx1[:, 0], kx1[:, 0]])
#
# rx_hom = r.solve(simtime, p0 - p_s, False)
# rx = r.getSolution(rx_hom, [p_s])
#
# flux = r.collar_flux(simtime, rx, [p_s])
# print("Transpiration", flux, "cm3/day")
#
# ax2.plot(rx, nodes[:, 2] , "r*")
# ax2.set_xlabel("Xylem pressure (cm)")
# ax2.set_ylabel("Depth (m)")
# ax2.set_title("Age dependent conductivities")

plt.show()
