import sys
sys.path.append("../../../build-cmake/rosi_benchmarking/python_solver/")

from solver.xylem_flux import XylemFluxPython  # Python hybrid solver
import solver.plantbox as pb
import solver.rsml_reader as rsml
from rosi_richards import RichardsSP  # C++ part (Dumux binding)
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

also works parallel with mpiexec (only slightly faster, due to overhead)
"""

""" Parameters """
sim_time = 0.5  # [day] for task b
trans = 6.4  # cm3 /day (sinusoidal)
wilting_point = -10000  # cm
loam = [0.08, 0.43, 0.04, 1.6, 50]
periodic = False  # for perioidc, we take

""" Root problem (a) or (b)"""
r = XylemFluxPython("../grids/RootSystem.rsml")
# print("Types:")
# print(r.rs.types)

# (a)
r.setKr([ 1.728e-4])  # [cm^3/day]
r.setKx([4.32e-2 ])  # [1/day]

# (b)
# kr0 = np.array([[0, 1.14e-03], [2, 1.09e-03], [4, 1.03e-03], [6, 9.83e-04], [8, 9.35e-04], [10, 8.90e-04], [12, 8.47e-04], [14, 8.06e-04], [16, 7.67e-04], [18, 7.30e-04], [20, 6.95e-04], [22, 6.62e-04], [24, 6.30e-04], [26, 5.99e-04], [28, 5.70e-04], [30, 5.43e-04], [32, 5.17e-04]])
# kr1 = np.array([[0, 4.11e-03], [1, 3.89e-03], [2, 3.67e-03], [3, 3.47e-03], [4, 3.28e-03], [5, 3.10e-03], [6, 2.93e-03], [7, 2.77e-03], [8, 2.62e-03], [9, 2.48e-03], [10, 2.34e-03], [11, 2.21e-03], [12, 2.09e-03], [13, 1.98e-03], [14, 1.87e-03], [15, 1.77e-03], [16, 1.67e-03], [17, 1.58e-03]])
# r.setKrTables([kr0[:, 1], kr1[:, 1], kr1[:, 1], kr1[:, 1]], [kr0[:, 0], kr1[:, 0], kr1[:, 0], kr1[:, 0]])
# kx0 = np.array([[0, 6.74e-02], [2, 7.48e-02], [4, 8.30e-02], [6, 9.21e-02], [8, 1.02e-01], [10, 1.13e-01], [12, 1.26e-01], [14, 1.40e-01], [16, 1.55e-01], [18, 1.72e-01], [20, 1.91e-01], [22, 2.12e-01], [24, 2.35e-01], [26, 2.61e-01], [28, 2.90e-01], [30, 3.21e-01], [32, 3.57e-01]])
# kx1 = np.array([[0, 4.07e-04], [1, 5.00e-04], [2, 6.15e-04], [3, 7.56e-04], [4, 9.30e-04], [5, 1.14e-03], [6, 1.41e-03], [7, 1.73e-03], [8, 2.12e-03], [9, 2.61e-03], [10, 3.21e-03], [11, 3.95e-03], [12, 4.86e-03], [13, 5.97e-03], [14, 7.34e-03], [15, 9.03e-03], [16, 1.11e-02], [17, 1.36e-02]])
# r.setKxTables([kx0[:, 1], kx1[:, 1], kx1[:, 1], kx1[:, 1]], [kx0[:, 0], kx1[:, 0], kx1[:, 0], kx1[:, 0]])

nodes = r.get_nodes()

""" Soil problem """
cpp_base = RichardsSP()
s = RichardsWrapper(cpp_base)
s.initialize()

s.createGrid([-4., -4., -20.], [4., 4., 0.], [16, 16, 40])  # [cm]
r.rs.setRectangularGrid(pb.Vector3d(-4., -4., -20.), pb.Vector3d(4., 4., 0.), pb.Vector3d(16, 16, 40))  # cut root segments to grid (segments are not mapped after)

# s.createGrid([-4., -4., -20.], [4., 4., 0.], [8, 8, 20])  # [cm]
# r.rs.setRectangularGrid(pb.Vector3d(-4., -4., -20.), pb.Vector3d(4., 4., 0.), pb.Vector3d(8, 8, 20))  # cut root segments to grid (segments are not mapped after)

s.setHomogeneousIC(-669.8 - 10, True)  # cm pressure head, equilibrium
s.setTopBC("noFlux")
s.setBotBC("noFlux")
s.setVGParameters([loam])
s.initializeProblem()
s.setCriticalPressure(wilting_point)
s.setRegularisation(1.e-4, 1.e-4)

""" Coupling (map indices) """
picker = lambda x, y, z : s.pick([x, y, z])
r.rs.setSoilGrid(picker)  # maps segments
cci = picker(nodes[0, 0], nodes[0, 1], nodes[0, 2])  # collar cell index

""" Numerical solution (a) """
start_time = timeit.default_timer()
x_, y_, w_, cpx, cps = [], [], [], [], []
sx = s.getSolutionHead()  # inital condition, solverbase.py

dt = 120. / (24 * 3600)  # [days] Time step must be very small
N = round(sim_time / dt)
t = 0.

for i in range(0, N):

    if rank == 0:  # Root part is not parallel
        rx = r.solve(t, -trans * sinusoidal(t), sx[cci], sx, True, wilting_point)  # xylem_flux.py
        fluxes = r.soilFluxes(t, rx, sx, approx=False)  # class XylemFlux is defined in MappedOrganism.h
        sum_flux = 0.
        for f in fluxes.values():
            sum_flux += f
        print("Fuxes ", sum_flux, "= prescribed", -trans * sinusoidal(t) , "= collar flux", r.collar_flux(0., rx, sx))
    else:
        fluxes = None

    fluxes = comm.bcast(fluxes, root=0)  # Soil part runs parallel
    # s.setSource(fluxes)  # richards.py
    sx = s.getSolutionHead()  # richards.py
    sx = s.applySource(dt, sx, fluxes, wilting_point)  # richards.py
    s.setInitialCondition(sx)    

    s.solve(dt)
    
    water = s.getWaterVolume()
    print("Minimum", np.min(sx), "cm", "Total  water: ", water, "cm3")

    if rank == 0:
        n = round(float(i) / float(N) * 100.)
        min_sx = np.min(sx)
        min_rx = np.min(rx)
        max_sx = np.max(sx)
        max_rx = np.max(rx)
        print("[" + ''.join(["*"]) * n + ''.join([" "]) * (100 - n) + "], [{:g}, {:g}] cm soil [{:g}, {:g}] cm root at {:g} days {:g}"
              .format(min_sx, max_sx, min_rx, max_rx, s.simTime, rx[0]))
        f = float(r.collar_flux(t, rx, sx))  # exact root collar flux
        x_.append(t)
        y_.append(f)
        w_.append(water)  # diff(watter)/dt/1000.
        cpx.append(rx[0])
        cps.append(float(sx[cci]))

        # print("Time:", t, ", collar flux", f, "cm^3/day at", rx[0], "cm xylem ", float(sx_old[cci]), "cm soil", "; domain water", s.getWaterVolume(), "cm3")

    t += dt

s.writeDumuxVTK("c12_final")

""" Plot """
if rank == 0:
    print ("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s")
    fig, ax1 = plt.subplots()
    ax1.plot(x_, trans * sinusoidal(x_), 'k')  # potential transpiration
    ax1.plot(x_, -np.array(y_), 'g')  # actual transpiration (neumann)
    ax2 = ax1.twinx()
     # ax2.plot(x_, np.cumsum(-np.array(y_) * dt), 'c--')  # cumulative transpiration (neumann)
    ax2.plot(x_[1:], np.diff(np.array(w_)) / (8 * 8 * 20) / dt, 'c--')  # cumulative transpiration (neumann)
    ax1.set_xlabel("Time [d]")
    ax1.set_ylabel("Transpiration $[cm^3 d^{-1}]$")
    ax1.legend(['Potential', 'Actual', 'Cumulative'], loc='upper left')
    plt.show()

