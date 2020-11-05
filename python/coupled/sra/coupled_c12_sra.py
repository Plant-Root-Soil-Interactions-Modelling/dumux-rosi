import sys
sys.path.append("../../../build-cmake/rosi_benchmarking/python_solver/")
sys.path.append("../solvers/")  # for pure python solvers

from xylem_flux import XylemFluxPython  # Python hybrid solver
import plantbox as pb
import rsml_reader as rsml
from rosi_richards import RichardsSP  # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part
import vtk_plot as vp
import van_genuchten as vg
from root_conductivities import *

from math import *
import numpy as np
import matplotlib.pyplot as plt
import timeit
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()


def sinusoidal(t):
    return np.sin(2. * pi * np.array(t) - 0.5 * pi) + 1.


def mfp(h, soil):
#     return vg.matric_flux_potential(h, soil)
    return vg.fast_mfp[soil](h)


def imfp(mfp, soil):
#     return vg.matric_potential_mfp(h, soil)
    return vg.fast_imfp[soil](mfp)

""" 
Benchmark M1.2 static root system in soil (with the classic sink)

also works parallel with mpiexec (only slightly faster, due to overhead)
"""

""" Parameters """
min_b = [-4., -4., -15.]
max_b = [4., 4., 0.]
cell_number = [8, 8, 15]  # [8, 8, 15]  # [16, 16, 30]  # [32, 32, 60]  # [8, 8, 15]
periodic = False

name = "DuMux_1cm"
sand = [0.045, 0.43, 0.15, 3, 1000]
loam = [0.08, 0.43, 0.04, 1.6, 50]
clay = [0.1, 0.4, 0.01, 1.1, 10]
soil = loam

sp = vg.Parameters(soil)
vg.create_mfp_lookup(sp)
initial = -659.8 + 7.5  # -659.8

trans = 6.4  # cm3 /day (sinusoidal)
wilting_point = -15000  # cm

sim_time = 1  # [day] for task b
age_dependent = False  # conductivities
dt = 120. / (24 * 3600)  # [days] Time step must be very small

""" Initialize macroscopic soil model """
cpp_base = RichardsSP()
s = RichardsWrapper(cpp_base)
s.initialize()
s.createGrid(min_b, max_b, cell_number, periodic)  # [cm]
s.setHomogeneousIC(initial, True)  # cm pressure head, equilibrium
s.setTopBC("noFlux")
s.setBotBC("noFlux")
s.setVGParameters([soil])
s.initializeProblem()
s.setCriticalPressure(wilting_point)

""" Initialize xylem model (a) or (b)"""
r = XylemFluxPython("../grids/RootSystem8.rsml")
r.rs.setRectangularGrid(pb.Vector3d(min_b[0], min_b[1], min_b[2]), pb.Vector3d(max_b[0], max_b[1], max_b[2]),
                        pb.Vector3d(cell_number[0], cell_number[1], cell_number[2]), True)
init_conductivities(r, age_dependent)
r.rs.sort()  # ensures segment is located at index s.y-1
r.test()  # sanity checks
nodes = r.get_nodes()
rs_age = np.max(r.get_ages())

""" Coupling (map indices) """
picker = lambda x, y, z : s.pick([x, y, z])
r.rs.setSoilGrid(picker)  # maps segments
cci = picker(nodes[0, 0], nodes[0, 1], nodes[0, 2])  # collar cell index

""" MFP for Schroeder """
mfp_ = lambda h: mfp(h, sp)
imfp_ = lambda mfp: imfp(mfp, sp)

""" Numerical solution (a) """
start_time = timeit.default_timer()
x_, y_, w_, cpx, cps = [], [], [], [], []
sx = s.getSolutionHead()  # inital condition, solverbase.py
rsx = []

N = round(sim_time / dt)
t = 0.

for i in range(0, N):

    if rank == 0:  # Root simulation is not parallel

        if rsx:
            rsx = r.segSRA(rs_age + t, rx, sx, mfp_, imfp_)  # ! more unstable in case of mai scenario
            rx = r.solve(rs_age + t, -trans * sinusoidal(t), sx[cci], rsx, False, wilting_point, [])  # update rsx to last solution
        else:  # first call
            rx = r.solve(rs_age + t, -trans * sinusoidal(t), sx[cci], sx, True, wilting_point, [])  # this works

#         rsx = r.segSRA(rs_age + t, rx, sx, mfp_, imfp_)
#         seg_fluxes = r.segFluxes(rs_age + t, rx, rsx, approx = False, cells = False)
        seg_fluxes = r.segFluxes(rs_age + t, rx, sx, approx = False, cells = True)  # classic sink
        fluxes = r.sumSoilFluxes(seg_fluxes)

        sum_flux = 0.
        for f in fluxes.values():
            sum_flux += f
        print("Summed fluxes ", sum_flux, "= collar flux", r.collar_flux(rs_age + t, rx, rsx, [], cells = False), "= prescribed", -trans * sinusoidal(t))
        # print("Summed fluxes ", sum_flux, "= collar flux", r.collar_flux(rs_age + t, rx, sx), "= prescribed", -trans * sinusoidal(t))

    else:
        fluxes = None

    fluxes = comm.bcast(fluxes, root = 0)  # Soil part runs parallel
    s.setSource(fluxes)  # richards.py

    # simualte soil (parallel)
    s.ddt = dt / 10
    s.solve(dt)

    sx = s.getSolutionHead()  # richards.py
    water = s.getWaterVolume()

    if rank == 0:
        n = round(float(i) / float(N) * 100.)
        min_sx = np.min(sx)
        min_rx = np.min(rx)
        max_sx = np.max(sx)
        max_rx = np.max(rx)
        print("[" + ''.join(["*"]) * n + ''.join([" "]) * (100 - n) + "], [{:g}, {:g}] cm soil [{:g}, {:g}] cm root at {:g} days {:g}"
              .format(min_sx, max_sx, min_rx, max_rx, s.simTime, rx[0]))
        # f = float(r.collar_flux(rs_age + t, rx, sx))  # exact root collar flux
        f = r.collar_flux(rs_age + t, rx, rsx, [], cells = False),
        x_.append(t)
        y_.append(sum_flux)
        w_.append(water)
        cpx.append(rx[0])
        cps.append(float(sx[cci]))

        # print("Time:", t, ", collar flux", f, "cm^3/day at", rx[0], "cm xylem ", float(sx_old[cci]), "cm soil", "; domain water", s.getWaterVolume(), "cm3")

    t += dt

s.writeDumuxVTK(name)

""" Plot """
if rank == 0:
    print ("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s")

    # VTK vizualisation
    vp.plot_roots_and_soil(r.rs, "pressure head", rx, s, periodic, min_b, max_b, cell_number, name)

    fig, ax1 = plt.subplots()
    ax1.plot(x_, trans * sinusoidal(x_), 'k')  # potential transpiration
    ax1.plot(x_, -np.array(y_), 'g')  # actual transpiration (neumann)
    ax2 = ax1.twinx()
    ax2.plot(x_, np.cumsum(-np.array(y_) * dt), 'c--')  # cumulative transpiration (neumann)
    ax1.set_xlabel("Time [d]")
    ax1.set_ylabel("Transpiration $[cm^3 d^{-1}]$")
    ax1.legend(['Potential', 'Actual', 'Cumulative'], loc = 'upper left')
    np.savetxt(name, np.vstack((x_, -np.array(y_))), delimiter = ';')
    plt.show()

