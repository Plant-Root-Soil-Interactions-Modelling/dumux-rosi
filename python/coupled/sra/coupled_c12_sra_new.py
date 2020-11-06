import sys; sys.path.append("../../modules/"); sys.path.append("../../../../CPlantBox/");  sys.path.append("../../../build-cmake/cpp/python_binding/")
sys.path.append("../")

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

name = "DuMux_1cm_schroeder_clay"
sand = [0.045, 0.43, 0.15, 3, 1000]
loam = [0.08, 0.43, 0.04, 1.6, 50]
clay = [0.1, 0.4, 0.01, 1.1, 10]
soil = clay

sp = vg.Parameters(soil)
vg.create_mfp_lookup(sp)
initial = -659.8 + 7.5  # -659.8 + 7.5  # -659.8

trans = 6.4  # cm3 /day (sinusoidal)
wilting_point = -15000  # cm

sim_time = 1  # [day] for task b
age_dependent = False  # conductivities
dt = 3600. / (24 * 3600)  # [days] Time step must be very small
dx = 1.e-2

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
# s.setRegularisation(1.e-4, 1.e-4)

""" Initialize xylem model (a) or (b)"""
r = XylemFluxPython("../../../grids/RootSystem8.rsml")
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
water_uptake = []

N = round(sim_time / dt)
t = 0.

k = vg.hydraulic_conductivity(wilting_point, sp)  # hydraulic conductivity at wilting point
rx = []

for i in range(0, N):

    sx = s.getSolutionHead()  # inital condition, solverbase.py

    if rank == 0:  # Root simulation is not parallel

#         if len(rx) > 0:
#             rsx = r.segSchroeder(rs_age + t, rx, sx, wilting_point, mfp_, imfp_)  # ! more unstable in case of mai scenario
#             rx = r.solve(rs_age + t, -trans * sinusoidal(t), sx[cci], rsx, False, wilting_point, [])  # update rsx to last solution
#         else:  # first time step

        rx = r.solve(rs_age + t, -trans * sinusoidal(t), sx[cci], sx, True, wilting_point, [])  # [cm] in xylem_flux.py, cells = True
        seg_nostress = np.array(r.segFluxes(rs_age + t, rx, sx, approx = False, cells = True))  # classic sink in case of no stress

        seg_stress = np.array(r.segSRAStressedFlux(sx, wilting_point, k, mfp_, imfp_, dx))  # steady rate approximation in case of stress
        print("stressed:", np.min(seg_stress), np.max(seg_stress), np.sum(seg_stress))
        print("nostress:", np.min(seg_nostress), np.max(seg_nostress), np.sum(seg_nostress), "at", -trans * sinusoidal(t))

        # seg_stress = np.minimum(np.zeros(seg_stress.shape), seg_stress)  # use only for inflow
        seg_stress = np.maximum(seg_nostress, seg_stress)  # limit by potential transpiration, ensure unstressed>stressed
        # seg_stress = np.minimum(np.zeros(seg_stress.shape), seg_stress)  # use only for inflow

        seg_head = np.array(r.segSRA(rs_age + t, rx, sx, mfp_, imfp_))  # to determine if stressed or not
        seg_fluxes = np.zeros(seg_nostress.shape)
        ii = seg_head > -1  # indices of stressed segments
        seg_fluxes[ii] = seg_stress[ii]
        seg_fluxes[~ii] = seg_nostress[~ii]  # ~ = boolean not
#         # loop version (of above)
#         for j in range(0, seg_fluxes.shape[0]):
#             seg_fluxes[j] = (seg_head[j] > -1) * seg_stress[j] + (seg_head[j] <= -1) * seg_nostress[j]

        fluxes = r.sumSoilFluxes(seg_fluxes)  # seg_fluxes

        sum_flux = 0.
        for f in fluxes.values():
            sum_flux += f
        water_uptake.append(sum_flux)
        print("Summed fluxes ", sum_flux, "= collar flux", r.collar_flux(rs_age + t, rx, sx, [], cells = True), "= prescribed", -trans * sinusoidal(t))

    else:
        fluxes = None

    fluxes = comm.bcast(fluxes, root = 0)  # Soil part runs parallel

    s.setSource(fluxes)  # richards.py
#     sx = s.applySource(dt, sx, fluxes, wilting_point)
#     s.setInitialCondition(sx)

    # simualte soil (parallel)

    try:
        s.ddt = dt
        s.solve(dt)
    except:
        print("TRYING AGAIN")
        s.setInitialCondition(sx)
        s.ddt = dt / 100
        s.solve(dt / 2)
        s.solve(dt / 2)

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
        x_.append(t)
        y_.append(sum_flux)  # sum_flux
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

    np.savetxt(name, np.vstack((x_, -np.array(y_), -np.array(water_uptake))), delimiter = ';')

    plt.show()

