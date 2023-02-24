""" 
Jan's new scenarios with the new SRA sink

2d soil, 2*2*2*cm^3 voxels, dynamic conductivities (see root_conductivities.py)
"""

import sys; sys.path.append("../../../modules/"); sys.path.append("../../../../../CPlantBox/");  sys.path.append("../../../../../CPlantBox/src/python_modules")
sys.path.append("../../../../build-cmake/cpp/python_binding/"); sys.path.append("../../../modules/fv/");
sys.path.append("../");

from xylem_flux import XylemFluxPython  # Python hybrid solver
import plantbox as pb
import rsml_reader as rsml
from rosi_richards import RichardsSP  # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part
import vtk_plot as vp
import van_genuchten as vg
from root_conductivities import *
from rhizo_models import plot_transpiration

from math import *
import numpy as np
import matplotlib.pyplot as plt
import timeit
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()
from scipy.optimize import fsolve
from scipy import sparse
import scipy.sparse.linalg as LA

from sra_table_lookup import *  # <------- new


def sinusoidal(t):
    return np.sin(2. * pi * np.array(t) - 0.5 * pi) + 1.


def soil_root_interface_table2(rx, sx, inner_kr_, rho_, f):
    assert rx.shape == sx.shape
    rsx = f((rx, sx, inner_kr_ , rho_))
    return rsx

""" 
Parameters  
"""

""" soil """
min_b = [-7.5, -37.5, -110.]
max_b = [7.5, 37.5, 0.]
cell_number = [8, 38, 55]  # 2*2*2 cm3
periodic = True  # check data first
fname = "../../../../grids/RootSystem_verysimple2.rsml"

name = "dry_3d_2"  # name to export resutls

alpha = 0.018;  # (cm-1)
n = 1.8;
Ks = 28.46;  # (cm d-1)
loam = [0.08, 0.43, alpha, n, Ks]
p_top = -5000  # -5000 (dry), -310 (wet)
p_bot = -200  #

soil_ = loam
soil = vg.Parameters(soil_)

""" root system """
trans = 0.5 * 15 * 75  # average per day [cm3 /day] (sinusoidal)
wilting_point = -15000  # [cm]
predefined_growth = False  # root growth by setting radial conductivities
rs_age = 78  # initial root system age

""" simulation time """
sim_time = 7.1  # 0.65  # 0.25  # [day]
dt = 60 / (24 * 3600)  # time step [day], 120 schwankt stark
NT = int(np.ceil(sim_time / dt))  # number of iterations
skip = 20  # for output and results, skip iteration

""" Initialize macroscopic soil model """
s = RichardsWrapper(RichardsSP())
s.initialize()
s.createGrid(min_b, max_b, cell_number, periodic)  # [cm]
s.setLinearIC(p_top, p_bot)  # cm pressure head, equilibrium
s.setTopBC("noFlux")
s.setBotBC("noFlux")
s.setVGParameters([soil_])
# s.setParameter("Newton.EnableChop", "True")
s.setParameter("Newton.EnableAbsoluteResidualCriterion", "True")
s.setParameter("Soil.SourceSlope", "1000")  # turns regularisation of the source term on
s.initializeProblem()
s.setCriticalPressure(wilting_point)

""" Initialize xylem model (a) or (b)"""
r = XylemFluxPython(fname)
r.rs.setRectangularGrid(pb.Vector3d(min_b[0], min_b[1], min_b[2]), pb.Vector3d(max_b[0], max_b[1], max_b[2]),
                        pb.Vector3d(cell_number[0], cell_number[1], cell_number[2]), True)  # cutting
init_conductivities_scenario_jan(r)

""" Coupling (map indices) """
picker = lambda x, y, z: s.pick([x, y, z])
cci = picker(r.rs.nodes[0].x, r.rs.nodes[0].y, r.rs.nodes[0].z)  # collar cell index
r.rs.setSoilGrid(picker)  # maps segment
outer_r = r.rs.segOuterRadii()

""" sanity checks """
# r.plot_conductivities()

types = r.rs.subTypes
types = (np.array(types) == 12) * 1  # all 0, only 12 are laterals
r.rs.subTypes = list(types)
r.test()  # sanity checks

seg_length = r.rs.segLength()
ns = len(seg_length)
# print("press any key"); input()
print("outer radii", np.min(outer_r) , np.max(outer_r))
radii = r.rs.radii
print("inner radii", np.min(radii) , np.max(radii))

rho_ = np.divide(outer_r, np.array(radii))
print("rho", np.min(rho_) , np.max(rho_), np.sum(rho_))
rho_ = np.maximum(np.minimum(rho_, 200 * np.ones(rho_.shape)), 10 * np.ones(rho_.shape))

radius_min, radius_max = np.min(radii) , np.max(radii)
kr_max = 0.000181
kr_min = 0.0000173
print("inner_r*kr", kr_min * radius_min, kr_max * radius_max)

sra_table_lookup = open_sra_lookup("../table_jan2")  # table_jan ()

# quick check
# rsx2 = soil_root_inerface(np.array([-15000]), np.array([-700]), r, sp, outer_r)
# print(r.rs.radii[0])
# print(outer_r[0])
# rsx3 = sra_table_lookup((-15000, -700, 0.1679960056208074, 0.6952332821448589))
# print(rsx2, rsx3)

""" for fixed root system """
inner_ = np.zeros((len(outer_r),))
outer_ = np.zeros(len(outer_r),)
for i in range(0, len(outer_r)):
    inner_[i] = max(min(radii[i] , 0.2), 0.01)
    outer_[i] = max(min(outer_r[i] , 20), 0.1)
mapping = np.array([r.rs.seg2cell[j] for j in range(0, ns)])

""" Numerical solution (a) """
start_time = timeit.default_timer()
x_, y_, w_, cpx, cps, cf = [], [], [], [], [], []
sink1d = []
sink1d2 = []
out_times = []  # days
sx = s.getSolutionHead()  # inital condition, solverbase.py
hsb = np.array([sx[mapping[j]][0] for j in range(0, ns)])  # soil bulk matric potential per segment
kr_ = np.zeros((ns,))
rsx = hsb.copy()  # initial values for fix point iteration

t = 0.

for i in range(0, NT):

    wall_iteration = timeit.default_timer()

    if rank == 0:  # root part is not parallel

        wall_fixpoint = timeit.default_timer()

        if i == 0:  # only first time
            rx = r.solve(rs_age + t, -trans * sinusoidal(t), 0., rsx, False, wilting_point, [])  # xylem_flux.py, cells = False
            rx_old = rx.copy()

        err = 1.e6
        c = 1
        while err > 1. and c < 100:

            """ interpolation """
            wall_interpolation = timeit.default_timer()
            for j in range(0, len(outer_r)):  # determine kr at this time step
                kr_[j] = r.kr_f(rs_age + t, types[j])
            inner_kr_ = np.multiply(radii, kr_)  # multiply for table look up
            rsx = soil_root_interface_table2(rx[1:], hsb, inner_kr_, rho_, sra_table_lookup)
            wall_interpolation = timeit.default_timer() - wall_interpolation

            """ xylem matric potential """
            wall_xylem = timeit.default_timer()
            rx = r.solve(rs_age + t, -trans * sinusoidal(t), 0., rsx, False, wilting_point, [])  # xylem_flux.py, cells = False
            err = np.linalg.norm(rx - rx_old)
            wall_xylem = timeit.default_timer() - wall_xylem
            # print(err)
            rx_old = rx.copy()
            c += 1

        # print(c, ": ", np.sum(rx[1:]), np.sum(hsb), np.sum(inner_kr_), np.sum(rho_))
        # print(c, "iterations", wall_interpolation / (wall_interpolation + wall_xylem), wall_xylem / (wall_interpolation + wall_xylem))

        wall_fixpoint = timeit.default_timer() - wall_fixpoint

        # rx = r.solve(rs_age + t, -trans * sinusoidal(t), 0., rsx, False, wilting_point, [])  # xylem_flux.py, cells = False
        fluxes = np.array(r.segFluxes(rs_age + t, rx, rsx, False))  # class XylemFlux is defined in MappedOrganism.h, approx = True
        min_rsx = np.min(rsx)  # for console output
        max_rsx = np.max(rsx)

#         sum_flux = 0.
#         for f in fluxes.values():
#             sum_flux += f
#         print(sum_flux, r.collar_flux(rs_age + t, rx, rsx))

    else:
        fluxes = None

    wall_soil = timeit.default_timer()
    fluxes = comm.bcast(r.sumSegFluxes(fluxes), root = 0)  # Soil part runs parallel
    s.setSource(fluxes.copy())  # richards.py
    s.solve(dt)
    sx = s.getSolutionHead()  # richards.py
    hsb = np.array([sx[mapping[j]][0] for j in range(0, ns)])  # soil bulk matric potential per segment
    water = s.getWaterVolume()
    wall_soil = timeit.default_timer() - wall_soil

    wall_iteration = timeit.default_timer() - wall_iteration

    if rank == 0 and i % skip == 0:
        min_sx = np.min(sx)
        max_sx = np.max(sx)
        min_rx = np.min(rx)
        max_rx = np.max(rx)
        x_.append(t)
        sum_flux = 0.
        for f in fluxes.values():
            sum_flux += f
        cf_ = r.collar_flux(rs_age + t, rx, rsx, k_soil = [], cells = False)
        print("Summed fluxes ", sum_flux, "= collar flux", cf_, "= prescribed", -trans * sinusoidal(t))
        y_.append(sum_flux)  # cm3/day
        w_.append(water)  # cm3
        cf.append(cf_)  # cm3/day
        cpx.append(rx[0])  # cm
        cps.append(float(sx[cci]))  # cm
        n = round(float(i) / float(NT) * 100.)
        print("[" + ''.join(["*"]) * n + ''.join([" "]) * (100 - n) + "], soil sx [{:g}, {:g}], interface [{:g}, {:g}] cm, root [{:g}, {:g}] cm, {:g} days"
              .format(min_sx, max_sx, min_rsx, max_rsx, min_rx, max_rx, s.simTime))
        print("Iteration {:g} took {:g} seconds [{:g} fixpoint iteration, {:g} soil] \n".format(i, wall_iteration, wall_fixpoint, wall_soil))

        """ Additional sink plot """
        if i % (60 * 6) == 0:  # every 6h
            ana = pb.SegmentAnalyser(r.rs)
            fluxes = r.segFluxes(rs_age + t, rx, rsx, False)
            ana.addData("fluxes", fluxes)  # cut off for vizualisation
            flux1d = ana.distribution("fluxes", max_b[2], min_b[2], cell_number[2], True)
            sink1d.append(np.array(flux1d))

    t += dt

s.writeDumuxVTK(name)

""" Plot """
if rank == 0:
    print ("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s")

    np.savetxt(name, np.vstack((x_, -np.array(y_))), delimiter = ';')
    sink1d = np.array(sink1d)
    np.save(name + "_sink", sink1d)

#     plot_transpiration(x_, y_, cf, lambda t: trans * sinusoidal(t))
#
#     ana = pb.SegmentAnalyser(r.rs)
#     ana.addData("pressure", rx)
#     vp.plot_roots(ana, "pressure")

    vp.plot_roots_and_soil(r.rs, "pressure head", rx, s, periodic, min_b, max_b, cell_number, name)  # VTK vizualisation
