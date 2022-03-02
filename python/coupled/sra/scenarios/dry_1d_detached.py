""" 
Jan's new scenarios with the new SRA sink

see aggregated rs
"""

import sys; sys.path.append("../../../modules/"); sys.path.append("../../../../../CPlantBox/");  sys.path.append("../../../../../CPlantBox/src/python_modules")
sys.path.append("../../../../build-cmake/cpp/python_binding/"); sys.path.append("../../../modules/fv/");
sys.path.append("../");

import plantbox as pb  # CPlantBox
from rosi_richards import RichardsSP  # C++ part (Dumux binding), macroscopic soil model
from richards import RichardsWrapper  # Python part, macroscopic soil model
from xylem_flux import *  # root system Python hybrid solver
from rhizo_models import *  # Helper class for cylindrical rhizosphere models
import aggregated_rs as agg

import vtk_plot as vp
import van_genuchten as vg
from sra_table_lookup import *

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import fsolve


def soil_root_interface(rx, sx, inner_kr, rho, sp):
    """
    rx             xylem matric potential [cm]
    sx             bulk soil matric potential [cm]
    inner_kr       root radius times hydraulic conductivity [cm/day] 
    rho            geometry factor [1]
    sp             soil van Genuchten parameters (type vg.Parameters)
    """
    k_soilfun = lambda hsoil, hint: (vg.fast_mfp[sp](hsoil) - vg.fast_mfp[sp](hint)) / (hsoil - hint)
    # rho = outer_r / inner_r  # Eqn [5]
    rho2 = rho * rho  # rho squared
    # b = 2 * (rho2 - 1) / (1 + 2 * rho2 * (np.log(rho) - 0.5))  # Eqn [4]
    b = 2 * (rho2 - 1) / (1 - 0.53 * 0.53 * rho2 + 2 * rho2 * (np.log(rho) + np.log(0.53)))  # Eqn [7]
    fun = lambda x: (inner_kr * rx + b * sx * k_soilfun(sx, x)) / (b * k_soilfun(sx, x) + inner_kr) - x
    rsx = fsolve(fun, rx)
    return rsx


def soil_root_interface_table2(rx, sx, inner_kr_, rho_, f):
    assert rx.shape == sx.shape
    try:
        rsx = f((rx, sx, inner_kr_ , rho_))
    except:
        print("failed:", rx, sx, inner_kr_ , rho_)
    return rsx


def double_(rsx):
    rsx2 = np.array([ 0. if i % 2 == 0 else rsx[int(i / 2)] for i in range(0, ns)])
    return rsx2

""" 
Parameters  
"""

""" soil """
name = "dry_agg"  # name to export resutls
min_b = [-7.5, -37.5, -110.]
max_b = [7.5, 37.5, 0.]
cell_number = [1, 1, 55]  # [8, 38, 55]  # 2cm3
periodic = True  # check data first
fname = "../../../../grids/RootSystem_verysimple2.rsml"
alpha = 0.018  # (cm-1)
n = 1.8
Ks = 28.46  # (cm d-1)
loam = [0.08, 0.43, alpha, n, Ks]
p_top = -5000  #  (dry), -310 (wet)
p_bot = -200  #
soil_ = loam
soil = vg.Parameters(soil_)

""" root system """
trans = 0.5 * 15 * 75  # average per day [cm3 /day] (sinusoidal)
wilting_point = -15000  # [cm]
predefined_growth = False  # root growth by setting radial conductivities
rs_age = 78  # initial root system age

""" simulation time """
sim_time = 0.1  # 0.65  # 0.25  # [day]
dt = 60 / (24 * 3600)  # time step [day], 120 schwankt stark
NT = int(np.ceil(sim_time / dt))  # number of iterations
skip = 20  # for output and results, skip iteration

""" Initialize macroscopic soil model """
s = RichardsWrapper(RichardsSP())
s.initialize()
s.createGrid(min_b, max_b, cell_number, periodic)  # [cm]
s.setLinearIC(p_top, p_bot)  # cm pressure head
s.setTopBC("noFlux")
s.setBotBC("noFlux")
s.setVGParameters([soil_])
s.setParameter("Newton.EnableAbsoluteResidualCriterion", "True")
s.setParameter("Soil.SourceSlope", "1000")  # turns regularisation of the source term on
s.initializeProblem()
s.setCriticalPressure(wilting_point)
s.ddt = 1.e-5  # [day] initial Dumux time step

""" 
Initialize xylem model 
"""
rr = XylemFluxPython(fname)
types = rr.rs.subTypes  # simplify root types
types = (np.array(types) >= 12) * 1  # all roots type 0, only >=12 are laterals type 1
rr.rs.subTypes = list(types)
agg.init_conductivities(rr)
r = agg.create_aggregated_rs(rr, rs_age, min_b, max_b, cell_number)

""" Coupling (map indices) """
picker = lambda x, y, z: s.pick([x, y, z])
r.rs.setSoilGrid(picker)  # maps segment
outer_r = r.rs.segOuterRadii()

""" sanity checks """
# r.plot_conductivities()
types = r.rs.subTypes
types = (np.array(types) == 12) * 1  # all 0, only 12 are laterals
r.rs.subTypes = list(types)
# r.test()  # sanity checks

seg_length = r.segLength()
ns = len(seg_length)
# print("press any key"); input()
print("outer radii", np.min(outer_r) , np.max(outer_r))
radii = r.rs.radii
print("inner radii", np.min(radii) , np.max(radii))

rho_ = np.divide(outer_r, np.array(radii))
print("rho", np.min(rho_) , np.max(rho_))

radius_min, radius_max = np.min(radii) , np.max(radii)
kr_max = 0.000181
kr_min = 0.0000173
print("inner_r*kr", kr_min * radius_min, kr_max * radius_max)

sra_table_lookup = open_sra_lookup("../table_jan2")

# quick check
# rsx2 = soil_root_inerface(np.array([-15000]), np.array([-700]predefined_growth), r, sp, outer_r)
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
        # print(c, "iterations", wall_interpolation / (wall_interpolation + wall_xylem), wall_xylem / (wall_interpolation + wall_xylem))

        wall_fixpoint = timeit.default_timer() - wall_fixpoint

        # rx = r.solve(rs_age + t, -trans * sinusoidal(t), 0., rsx, False, wilting_point, [])  # xylem_flux.py, cells = False
        fluxes = r.segFluxes(rs_age + t, rx, rsx, False)  # class XylemFlux is defined in MappedOrganism.h, approx = True
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
