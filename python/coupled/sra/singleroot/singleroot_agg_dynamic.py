""" 
Single root scenario: Soil depletion due to sinusoidal transpiration

using an aggregated root system, and steady rate approach and fix point iteration
"""

import sys; sys.path.append("../../../modules/"); sys.path.append("../../../../../CPlantBox/");  sys.path.append("../../../../../CPlantBox/src/python_modules")
sys.path.append("../../../../build-cmake/cpp/python_binding/"); sys.path.append("../../../modules/fv/");
sys.path.append("../"); sys.path.append("../scenarios/");

import plantbox as pb  # CPlantBox
from rosi_richards import RichardsSP  # C++ part (Dumux binding), macroscopic soil model
from richards import RichardsWrapper  # Python part, macroscopic soil model
from xylem_flux import *  # root system Python hybrid solver
from rhizo_models import *  # Helper class for cylindrical rhizosphere models

import vtk_plot as vp
import van_genuchten as vg
from detach2 import *
from root_conductivities import *
from sra_table_lookup import *
import aggregated_rs as agg

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
p_top = -5000  # -5000 (_dry), -1000 (_wet)
p_bot = -200  #
sstr = "_dry"  # <---------------------------------------------------------- (dry or wet)

min_b = [-0.5, -0.5, -50.]  # domain
max_b = [0.5, 0.5, 0.]
cell_number = [1, 1, 100]
periodic = False
domain_volume = np.prod(np.array(max_b) - np.array(min_b))

alpha = 0.018  # (cm-1) soil
n = 1.8
Ks = 28.46  # (cm d-1)
loam = [0.08, 0.43, alpha, n, Ks]
soil_ = loam
soil = vg.Parameters(soil_)
vg.create_mfp_lookup(soil, -1.e5, 1000)  # creates the matrix flux potential look up table (in case for exact)
sra_table_lookup = open_sra_lookup("../table_jan2")  # opens the precomputed soil root interface potentials

""" root system """
trans = 0.5  # cm3/day
radius = 0.05  # cm
wilting_point = -10000

""" simulation time """
sim_time = 7.1  # 0.65  # 0.25  # [day]
dt = 60 / (24 * 3600)  # time step [day], 120 schwankt stark
NT = int(np.ceil(sim_time / dt))  # number of iterations
skip = 1  # for output and results, skip iteration

""" 
Initialize macroscopic soil model (Dumux binding)
"""
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
s.setCriticalPressure(wilting_point)  # new source term regularisation
s.ddt = 1.e-5  # [day] initial Dumux time step

""" 
Initialize xylem model 
"""
""" normal root system for krs, suf """

r = agg.create_singleroot()

agg.init_conductivities_const(r)

r = agg.create_aggregated_rs(r, 0., min_b, max_b, cell_number)
picker = lambda x, y, z: s.pick([0., 0., z])  #  function that return the index of a given position in the soil grid (should work for any grid - needs testing)
r.rs.setSoilGrid(picker)  # maps segments, maps root segements and soil grid indices to each other in both directions
# vp.plot_roots(pb.SegmentAnalyser(rs), "radius")
nodes = r.rs.nodes

""" CHECK """
# print(r.rs.nodes[0], r.rs.nodes[1], r.rs.nodes[2], r.rs.nodes[3], r.rs.nodes[-4], r.rs.nodes[-3], r.rs.nodes[-2], r.rs.nodes[-1])
# print(r.rs.segments[0], r.rs.segments[1], r.rs.segments[2], r.rs.segments[3], r.rs.segments[-4], r.rs.segments[-3], r.rs.segments[-2], r.rs.segments[-1])

""" sanity checks """
ns = len(r.rs.segments)
r.test()  # sanity checks
mapping = np.array([r.rs.seg2cell[j] for j in range(0, ns)])
outer_r = r.rs.segOuterRadii()
inner_r = r.rs.radii
types = r.rs.subTypes
rho_ = np.divide(outer_r, np.array(inner_r))

# print(outer_r)
# print(inner_r)
# print(types)
# print(rho_)

""" Numerical solution (a) """
start_time = timeit.default_timer()
# print(s.getCellCenters())

# for post processing
out_times = []  # days
psi_x_ = []
psi_s_ = []
psi_s2_ = []
sink_ = []
collar_vfr = []
sink_sum = []
x_, y_ = [], []

sx = s.getSolutionHead()  # inital condition, solverbase.py
cell_centers = s.getCellCenters()
hsb = np.array([sx[mapping[2 * j + 1]][0] for j in range(0, int(ns / 2))])  # soil bulk matric potential per segment
cell_centers_z = np.array([cell_centers[mapping[2 * j + 1]][2] for j in range(0, int(ns / 2))])
kr_ = np.zeros((ns,))
rsx = hsb.copy()  # initial values for fix point iteration
# print(s.getCellCenters())

t = 0.
rs_age = 0.

for i in range(0, NT):

    t = i * dt  # current simulation time

    wall_iteration = timeit.default_timer()
    wall_fixpoint = timeit.default_timer()

    if i == 0:  # only first time
        # rx = r.solve_dirichlet(rs_age + t, [collar], 0., rsx, cells = False, soil_k = [])
        rx = r.solve(rs_age + t, -trans * sinusoidal(t), 0., double_(rsx), False, wilting_point, soil_k=[])
        rx_old = rx.copy()

    for j in range(0, len(outer_r)):  # determine kr at this time step
        kr_[j] = r.kr_f(rs_age + t, types[j], 2, 2, j)
    inner_kr_ = np.multiply(inner_r, kr_)  # multiply for table look up

    err = 1.e6
    c = 1
    while err > 1 and c < 100:

        """ interpolation """
        # wall_interpolation = timeit.default_timer()
        rx_ = rx[1::2] - np.array([nodes[ii].z for ii in range(1, ns, 2)])  # from total matric potential to matric potential
        hsb_ = hsb - cell_centers_z  # from total matric potential to matric potential
        rsx = soil_root_interface_table2(rx_, hsb_, inner_kr_[1::2], rho_[1::2], sra_table_lookup)  # [1::2] every second entry, starting from 1
        rsx = rsx + np.array([nodes[ii].z for ii in range(1, ns, 2)])  # from matric potential to total matric potential
        # wall_interpolation = timeit.default_timer() - wall_interpolation

        """ xylem matric potential """
        # wall_xylem = timeit.default_timer()
        rx = r.solve(rs_age + t, -trans * sinusoidal(t), 0., double_(rsx), False, wilting_point, soil_k=[])  # xylem_flux.py, cells = False
        err = np.linalg.norm(rx - rx_old)
        # wall_xylem = timeit.default_timer() - wall_xylem
        rx_old = rx.copy()
        c += 1
#         print(c, ": ", np.sum(rx[1:]), np.sum(hsb), np.sum(inner_kr_), np.sum(rho_))
#        print(c, "iterations", wall_interpolation / (wall_interpolation + wall_xylem), wall_xylem / (wall_interpolation + wall_xylem))

#    wall_fixpoint = timeit.default_timer() - wall_fixpoint

    fluxes = r.segFluxes(rs_age + t, rx, double_(rsx), approx=False, cells=False)

#     min_rsx = np.min(rsx)  # for console output
#     max_rsx = np.max(rsx)
#     print("from", min_rsx, "to", max_rsx)

    wall_soil = timeit.default_timer()
    soil_fluxes = r.sumSegFluxes(fluxes)
    s.setSource(soil_fluxes.copy())  # richards.py
    s.solve(dt)
    sx = s.getSolutionHead()  # richards.py
    hsb = np.array([sx[mapping[2 * j + 1]][0] for j in range(0, int(ns / 2))])
    wall_soil = timeit.default_timer() - wall_soil

    wall_iteration = timeit.default_timer() - wall_iteration

    """ remember results ... """
    if i % skip == 0:
        print(i / skip)
        x_.append(t)
        sum_flux = 0.
        for f in soil_fluxes.values():
            sum_flux += f
        y_.append(sum_flux)  # cm3/day        
        rx_ = rx[1:]  # 0.5 * (rx[0:-1] + rx[1:])  # psix is given per node, converted to per segment
        psi_x_.append(rx_)
        psi_s_.append(rsx.copy())
        dd = np.array(sx)
        psi_s2_.append(dd[:, 0])
        sink_.append(fluxes[1::2])
        collar_vfr.append(r.collar_flux(0, rx.copy(), rsx.copy(), k_soil=[], cells=False))  # def collar_flux(self, sim_time, rx, sxx, k_soil=[], cells=True):
        sink_sum.append(np.sum(fluxes))

""" xls file output """
print("write xls")

file1 = 'results/psix_singleroot_agg_dynamic_constkrkx' + sstr + '.xls'
df1 = pd.DataFrame(np.array(psi_x_))
df1.to_excel(file1, index=False, header=False)

file2 = 'results/psiinterface_singleroot_agg_dynamic_constkrkx' + sstr + '.xls'
df2 = pd.DataFrame(np.array(psi_s_))
df2.to_excel(file2, index=False, header=False)

file3 = 'results/sink_singleroot_agg_dynamic_constkrkx' + sstr + '.xls'
df3 = pd.DataFrame(-np.array(sink_))
df3.to_excel(file3, index=False, header=False)

file4 = 'results/transpiration_singleroot_agg_dynamic_constkrkx' + sstr
np.savetxt(file4, np.vstack((x_, -np.array(y_))), delimiter=';')

file5 = 'results/soil_singleroot_agg_dynamic_constkrkx' + sstr + '.xls'
df5 = pd.DataFrame(np.array(psi_s2_))
df5.to_excel(file5, index=False, header=False)

print("fin")
# print(collar_vfr)
# print(sink_sum)

