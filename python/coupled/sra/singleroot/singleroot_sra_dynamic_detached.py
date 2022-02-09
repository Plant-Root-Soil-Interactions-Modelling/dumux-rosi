""" 
Jan's single root scenario

steady rate approach (to check if numbers agree to Matlab) 
"""

import sys; sys.path.append("../../../modules/"); sys.path.append("../../../../../CPlantBox/");  sys.path.append("../../../../../CPlantBox/src/python_modules")
sys.path.append("../../../../build-cmake/cpp/python_binding/"); sys.path.append("../../../modules/fv/");
sys.path.append("../");

import plantbox as pb  # CPlantBox
from rosi_richards import RichardsSP  # C++ part (Dumux binding), macroscopic soil model
from richards import RichardsWrapper  # Python part, macroscopic soil model
from xylem_flux_detached import *  # root system Python hybrid solver
from rhizo_models import *  # Helper class for cylindrical rhizosphere models

import vtk_plot as vp
import van_genuchten as vg
from root_conductivities import *
from detach import *  # detached root conductivities
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

""" 
Parameters  
"""

""" soil """
name = "singleroot"  # name to export resutls
ns = 100  # 50 cm root, 100 segments, 0.5 cm each
L = 50.  # cm
min_b = [-0.5, -0.5, -L]
max_b = [0.5, 0.5, 0.]
cell_number = [1, 1, ns]  # # full is very slow
periodic = False  # check data first
domain_volume = np.prod(np.array(max_b) - np.array(min_b))
alpha = 0.018;  # (cm-1)
n = 1.8;
Ks = 28.46;  # (cm d-1)
loam = [0.08, 0.43, alpha, n, Ks]
p_top = -1000  # -5000 (_dry), -1000 (_wet)
p_bot = -200  #
sstr = "_wet"  # <---------------------------------------------------------- (dry or wet)
soil_ = loam
soil = vg.Parameters(soil_)
vg.create_mfp_lookup(soil, -1.e5, 1000)  # creates the matrix flux potential look up table (in case for exact)
sra_table_lookup = open_sra_lookup("../table_jan2")

""" root system """
# collar = -8000  # dirichlet
trans = 0.5  # 0.5
radius = 0.05  # cm
wilting_point = -10000

""" simulation time """
sim_time = 7.1  # 0.65  # 0.25  # [day]
dt = 60 / (24 * 3600)  # time step [day], 120 schwankt stark
NT = int(np.ceil(sim_time / dt))  # number of iterations
skip = 1 * 60 * 6  # for output and results, skip iteration

""" 
Initialize macroscopic soil model (Dumux binding)
"""
s = RichardsWrapper(RichardsSP())
s.initialize()
s.createGrid(min_b, max_b, cell_number, periodic)  # [cm]
s.setLinearIC(p_top, p_bot)  # cm pressure head, equilibrium
s.setTopBC("noFlux")
s.setBotBC("noFlux")
s.setVGParameters([soil_])
s.setParameter("Newton.EnableAbsoluteResidualCriterion", "True")
# s.setParameter("Soil.SourceSlope", "1000")  # turns regularisation of the source term on
s.initializeProblem()
s.setCriticalPressure(wilting_point)  # new source term regularisation
s.ddt = 1.e-5  # [day] initial Dumux time step

""" 
Initialize xylem model 
"""
radii = np.array([radius] * ns)
nodes = [pb.Vector3d(0, 0, 0)]
segs = []
for i in range(0, ns): 
    # print(((i + 1) / ns) * L)
    nodes.append(pb.Vector3d(0, 0, -((i + 1) / ns) * L))  # node i+1
    segs.append(pb.Vector2i(i, i + 1))

rs = pb.MappedSegments(nodes, segs, radii)
rs.setRectangularGrid(pb.Vector3d(min_b[0], min_b[1], min_b[2]), pb.Vector3d(max_b[0], max_b[1], max_b[2]),
                        pb.Vector3d(cell_number[0], cell_number[1], cell_number[2]), cut=False)
# connected root system
r = XylemFluxPython(rs)
init_singleroot_contkrkx(r)
suf = r.get_suf(sim_time=0.)  
krs, _ = r.get_krs(sim_time=0.)  
# detached root system
r = XylemFluxDetached(rs)  # wrap the xylem model around the MappedSegments
init_singleroot_contkrkx(r)
detached_conductivities(r, suf, krs)
picker = lambda x, y, z: s.pick([x, y, z])  #  function that return the index of a given position in the soil grid (should work for any grid - needs testing)
rs.setSoilGrid(picker)  # maps segments, maps root segements and soil grid indices to each other in both directions
# print("index collar", rs.segments[0].x)

""" sanity checks """
r.test()  # sanity checks
mapping = np.array([r.rs.seg2cell[j] for j in range(0, ns)])
outer_r = r.rs.segOuterRadii()
inner_r = r.rs.radii
types = r.rs.subTypes
rho_ = np.divide(outer_r, np.array(inner_r))

""" Numerical solution (a) """
start_time = timeit.default_timer()

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
hsb = np.array([sx[mapping[j]][0] for j in range(0, ns)])  # soil bulk matric potential per segment
kr_ = np.zeros((ns,))
rsx = hsb.copy()  # initial values for fix point iteration

t = 0.
rs_age = 0.

for i in range(0, NT):

    t = i * dt  # current simulation time

    wall_iteration = timeit.default_timer()
    wall_fixpoint = timeit.default_timer()

    if i == 0:  # only first time
        # rx = r.solve_dirichlet(rs_age + t, [collar], 0., rsx, cells = False, soil_k = [])
        rx = r.solve(rs_age + t, -trans * sinusoidal(t), 0., rsx, False, wilting_point, soil_k=[])
        rx_old = rx.copy()

    err = 1.e6
    c = 1

    for j in range(0, len(outer_r)):  # determine kr at this time step
        kr_[j] = r.kr_f(rs_age + t, types[j], 2, 2, j)
    inner_kr_ = np.multiply(inner_r, kr_)  # multiply for table look up

    while err > 1 and c < 100:

        """ interpolation """
        wall_interpolation = timeit.default_timer()

        mean_p2 = np.array([r.mean_xylem_pressure(ii, 0., rx, hsb, cells=False) for ii in range(0, len(segs))])
        # mean_p = 0.5 * (rx[1:] + rx[0] * np.ones(rx[1:].shape))
        # print(mean_p[0], mean_p2[0], rx[0], rx[1])
        rsx = soil_root_interface_table2(mean_p2, hsb, inner_kr_, rho_, sra_table_lookup)
        # rsx = soil_root_interface(rx[1:] , hsb, inner_kr_, rho_, soil)
        wall_interpolation = timeit.default_timer() - wall_interpolation

        """ xylem matric potential """
        wall_xylem = timeit.default_timer()
        # rx = r.solve_dirichlet(rs_age + t, [collar], 0., rsx, cells = False, soil_k = [])
        rx = r.solve(rs_age + t, -trans * sinusoidal(t), 0., rsx, False, wilting_point, soil_k=[])  # xylem_flux.py, cells = False
        err = np.linalg.norm(rx - rx_old)
        wall_xylem = timeit.default_timer() - wall_xylem
        # print(err)
        rx_old = rx.copy()
        c += 1
    # print(c, ": ", rx[0], np.min(rsx))  # np.sum(rx[1:]), np.sum(hsb), np.sum(inner_kr_), np.sum(rho_))
    # print(c, "iterations", rx[0])  # wall_interpolation / (wall_interpolation + wall_xylem), wall_xylem / (wall_interpolation + wall_xylem)
    wall_fixpoint = timeit.default_timer() - wall_fixpoint

    fluxes = r.segFluxes_detached(rs_age + t, rx, rsx, approx=False, cells=False)

    min_rsx = np.min(rsx)  # for console output
    max_rsx = np.max(rsx)
    # print("from", min_rsx, "to", max_rsx)

    wall_soil = timeit.default_timer()

    soil_fluxes = r.sumSegFluxes(fluxes)
    s.setSource(soil_fluxes.copy())  # richards.py
    s.solve(dt)
    sx = s.getSolutionHead()  # richards.py
    hsb = np.array([sx[mapping[j]][0] - nodes[segs[j].y].z for j in range(0, ns)])  # soil bulk matric potential per segment #
    water = s.getWaterVolume()
    wall_soil = timeit.default_timer() - wall_soil

    wall_iteration = timeit.default_timer() - wall_iteration

    """ remember results ... """
    if i % (skip / 60) == 0:
        x_.append(t)
        sum_flux = 0.
        min_flux = 1.e6
        max_flux = -1.e6
        for f in soil_fluxes.values():
            sum_flux += f  # TODO MINMAX
            min_flux = min(f, min_flux)  # for console output
            max_flux = max(f, max_flux)
        y_.append(sum_flux)  # cm3/day
        print("target", -trans * sinusoidal(t), c, "real sink", y_[-1], r.last, rx[0], rx[1], rx[-1])  #  "real collar", r.collar_flux(0, rx.copy(), rsx.copy(), k_soil=[], cells=False),
        # print("min flux", min_flux, "max flux", max_flux)
    if i % skip == 0:
        print(i / skip)
        rx_ = rx[1:]  # 0.5 * (rx[0:-1] + rx[1:])  # psix is given per node, converted to per segment
        psi_x_.append(0.5 * (rx[1:] + rx[0] * np.ones(rx[1:].shape)))  #  XXXXXXXXXXXXXXXXXXXXXXXx3
        psi_s_.append(rsx.copy())
        dd = np.array(s.getWaterContent())
#         print(s.getDofCoordinates())
        psi_s2_.append(dd[:, 0])
        sink_.append(fluxes.copy())
        # collar_vfr.append(r.collar_flux(0, rx.copy(), rsx.copy(), k_soil=[], cells=False))  # def collar_flux(self, sim_time, rx, sxx, k_soil=[], cells=True):
        # TODO collar flux is currently not working for detached
        sink_sum.append(np.sum(fluxes))

""" xls file output """

file1 = 'results/psix_singleroot_sra_dynamicA_constkrkx' + sstr + '.xls'
df1 = pd.DataFrame(np.transpose(np.array(psi_x_)))
df1.to_excel(file1, index=False, header=False)

file2 = 'results/psiinterface_singleroot_sra_dynamicA_constkrkx' + sstr + '.xls'
df2 = pd.DataFrame(np.transpose(np.array(psi_s_)))
df2.to_excel(file2, index=False, header=False)

file3 = 'results/sink_singleroot_sra_dynamicA_constkrkx' + sstr + '.xls'
df3 = pd.DataFrame(-np.transpose(np.array(sink_)))
df3.to_excel(file3, index=False, header=False)

file4 = 'results/transpiration_singleroot_sra_dynamicA_constkrkx' + sstr
np.savetxt(file4, np.vstack((x_, -np.array(y_))), delimiter=';')

file5 = 'results/soil_singleroot_sra_dynamicA_constkrkx' + sstr + '.xls'
df5 = pd.DataFrame(np.transpose(np.array(psi_s2_)))
df5.to_excel(file5, index=False, header=False)

print(sink_sum)

