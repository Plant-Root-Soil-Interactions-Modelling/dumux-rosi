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
from xylem_flux import *  # root system Python hybrid solver
from rhizo_models import *  # Helper class for cylindrical rhizosphere models

import vtk_plot as vp
import van_genuchten as vg
from root_conductivities import *
from sra_table_lookup import *

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def soil_root_interface_table2(rx, sx, inner_kr_, rho_, f):
    assert rx.shape == sx.shape
    rsx = f((rx, sx, inner_kr_ , rho_))
    return rsx

""" 
Parameters  
"""

""" soil """
name = "singleroot"  # name to export resutls
min_b = [-0.5, -0.5, -50.]
max_b = [0.5, 0.5, 0.]
cell_number = [1, 1, 50]  # # full is very slow
periodic = False  # check data first
domain_volume = np.prod(np.array(max_b) - np.array(min_b))
alpha = 0.018;  # (cm-1)
n = 1.8;
Ks = 28.46;  # (cm d-1)
loam = [0.08, 0.43, alpha, n, Ks]
p_top = -300  # -5000 (dry), -310 (wet)
p_bot = -200  #
soil_ = loam
soil = vg.Parameters(soil_)
sra_table_lookup = open_sra_lookup("../table_jan2")

""" root system """
collar = -8000  # dirichlet
radius = 0.05  # cm
wilting_point = -15000

""" rhizosphere models """
mode = "dumux"
NC = 10  # dof+1
logbase = 0.5  # according to Mai et al. (2019)
split_type = 0  # type 0 == volume, type 1 == surface, type 2 == length

""" simulation time """
sim_time = 1.  # 0.65  # 0.25  # [day]
dt = 60 / (24 * 3600)  # time step [day], 120 schwankt stark
NT = int(np.ceil(sim_time / dt))  # number of iterations
skip = 6 * 60  # for output and results, skip iteration

""" 
Initialize macroscopic soil model (Dumux binding)
"""
s = RichardsWrapper(RichardsSP())
s.initialize()
s.createGrid(min_b, max_b, cell_number, periodic)  # [cm]
# s.setHomogeneousIC(initial_sp)
s.setLinearIC(p_top, p_bot)  # cm pressure head, equilibrium
s.setTopBC("noFlux")
s.setBotBC("noFlux")
s.setVGParameters([soil_])
s.setParameter("Newton.EnableAbsoluteResidualCriterion", "True")
# s.setParameter("Soil.SourceSlope", "1000")  # turns regularisation of the source term on
s.initializeProblem()
s.setCriticalPressure(wilting_point)  # new source term regularisation
s.ddt = 1.e-5  # [day] initial Dumux time step
print()

""" 
Initialize xylem model 
"""
ns = 100  # 50 cm root, 100 segments, 0.5 cm each
radii = np.array([radius] * ns)
nodes = [pb.Vector3d(0, 0, 0)]
segs = []
for i in range(0, 100):
    nodes.append(pb.Vector3d(0, 0, -(i + 1) * 0.5))
    segs.append(pb.Vector2i(i, i + 1))

ms = pb.MappedSegments(nodes, segs, radii)
rs = RhizoMappedSegments(ms, wilting_point, NC, logbase, mode)
rs.setRectangularGrid(pb.Vector3d(min_b[0], min_b[1], min_b[2]), pb.Vector3d(max_b[0], max_b[1], max_b[2]),
                        pb.Vector3d(cell_number[0], cell_number[1], cell_number[2]), True)
r = XylemFluxPython(rs)  # wrap the xylem model around the MappedSegments
init_singleroot_contkrkx(r)
picker = lambda x, y, z: s.pick([x, y, z])  #  function that return the index of a given position in the soil grid (should work for any grid - needs testing)
rs.setSoilGrid(picker)  # maps segments, maps root segements and soil grid indices to each other in both directions
rs.set_xylem_flux(r)

mapping = np.array([r.rs.seg2cell[j] for j in range(0, ns)])
outer_r = r.rs.segOuterRadii()
inner_r = r.rs.radii
types = r.rs.subTypes
rho_ = np.divide(outer_r, np.array(inner_r))

""" sanity checks """
r.test()  # sanity checks

""" 
Initialize local soil models (around each root segment) 
"""
start_time = timeit.default_timer()
x = s.getSolutionHead()  # initial condition of soil [cm]
rs.initialize(soil_, x)
# print("press any key"); input()

""" Numerical solution (a) """
start_time = timeit.default_timer()
x_, y_, w_, cpx, cps, cf = [], [], [], [], [], []
sink1d = []
sink1d2 = []

# for post processing
out_times = []  # days
psi_x_ = []
psi_s_ = []
sink_ = []

sx = s.getSolutionHead()  # inital condition, solverbase.py
hsb = np.array([sx[mapping[j]][0] for j in range(0, ns)])  # soil bulk matric potential per segment
kr_ = np.zeros((ns,))
rsx = hsb.copy()  # initial values for fix point iteration

t = 0.
rs_age = 0.

for i in range(0, NT):

    wall_iteration = timeit.default_timer()
    wall_fixpoint = timeit.default_timer()

    if i == 0:  # only first time
        rx = r.solve_dirichlet(rs_age + t, [collar], 0., rsx, cells = False, soil_k = [])
        # rx = r.solve(rs_age + t, -trans * sinusoidal(t), 0., rsx, False, wilting_point, [])  # xylem_flux.py, cells = False
        rx_old = rx.copy()

    err = 1.e6
    c = 1
    while err > 1. and c < 100:

        """ interpolation """
        wall_interpolation = timeit.default_timer()

        for j in range(0, len(outer_r)):  # determine kr at this time step
            kr_[j] = r.kr_f(rs_age + t, types[j])

        inner_kr_ = np.multiply(inner_r, kr_)  # multiply for table look up
        rsx = soil_root_interface_table2(rx[1:], hsb, inner_kr_, rho_, sra_table_lookup)
        wall_interpolation = timeit.default_timer() - wall_interpolation

        """ xylem matric potential """
        wall_xylem = timeit.default_timer()
        rx = r.solve_dirichlet(rs_age + t, [collar], 0., rsx, cells = False, soil_k = [])
        # rx = r.solve(rs_age + t, -trans * sinusoidal(t), 0., rsx, False, wilting_point, [])  # xylem_flux.py, cells = False
        err = np.linalg.norm(rx - rx_old)
        wall_xylem = timeit.default_timer() - wall_xylem
        # print(err)
        rx_old = rx.copy()
        c += 1
#         print(c, ": ", np.sum(rx[1:]), np.sum(hsb), np.sum(inner_kr_), np.sum(rho_))
        print(c, "iterations", wall_interpolation / (wall_interpolation + wall_xylem), wall_xylem / (wall_interpolation + wall_xylem))

        wall_fixpoint = timeit.default_timer() - wall_fixpoint

    fluxes = r.segFluxes(rs_age + t, rx, rsx, False)
    min_rsx = np.min(rsx)  # for console output
    max_rsx = np.max(rsx)

    wall_soil = timeit.default_timer()

    soil_fluxes = r.sumSegFluxes(fluxes)
    s.setSource(soil_fluxes.copy())  # richards.py
    s.solve(dt)
    sx = s.getSolutionHead()  # richards.py
    hsb = np.array([sx[mapping[j]][0] for j in range(0, ns)])  # soil bulk matric potential per segment
    water = s.getWaterVolume()
    wall_soil = timeit.default_timer() - wall_soil

    wall_iteration = timeit.default_timer() - wall_iteration

    """ remember results ... """
    if i % skip == 0:
        print(i / skip)

        rx_ = 0.5 * (rx[0:-1] + rx[1:])  # psix is given per node, i convert to per segment
        psi_x_.append(rx_)

        psi_s_.append(rsx.copy())

        # sink_.append(np.array(rs.get_inner_fluxes()))
        sink_.append(fluxes)

""" xls file output """
if rank == 0:

    file1 = 'psix_singleroot_sra_constkrkx_wet.xls'
    df1 = pd.DataFrame(np.transpose(np.array(psi_x_)))
    df1.to_excel(file1, index = False)

    file2 = 'psiinterface_singleroot_sra_constkrkx_wet.xls'
    df2 = pd.DataFrame(np.transpose(np.array(psi_s_)))
    df2.to_excel(file2, index = False)

    file3 = 'sink_singleroot_sra_constkrkx_wet.xls'
    df3 = pd.DataFrame(-np.transpose(np.array(sink_)))
    df3.to_excel(file3, index = False)

