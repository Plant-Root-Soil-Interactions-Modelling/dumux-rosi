""" 
Jan's new scenario

Coupled to cylindrical rhizosphere models using 1d richards equation (DUMUX solver)

complicated MPI support (a non-mpi version of richards_cyl is needed, see script dumux3_nompi.sh)
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

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

""" 
Parameters  
"""

""" soil """
name = "singleroot"  # name to export resutls
min_b = [-0.5, -0.5, -50.]
max_b = [0.5, 0.5, 0.]
cell_number = [1, 1, 100]  # # full is very slow
periodic = False  # check data first
domain_volume = np.prod(np.array(max_b) - np.array(min_b))
alpha = 0.018;  # (cm-1)
n = 1.8;
Ks = 28.46;  # (cm d-1)
loam = [0.08, 0.43, alpha, n, Ks]
p_top = -5000  # -5000 (dry), -310 (wet)
p_bot = -200  #
sstr = "_dry"  # <---------------------------------------------------------- (dry or wet)
soil_ = loam
soil = vg.Parameters(soil_)

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
sim_time = 0.51  # 0.65  # 0.25  # [day]
dt = 2 / (24 * 3600)  # time step [day], 120 schwankt stark
NT = int(np.ceil(sim_time / dt))  # number of iterations
skip = 30 * 1 * 60  # for output and results, skip iteration

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
n = 100  # 50 cm root, 100 segments, 0.5 cm each
radii = np.array([radius] * n)
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

""" sanity checks """
r.test()  # sanity checks

""" 
Initialize local soil models (around each root segment) 
"""
start_time = timeit.default_timer()
x = s.getSolutionHead()  # initial condition of soil [cm]
rs.initialize(soil_, x)
# print("press any key"); input()

""" 
Simulation 

loop
1. xylem model
2. local soil models
3. macroscopic soil model 
"""
print("Starting simulation")
start_time = timeit.default_timer()

# for post processing
out_times = []  # days
psi_x_ = []
psi_s_ = []
sink_ = []

water_uptake, water_collar_cell, water_cyl, water_domain = [], [], [], []  # cm3

cci = picker(rs.nodes[0].x, rs.nodes[0].y, rs.nodes[0].z)  # collar cell index
cell_volumes = s.getCellVolumes()  # cm3
cell_volumes = comm.bcast(cell_volumes, root = 0)
net_flux = np.zeros(cell_volumes.shape)

for i in range(0, NT + 1):

    wall_iteration = timeit.default_timer()

    """ 1. xylem model """
    wall_root_model = timeit.default_timer()
    rsx = rs.get_inner_heads()  # matric potential at the root soil interface, i.e. inner values of the cylindric models (not extrapolation to the interface!) [cm]
    rsx_old = rsx.copy()
    soil_k = np.divide(vg.hydraulic_conductivity(rsx, soil), rs.radii)  # only valid for homogenous soil
    rx = r.solve_dirichlet(0., [collar], 0., rsx, cells = False, soil_k = soil_k)
    proposed_inner_fluxes = r.segFluxes(0., rx, rsx, approx = False, cells = False, soil_k = soil_k)  # [cm3/day]
    wall_root_model = timeit.default_timer() - wall_root_model

    """ 2. local soil models """
    wall_rhizo_models = timeit.default_timer()

    proposed_outer_fluxes = r.splitSoilFluxes(net_flux / dt, split_type)
    rs.solve(dt, proposed_inner_fluxes, proposed_outer_fluxes)  # left and right neumann fluxes
    realized_inner_fluxes = rs.get_inner_fluxes()  # identical for mode = "dumux"

    wall_rhizo_models = timeit.default_timer() - wall_rhizo_models

    """ 3a. macroscopic soil model """
    wall_soil_model = timeit.default_timer()

    water_content = np.array(s.getWaterContent())  # theta per cell [1]
    soil_water = np.multiply(water_content, cell_volumes)  # water per cell [cm3]
    soil_fluxes = r.sumSegFluxes(realized_inner_fluxes)  # [cm3/day]  per soil cell
    s.setSource(soil_fluxes.copy())  # [cm3/day], in moduels/richards.py
    s.solve(dt)  # in modules/solverbase.py

    """ 3b. calculate net fluxes """
    water_content = np.array(s.getWaterContent())
    new_soil_water = np.multiply(water_content, cell_volumes)  # calculate net flux
    net_flux = new_soil_water - soil_water  # change in water per cell [cm3]
    for k, root_flux in soil_fluxes.items():
        net_flux[k] -= root_flux * dt
    soil_water = new_soil_water

    wall_soil_model = timeit.default_timer() - wall_soil_model
    wall_iteration = timeit.default_timer() - wall_iteration

    """ remember results ... """
    if i % skip == 0:
        print(i / skip)
        rx_ = rx[1:]
        psi_x_.append(rx_)
        psi_s_.append(np.array(rsx_old))
        sink_.append(np.array(rs.get_inner_fluxes()))

print ("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s")

""" xls file output """

file1 = 'results/psix_singleroot_cyl_constkrkx' + sstr + '.xls'
df1 = pd.DataFrame(np.transpose(np.array(psi_x_)))
df1.to_excel(file1, index = False, header = False)

file2 = 'results/psiinterface_singleroot_cyl_constkrkx' + sstr + '.xls'
df2 = pd.DataFrame(np.transpose(np.array(psi_s_)))
df2.to_excel(file2, index = False, header = False)

file3 = 'results/sink_singleroot_cyl_constkrkx' + sstr + '.xls'
df3 = pd.DataFrame(-np.transpose(np.array(sink_)))
df3.to_excel(file3, index = False, header = False)

