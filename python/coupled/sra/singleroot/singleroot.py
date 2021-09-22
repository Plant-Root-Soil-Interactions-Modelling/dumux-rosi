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

""" root system """
rs_age = 0.  # in case of age dep conductivitities
# we are using -8000 dirichle
trans = 0.5  # * 15 * 75  # average per day [cm3 /day] (sinusoidal),
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

    t = i * dt  # current simulation time

    """ 1. xylem model """
    wall_root_model = timeit.default_timer()
    rsx = rs.get_inner_heads()  # matric potential at the root soil interface, i.e. inner values of the cylindric models (not extrapolation to the interface!) [cm]
    soil_k = np.divide(vg.hydraulic_conductivity(rsx, soil), rs.radii)  # only valid for homogenous soil
    # rx = r.solve(rs_age + t, -trans, 0., rsx, False, wilting_point, soil_k)  # [cm]   False means that rsx is given per root segment not per soil cell
    rx = r.solve_dirichlet(rs_age + t, [-8000.], 0., rsx, cells = False, soil_k = soil_k)
    proposed_inner_fluxes = r.segFluxes(rs_age + t, rx, rsx, approx = False, cells = False, soil_k = soil_k)  # [cm3/day]
    wall_root_model = timeit.default_timer() - wall_root_model
    # print(proposed_inner_fluxes)

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

        rsx_ = 0.5 * (rx[0:-1] + rx[1:])  # psix is given per node, i convert to per segment
        psi_x_.append(rsx_)

        psi_s = rs.get_inner_heads()
        psi_s_.append(psi_s)

        # sink_.append(np.array(rs.get_inner_fluxes()))
        sink_.append(proposed_inner_fluxes)

""" plots and output """
if rank == 0:
    print ("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s")

    # np.savetxt(name, np.vstack((out_times, -np.array(collar_flux), -np.array(water_uptake))), delimiter = ';')
    # sink1d = np.array(sink1d)
    # np.save(name + "_sink", sink1d)
    #
    # rsx = rs.get_inner_heads()  # matric potential at the root soil interface, i.e. inner values of the cylindric models [cm]
    # soil_k = np.divide(vg.hydraulic_conductivity(rsx, soil), rs.radii)  # only valid for homogenous soil
    # rx = r.solve(rs_age + t, -trans * sinusoidal(t), 0., rsx, False, wilting_point, soil_k)
    # fluxes = r.segFluxes(rs_age + t, rx, rsx, approx = False, cells = False, soil_k = soil_k)
    # # ana = pb.SegmentAnalyser(r.rs)
    # # ana.addData("rsx", rsx)
    # # ana.addData("rx", rx)
    # # ana.addData("fluxes", fluxes)
    # # vp.plot_roots(ana, "rsx")  # VTK vizualisation
    # # vp.plot_roots(ana, "fluxes")  # VTK vizualisation
    #
    # crit_min_i, crit_max_i, crit_min_o, crit_max_o = rs.plot_cylinders()
    # print(crit_min_i)
    # rs.plot_cylinder(crit_min_i)
    # print(rs.radii[crit_min_i])
    #
    # plot_transpiration(out_times, water_uptake, collar_flux, lambda t: trans * sinusoidal(t))  # in rhizo_models.py
    # plot_info(out_times, water_collar_cell, water_cyl, collar_sx, min_sx, min_rx, min_rsx, water_uptake, water_domain)  # in rhizo_models.py
    #
    # # vp.plot_roots_and_soil(r.rs, "pressure head", rsx, s, periodic, min_b, max_b, cell_number, name)  # VTK vizualisation

    # # convert your array into a dataframe
    file1 = 'psix_singleroot_constkrkx_wet.xls'
    df1 = pd.DataFrame(np.transpose(np.array(psi_x_)))
    df1.to_excel(file1, index = False)

    file2 = 'psiinterface_singleroot_constkrkx_wet.xls'
    df2 = pd.DataFrame(np.transpose(np.array(psi_s_)))
    df2.to_excel(file2, index = False)

    file3 = 'sink_singleroot_constkrkx_wet2.xls'
    df3 = pd.DataFrame(-np.transpose(np.array(sink_)))
    df3.to_excel(file3, index = False)

