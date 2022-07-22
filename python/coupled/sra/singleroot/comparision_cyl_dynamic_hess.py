""" 
Single root scenario: Soil depletion due to sinusoidal transpiration

coupled to cylindrical rhizosphere models using 1d axi-symetric richards equation (DUMUX solver)
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
import aggregated_rs as agg

import matplotlib.pyplot as plt
import numpy as np

""" 
Parameters  
"""

""" soil """
p_top = -330
p_bot = -180
sstr = "_hess"

min_b = [-1, -1, -150.]  # domain
max_b = [1, 1, 0.]
cell_number = [1, 1, 150]
periodic = False

alpha = 0.0383  # (cm-1) soil
n = 1.3774
Ks = 60.  # (cm d-1)
loam = [0.025, 0.403, alpha, n, Ks]
soil_ = loam
soil = vg.Parameters(soil_)
vg.create_mfp_lookup(soil, -1.e5, 1000)  # creates the matrix flux potential look up table (in case for exact)

""" root system """
trans = 0.6 * 4  # cm3/day
radius = 0.05  # cm
wilting_point = -15000

""" rhizosphere models """
mode = "dumux"
NC = 10  # dof+1, i.e. dof = 9
logbase = 0.5  # according to Mai et al. (2019)
split_type = 0  # type 0 == volume, type 1 == surface, type 2 == length

""" simulation time """
sim_time = 21  # 0.65  # 0.25  # [day]
dt = 20 / (24 * 3600)  # time step [day], 120 schwankt stark
NT = int(np.ceil(sim_time / dt))  # number of iterations
skip = 3  # for output and results, skip iteration

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
ns = 100
r = agg.create_singleroot(ns, 100, 0.05)
r.rs.setRectangularGrid(pb.Vector3d(min_b[0], min_b[1], min_b[2]), pb.Vector3d(max_b[0], max_b[1], max_b[2]),
                        pb.Vector3d(cell_number[0], cell_number[1], cell_number[2]), cut = False)
picker = lambda x, y, z: s.pick([x, y, z])  #  function that return the index of a given position in the soil grid (should work for any grid - needs testing)
r.rs.setSoilGrid(picker)  # maps segments, maps root segements and soil grid indices to each other in both directions
nodes = r.rs.nodes

rs = RhizoMappedSegments(r.rs, wilting_point, NC, logbase, mode)
rs.setRectangularGrid(pb.Vector3d(min_b[0], min_b[1], min_b[2]), pb.Vector3d(max_b[0], max_b[1], max_b[2]),
                        pb.Vector3d(cell_number[0], cell_number[1], cell_number[2]), True)
picker = lambda x, y, z: s.pick([x, y, z])  #  function that return the index of a given position in the soil grid (should work for any grid - needs testing)
rs.setSoilGrid(picker)  # maps segments, maps root segements and soil grid indices to each other in both directions
rs.set_xylem_flux(r)
agg.init_comp_conductivities_const(r)

r.linearSystem = r.linearSystem_doussan  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11

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
psi_s_ = []  # root soil interface
psi_s2_ = []  # soil
sink_ = []
collar_vfr = []
sink_sum = []
x_, y_ = [], []

water_uptake, water_collar_cell, water_cyl, water_domain = [], [], [], []  # cm3

cci = picker(rs.nodes[0].x, rs.nodes[0].y, rs.nodes[0].z)  # collar cell index
cell_volumes = s.getCellVolumes()  # cm3
cell_volumes = comm.bcast(cell_volumes, root = 0)
net_flux = np.zeros(cell_volumes.shape)

t = 0.
rs_age = 0.

for i in range(0, NT + 1):

    t = i * dt  # current simulation time

    """ 1. xylem model """
    rsx = rs.get_inner_heads()  # matric potential at the root soil interface, i.e. inner values of the cylindric models (not extrapolation to the interface!) [cm]
    if i == 0:
        # rx = r.solve_dirichlet(0., [collar], 0., rsx.copy(), cells = False, soil_k = [])
        rx = r.solve(rs_age + t, -trans * sinusoidal2(t, dt), 0., rsx, cells = False, wilting_point = wilting_point, soil_k = [])

    # exact
    soil_k0 = np.zeros(rsx.shape)
    for j in range(0, rsx.shape[0]):
        hsoil = rsx[j]
        hint = rx[j + 1]
        soil_k0[j] = (vg.fast_mfp[soil](hsoil) - vg.fast_mfp[soil](hint)) / (hsoil - hint)
    soil_k00 = np.divide(soil_k0, rs.get_dx2())  # only valid for homogenous soil
    soil_k000 = rs.get_soil_k(rx)

    # approx
    soil_k1 = np.divide(vg.hydraulic_conductivity(rsx, soil), rs.get_dx2())  # only valid for homogenous soil

    # approx2 rs.get_dx2()
    soil_k2 = np.divide(vg.hydraulic_conductivity(rx[1:], soil), rs.get_dx2())  # only valid for homogenous soil

    # old
    soil_k = np.divide(vg.hydraulic_conductivity(rsx, soil), rs.radii)  # only valid for homogenous soil

    # if i % skip == 0:
    #     print("exact  ", np.min(soil_k000), np.max(soil_k000))
    #     print("exact  ", np.min(soil_k00), np.max(soil_k00))
    #     print("approx1", np.min(soil_k1), np.max(soil_k1))
    #     print("approx2", np.min(soil_k2), np.max(soil_k2))
    #     print("old    ", np.min(soil_k), np.max(soil_k), "\n")

    soil_k = soil_k00

    # rx = r.solve_dirichlet(0., [collar], 0., rsx.copy(), cells = False, soil_k = soil_k)
    rx = r.solve(rs_age + t, -trans * sinusoidal2(t, dt), 0., rsx, cells = False, wilting_point = wilting_point, soil_k = soil_k)

    proposed_inner_fluxes = r.segFluxes(0., rx.copy(), rsx.copy(), approx = True, cells = False, soil_k = soil_k.copy())  # [cm3/day]
    # collar_flux = r.collar_flux(0., rx.copy(), rsx.copy(), k_soil = soil_k.copy(), cells = False)
    #
    # err = np.linalg.norm(np.sum(proposed_inner_fluxes) - collar_flux)
    #
    # if err > 1.e-10:
    #     print("error" , err)
    #     print(rx)
    #     print(np.min(rx), np.max(rx))
    #     print(np.min(rsx), np.max(rsx))
    #     print(np.min(soil_k), np.max(soil_k))
    #     print(np.min(proposed_inner_fluxes), np.max(proposed_inner_fluxes))
    #     print(collar_flux)

    """ 2. local soil models """
    proposed_outer_fluxes = r.splitSoilFluxes(net_flux / dt, split_type)
    # print(np.min(proposed_outer_fluxes), np.max(proposed_outer_fluxes), np.sum(proposed_outer_fluxes))
    rs.solve(dt, proposed_inner_fluxes, proposed_outer_fluxes)  # left and right neumann fluxes
    realized_inner_fluxes = rs.get_inner_fluxes()  # identical for mode = "dumux"

    err = np.linalg.norm(np.array(proposed_inner_fluxes) - np.array(realized_inner_fluxes))
    if err > 1.e-15:
        print("relative error" , err)

    """ 3a. macroscopic soil model """
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
        # net_flux[k] = 0  # same as SRA approach (currently no-flux at boundary)
    soil_water = new_soil_water

    """ remember results ... """
    if i % skip == 0:
        x_.append(t)
        sum_flux = 0.
        for f in soil_fluxes.values():
            sum_flux += f
        y_.append(sum_flux)  # cm3/day
        fluxes = np.array(proposed_inner_fluxes)
        collar_flux = r.collar_flux(0., rx, rsx, k_soil = soil_k, cells = False)
        print(i / skip, collar_flux, np.sum(fluxes), np.min(fluxes), np.max(fluxes))
        rx_ = rx[1:]
        psi_x_.append(rx_)
        psi_s_.append(np.array(rsx.copy()))
        dd = np.array(s.getSolutionHead())
        psi_s2_.append(dd[:, 0])
        sink_.append(fluxes.copy())
        collar_vfr.append(collar_flux)  # def collar_flux(self, sim_time, rx, sxx, k_soil=[], cells=True):
        sink_sum.append(np.sum(fluxes))

print ("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s")

""" xls file output """
file1 = 'results/psix_singleroot_cyl_dynamic_constkrkx' + sstr
np.save(file1, np.array(psi_x_))

file2 = 'results/psiinterface_singleroot_cyl_dynamic_constkrkx' + sstr
np.save(file2, np.array(psi_s_))

file3 = 'results/sink_singleroot_cyl_dynamic_constkrkx' + sstr
np.save(file3, -np.array(sink_))

file4 = 'results/transpiration_singleroot_cyl_dynamic_constkrkx' + sstr
np.save(file4, np.vstack((x_, -np.array(y_))), delimiter = ';')

file5 = 'results/soil_singleroot_cyl_dynamic_constkrkx' + sstr
np.save(file5, np.array(psi_s2_))

print("fin")
