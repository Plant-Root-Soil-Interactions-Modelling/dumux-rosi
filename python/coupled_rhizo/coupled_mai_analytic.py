import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src")

import os
import plantbox as pb  # CPlantBox
from rosi_richardsnc import RichardsNCSPanalytic as RichardsNCSP
from richards import RichardsWrapper  # Python part, macroscopic soil model
from functional.xylem_flux import *  # root system Python hybrid solver
from rhizo_models import *  # Helper class for cylindrical rhizosphere models

import visualisation.vtk_plot as vp
import functional.van_genuchten as vg
from functional.root_conductivities import *

import numpy as np
import timeit
import matplotlib.pyplot as plt
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank(); max_rank = comm.Get_size()


def write_file_array(name, data, space =","):
    name2 = './results/coupled_mai_analytic/'+ name+ '.txt'
    with open(name2, 'a') as log:
        log.write(space.join([num for num in map(str, data)])  +'\n')
results_dir="./results/coupled_mai_analytic/"
if not os.path.exists(results_dir):
    os.makedirs(results_dir)
else:
    test = os.listdir(results_dir)
    for item in test:
        try:
            os.remove(results_dir+item)
        except:
            pass
""" 
Parameters  
"""
name = "bauw_coupled_rhizo"

""" soil """
min_b = [-4., -4., -15.]  # cm
max_b = [4., 4., 0.]  # cm
domain_volume = np.prod(np.array(max_b) - np.array(min_b))
cell_number = [2, 2, 2]  
#[7, 7, 15]  # [8, 8, 15]  # [16, 16, 30]  # [32, 32, 60]  # [8, 8, 15] # [1]
periodic = False
loam = [0.08, 0.43, 0.04, 1.6, 50]
# loam = [0.03, 0.345, 0.01, 2.5, 28.6]
soil_ = loam
soil = vg.Parameters(soil_)
initial = -659.8 + (max_b[2] - min_b[2]) / 2  # -659.8 + 7.5 because -659.8 is the value at the top, but we need the average value in the domain

""" root system """
trans = 6.4  # average per day [cm3 /day] (sinusoidal)
wilting_point = -15000  # [cm]
age_dependent = False  # conductivities
predefined_growth = False  # root growth by setting radial conductivities
rs_age = 8 * (not predefined_growth) + 1 * predefined_growth  # rs_age = 0 in case of growth, else 8 days

""" rhizosphere models """
mode = "dumux_nc"  # or "dumux"
NC = 10  # dof+1
logbase = 1.5  # according to Mai et al. (2019)
split_type = 0  # type 0 == volume, type 1 == surface, type 2 == length

""" simulation time """
sim_time = 1  # 0.65  # 0.25  # [day]
dt = 30 / (24 * 3600)  # time step [day]
NT = int(np.ceil(sim_time / dt))  # number of iterations
skip = 1  # for output and results, skip iteration

""" 
Initialize macroscopic soil model (Dumux binding)
"""
s = RichardsWrapper(RichardsNCSP())
s.initialize()

s.setParameter("Component.MolarMass", "3.1e-2")  # TODO no idea, where this is neeeded, i don't want to use moles ever
s.setParameter("Component.LiquidDiffusionCoefficient", "6.e-10")  # m^2 s-1
s.setParameter("Component.freundlichN_", "124.8.")
s.setParameter("Component.freundlichK_", ".4")
# s.setParameter("Component.BufferPower", "140")  # buffer power = \rho * Kd [1]
s.setParameter("Soil.IC.C", "0.01")  # (mol)g / cm3  # TODO specialised setter?
s.setParameter("Soil.BC.Top.SType", "2")  # michaelisMenten=8 (SType = Solute Type)
s.setParameter("Soil.BC.Top.CValue", "0.")  # michaelisMenten=8 (SType = Solute Type)
# s.setParameter("Soil.BC.Top.SType", "1")  # michaelisMenten=8 (SType = Solute Type)
# s.setParameter("Soil.BC.Top.CValue", "0.007")  # michaelisMenten=8 (SType = Solute Type)
s.setParameter("Soil.BC.Bot.SType", "2")  # michaelisMenten=8 (SType = Solute Type)
s.setParameter("Soil.BC.Bot.CValue", "0.")
# s.setParameter("Soil.BC.Bot.SType", "8")  # michaelisMenten (SType = Solute Type)
# s.setParameter("RootSystem.Uptake.Vmax", s.dumux_str(3.26e-6 * 24 * 3600))  # (mol)g /cm^2 / s - > (mol)g /cm^2 / day
# s.setParameter("RootSystem.Uptake.Km", s.dumux_str(5.8e-3))  # (mol)g / cm3

s.createGrid(min_b, max_b, cell_number, periodic)  # [cm]
s.setHomogeneousIC(initial, True)  # cm pressure head, equilibrium
s.setTopBC("noFlux")
s.setBotBC("noFlux")
s.setVGParameters([soil_])
s.setParameter("Newton.EnableAbsoluteResidualCriterion", "True")
s.setParameter("Soil.SourceSlope", "1000")  # turns regularisation of the source term on
s.initializeProblem()
s.setCriticalPressure(wilting_point)  # new source term regularisation
s.ddt = 1.e-5  # [day] initial Dumux time step
print()

""" 
Initialize xylem model 
"""
rs = RhizoMappedSegments("../../grids/RootSystem8.rsml", wilting_point, NC, logbase, mode)
rs.setRectangularGrid(pb.Vector3d(min_b[0], min_b[1], min_b[2]), pb.Vector3d(max_b[0], max_b[1], max_b[2]),
                        pb.Vector3d(cell_number[0], cell_number[1], cell_number[2]), True)
# True: root segments are cut  to the soil grid so that each segment is completely within one soil control element, this works only for rectangular grids so far
picker = lambda x, y, z: s.pick([x, y, z])  #  function that return the index of a given position in the soil grid (should work for any grid - needs testing)
rs.setSoilGrid(picker)  # maps segments, maps root segements and soil grid indices to each other in both directions
r = XylemFluxPython(rs)  # wrap the xylem model around the MappedSegments
init_conductivities(r, age_dependent)  # age_dependent is a boolean, root conductivies are given in the file /root_conductivities.py
rs.set_xylem_flux(r)
print()

# # For debugging
# r.plot_conductivities()
# r.test()  # sanity checks (todo need improvements...)

""" 
Initialize local soil models (around each root segment) 
"""
start_time = timeit.default_timer()
x = s.getSolutionHead()  # initial condition of soil [cm]
x = comm.bcast(x, root = 0)  # Soil part runs parallel
ns = len(rs.segments)
dcyl = int(np.floor(ns / max_rank))
if rank + 1 == max_rank:
    rs.initialize(soil_, x, np.array(range(rank * dcyl, ns)))
    print ("Initialized rank {:g}/{:g} [{:g}-{:g}] in {:g} s".format(rank + 1, max_rank, rank * dcyl, ns, timeit.default_timer() - start_time))
else:
    rs.initialize(soil_, x, np.array(range(rank * dcyl, (rank + 1) * dcyl)))
    print ("Initialized rank {:g}/{:g} [{:g}-{:g}] in {:g} s".format(rank + 1, max_rank, rank * dcyl, (rank + 1) * dcyl, timeit.default_timer() - start_time))
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
min_sx, min_rx, min_rsx, collar_sx, collar_flux = [], [], [], [], []  # cm
water_uptake, water_collar_cell, water_cyl, water_domain, solute_uptake = [], [], [], [], []  # cm3
out_times = []  # days
cci = picker(rs.nodes[0].x, rs.nodes[0].y, rs.nodes[0].z)  # collar cell index
cell_volumes = s.getCellVolumes()  # cm3
cell_volumes = comm.bcast(cell_volumes, root = 0)
net_flux = np.zeros(cell_volumes.shape)

for i in range(0, NT):

    wall_iteration = timeit.default_timer()

    t = i * dt  # current simulation time

    """ 1. xylem model """
    wall_root_model = timeit.default_timer()
    rsx = rs.get_inner_heads()  # matric potential at the root soil interface, i.e. inner values of the cylindric models [cm]
    soil_k = np.divide(vg.hydraulic_conductivity(rsx, soil), rs.radii)  # only valid for homogenous soil
    if rank == 0:
        rx = r.solve(rs_age + t, -trans * sinusoidal(t), 0., rsx, False, wilting_point, soil_k)  # [cm]   False means that rsx is given per root segment not per soil cell
        proposed_inner_fluxes = r.segFluxes(rs_age + t, rx, rsx, approx = False, cells = False, soil_k = soil_k)  # [cm3/day]
    else:
        proposed_inner_fluxes = None
        rx = None
    wall_root_model = timeit.default_timer() - wall_root_model

    """ 2. local soil models """
    wall_rhizo_models = timeit.default_timer()
    proposed_outer_fluxes = r.splitSoilFluxes(net_flux / dt, split_type)
    #
    # TODO: outer solute fluxes
    #
    if rs.mode == "dumux_nc":
        proposed_inner_fluxes = comm.bcast(proposed_inner_fluxes, root = 0)
        rs.solve(dt, proposed_inner_fluxes, proposed_outer_fluxes)  # mass fluxes?
    else:
        print(rs.mode)
        raise Exception("this script is for dumux_nc only")
    realized_inner_fluxes = rs.get_inner_fluxes()
    realized_inner_fluxes = comm.bcast(realized_inner_fluxes, root = 0)
    realized_mass_fluxes = rs.get_inner_mass_fluxes()
    realized_mass_fluxes = comm.bcast(realized_mass_fluxes, root = 0)
    wall_rhizo_models = timeit.default_timer() - wall_rhizo_models

    """ 3a. macroscopic soil model """
    wall_soil_model = timeit.default_timer()
    water_content = np.array(s.getWaterContent())
    water_content = comm.bcast(water_content, root = 0)
    soil_water = np.multiply(water_content, cell_volumes)
    solute_conc = np.array(s.getSolution(1))  # solute concentration
    solutes = np.multiply(solute_conc, cell_volumes)  # water per cell [cm3]
    soil_fluxes = r.sumSegFluxes(realized_inner_fluxes)  # [cm3/day]  per soil cell
    s.setSource(soil_fluxes.copy())  # [cm3/day], in richards.py
    soil_mass_fluxes = r.sumSegFluxes(realized_mass_fluxes)  # [cm3/day]  per soil cell
    s.setSource(soil_mass_fluxes.copy(), 1)  # TODO UNITS ???? [cm3/day], in richards.py
    s.solve(dt)  # in solverbase.py

    """ 3b. calculate net fluxes """
    water_content = np.array(s.getWaterContent())
    water_content = comm.bcast(water_content, root = 0)
    new_soil_water = np.multiply(water_content, cell_volumes)  # calculate net flux
    net_flux = new_soil_water - soil_water  # change in water per cell [cm3]
    for k, root_flux in soil_fluxes.items():
        net_flux[k] -= root_flux * dt

    solute_conc = np.array(s.getSolution(1))  # solute concentration
    solute_conc = comm.bcast(solute_conc, root = 0)
    new_solutes = np.multiply(solute_conc, cell_volumes)  # water per cell [cm3]
    net_mass_flux = new_solutes - solutes
    for k, root_flux in soil_mass_fluxes.items():
        net_mass_flux[k] -= root_flux * dt

    soil_water = new_soil_water
    solutes = new_solutes

    wall_soil_model = timeit.default_timer() - wall_soil_model

    wall_iteration = timeit.default_timer() - wall_iteration

    if i % skip == 0:
        collar_sx.append(s.getSolutionHeadAt(cci))
        min_sx.append(np.min(s.getSolutionHead()))
        water_cyl.append(np.sum(rs.get_water_volume()))  # cm3

        if rank == 0:
            write_file_array("solute_conc",solute_conc)
            write_file_array("water_content",water_content)
            out_times.append(t)
            collar_flux.append(r.collar_flux(rs_age + t, rx, rsx, soil_k, False))
            min_rsx.append(np.min(np.array(rsx)))
            min_rx.append(np.min(np.array(rx)))
            print("Cylindrical model: minimum root soil interface {:g} cm, soil {:g} cm, root xylem {:g} cm".format(min_rsx[-1], min_sx[-1], min_rx[-1]))
            min_soil_fluxes, max_soil_fluxes, summed_soil_fluxes = 1.e9, -1.e9, 0.
            for k, v in soil_fluxes.items():
                summed_soil_fluxes += v
                if max_soil_fluxes < v:
                    max_soil_fluxes = v
                if min_soil_fluxes > v:
                    min_soil_fluxes = v
            print("Fluxes: summed local fluxes {:g}, collar flux {:g}, predescribed {:g}".format(summed_soil_fluxes, collar_flux[-1], -trans * sinusoidal(t)))
            water_domain.append(np.min(soil_water))  # cm3
            water_collar_cell.append(soil_water[cci])  # cm3
            water_uptake.append(summed_soil_fluxes)  # cm3/day
            summed_mass_fluxes = 0.
            for k, v in soil_mass_fluxes.items():
                summed_mass_fluxes += v
            solute_uptake.append(summed_mass_fluxes)
            n = round(float(i) / float(NT) * 100.)
            print("[" + ''.join(["*"]) * n + ''.join([" "]) * (100 - n) + "], {:g} days".format(s.simTime))
            print("Iteration {:g} took {:g} seconds [{:g}% root, {:g}% rhizo {:g}% soil ]\n".
                  format(i, wall_iteration, wall_root_model / wall_iteration, wall_rhizo_models / wall_iteration, wall_soil_model / wall_iteration))
#             if min_rsx[-1] < -16000:
#                 print("breaksim time ", sim_time)
#                 break

""" plots and output """
if rank == 0:
    print ("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s")
    vp.plot_roots_and_soil(r.rs, "pressure head", rsx, s, periodic, min_b, max_b, cell_number, name, 1)  # VTK vizualisation
#     rsx = rs.get_inner_heads()  # matric potential at the root soil interface, i.e. inner values of the cylindric models [cm]
#     soil_k = np.divide(vg.hydraulic_conductivity(rsx, soil), rs.radii)  # only valid for homogenous soil
#     rx = r.solve(rs_age + t, -trans * sinusoidal(t), 0., rsx, False, wilting_point, soil_k)
#     fluxes = r.segFluxes(rs_age + t, rx, rsx, approx=False, cells=False, soil_k=soil_k)
#     ana = pb.SegmentAnalyser(r.rs)
#     ana.addData("rsx", rsx)
#     ana.addData("rx", rx)
#     ana.addData("fluxes", fluxes)
#     vp.plot_roots(ana, "rsx")  # VTK vizualisation
#     vp.plot_roots(ana, "fluxes")  # VTK vizualisation
#     crit_i = np.argmin(rsx)
#     print("critical segment", crit_i)
#     cidx = rs.seg2cell[crit_i]
#     print("mapped to cell",)
#     print("cell water content", water_content[cidx], "matric potential", s.getSolutionHeadAt(cidx))
#     rs.plot_cylinder(crit_i)
    vp.write_soil("mai", s, min_b, max_b, cell_number, ["phosphate"])

    fig, (ax1) = plt.subplots(1, 1)
    ax1.set_title("phosphate uptake")
    uptake = -np.array(solute_uptake)
    ax1.plot(out_times, uptake, label = "phosphate uptake")
    ax2 = ax1.twinx()
    ax2.plot(out_times, np.cumsum(uptake), 'c--', label = "cumulative uptake")
    ax1.legend()
    ax1.set_xlabel("Time (days)")
    ax1.set_ylabel("g/cm3")
    plt.show()

    plot_transpiration(out_times, water_uptake, collar_flux, lambda t: trans * sinusoidal(t))  # in rhizo_models.py
    plot_info(out_times, water_collar_cell, water_cyl, collar_sx, min_sx, min_rx, min_rsx, water_uptake, water_domain)  # in rhizo_models.py
    np.savetxt(name, np.vstack((out_times, -np.array(collar_flux), -np.array(water_uptake))), delimiter = ';')

