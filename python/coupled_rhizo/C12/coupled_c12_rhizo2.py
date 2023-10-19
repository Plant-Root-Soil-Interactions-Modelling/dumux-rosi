import sys; sys.path.append("../../modules"); sys.path.append("../../../build-cmake/cpp/python_binding/");
sys.path.append("../../../../CPlantBox");  sys.path.append("../../../../CPlantBox/src");

import plantbox as pb  # CPlantBox
from rosi_richards import RichardsSP  # C++ part (Dumux binding), macroscopic soil model
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

from scenario_setup import *

""" 
Benchmark M1.2 static root system in soil

Coupled to cylindrical rhizosphere models using 1d richards equation (DUMUX solver)

Method description:
We use Neumann boundary condition, and prescribe the radial flux that the root hydraulic model suggests. 
This method requires a very small time step. Additionally, for regularization we used an effective radial conductivity 
$r_{r,eff} := min(k_r, k_s)$ for the root hydraulic model, where $k_s$ is the soil conductivity at the root soil interface. 

complicated MPI support (a non-mpi version of richards_cyl is needed, see script dumux3_nompi.sh)
"""

""" 
Parameters  
"""

r, rs, rs_age, trans, wilting_point, soil_, s, sra_table_lookup, mapping, sim_time, dt, skip = set_scenario("3D")
name = "c12_rhizo_1cm"  # scenario name, to save results
split_type = 0  # type 0 == volume, type 1 == surface, type 2 == length
dt = 30 / (24 * 3600)  # time step [day], 120 schwankt stark
NT = int(np.ceil(sim_time / dt))  # number of iterations
soil = vg.Parameters(soil_)

""" 
Initialize local soil models (around each root segment) 
"""
print("types")
print(type(r))
print(type(r.rs))
print(type(rs))
print("rs_age", rs_age)

rs = r.rs
start_time = timeit.default_timer()
x = s.getSolutionHead()[:]  # initial condition of soil [cm]
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
water_uptake, water_collar_cell, water_cyl, water_domain = [], [], [], []  # cm3
sink1d = []
sink1d2 = []
out_times = []  # days
cci = s.pick([rs.nodes[0].x, rs.nodes[0].y, rs.nodes[0].z])  # collar cell index
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
    if rs.mode == "dumux":
        proposed_inner_fluxes = comm.bcast(proposed_inner_fluxes, root = 0)
        rs.solve(dt, proposed_inner_fluxes, proposed_outer_fluxes)
    elif rs.mode == "dumux_exact":
        rx = comm.bcast(rx, root = 0)
        soil_k = comm.bcast(soil_k, root = 0)
        rs.solve(dt, rx, proposed_outer_fluxes, rsx, soil_k)
    realized_inner_fluxes = rs.get_inner_fluxes()
    realized_inner_fluxes = comm.bcast(realized_inner_fluxes, root = 0)
    wall_rhizo_models = timeit.default_timer() - wall_rhizo_models

    """ 3a. macroscopic soil model """
    wall_soil_model = timeit.default_timer()
    water_content = np.array(s.getWaterContent())
    water_content = comm.bcast(water_content, root = 0)
    soil_water = np.multiply(water_content, cell_volumes)  # water per cell [cm3]
    soil_fluxes = r.sumSegFluxes(realized_inner_fluxes)  # [cm3/day]  per soil cell
    s.setSource(soil_fluxes.copy())  # [cm3/day], in richards.py
    s.solve(dt)  # in solverbase.py

    """ 3b. calculate net fluxes """
    water_content = np.array(s.getWaterContent())
    water_content = comm.bcast(water_content, root = 0)
    new_soil_water = np.multiply(water_content, cell_volumes)  # calculate net flux
    net_flux = new_soil_water - soil_water  # change in water per cell [cm3]
    for k, root_flux in soil_fluxes.items():
        net_flux[k] -= root_flux * dt
    soil_water = new_soil_water
    wall_soil_model = timeit.default_timer() - wall_soil_model

    wall_iteration = timeit.default_timer() - wall_iteration

    if i % skip == 0:
        collar_sx.append(s.getSolutionHeadAt(cci))
        min_sx.append(np.min(s.getSolutionHead()))
        water_cyl.append(np.sum(rs.get_water_volume()))  # cm3

        if rank == 0:
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
            n = round(float(i) / float(NT) * 100.)
            print("[" + ''.join(["*"]) * n + ''.join([" "]) * (100 - n) + "], {:g} days".format(s.simTime))
            print("Iteration {:g} took {:g} seconds [{:g}% root, {:g}% rhizo {:g}% soil ]\n".
                  format(i, wall_iteration, wall_root_model / wall_iteration, wall_rhizo_models / wall_iteration, wall_soil_model / wall_iteration))
#             if min_rsx[-1] < -16000:
#                 print("breaksim time ", sim_time)
#                 break

            # """ Additional sink plot """
            # if i % (60 * 12) == 0:  # every 6h
            #     ana = pb.SegmentAnalyser(r.rs)
            #     fluxes = r.segFluxes(rs_age + t, rx, rsx, approx = False, cells = False, soil_k = soil_k)  # cm3/day
            #     ana.addData("fluxes", fluxes)  # cut off for vizualisation
            #     ana.addData("fluxes2", realized_inner_fluxes)  # cut off for vizualisation
            #     # flux1d = ana.distribution("fluxes", max_b[2], min_b[2], 15, False)
            #     # flux1d2 = ana.distribution("fluxes2", max_b[2], min_b[2], 15, False)
            #     # sink1d.append(np.array(flux1d))
            #     # sink1d2.append(np.array(flux1d2))
            #     # realized_inner_fluxes!!!!!!!!!!!!

""" plots and output """
if rank == 0:
    print ("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s")

    vp.plot_roots_and_soil(r.rs, "pressure head", rsx, s, periodic, min_b, max_b, cell_number, name)  # VTK vizualisation

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

    crit_min_i, crit_max_i, crit_min_o, crit_max_o = rs.plot_cylinders()
    # print(crit_min_i)
    rs.plot_cylinder(crit_min_i)
    print(rs.radii[crit_min_i])

    plot_transpiration(out_times, water_uptake, collar_flux, lambda t: trans * sinusoidal(t))  # in rhizo_models.py
    # plot_info(out_times, water_collar_cell, water_cyl, collar_sx, min_sx, min_rx, min_rsx, water_uptake, water_domain)  # in rhizo_models.py

    np.savetxt(name, np.vstack((out_times, -np.array(collar_flux), -np.array(water_uptake))), delimiter = ';')

    sink1d = np.array(sink1d)
    np.save("sink1d_rhizo", sink1d)
    np.save("sink1d2_rhizo", sink1d2)

    print(sink1d.shape)
