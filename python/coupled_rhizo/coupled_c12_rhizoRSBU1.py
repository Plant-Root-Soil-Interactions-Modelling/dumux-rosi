import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src")

import plantbox as pb  # CPlantBox
from rosi_richards import RichardsSP  # C++ part (Dumux binding), macroscopic soil model
from richards import RichardsWrapper  # Python part, macroscopic soil model
from functional.xylem_flux import *  # root system Python hybrid solver
from rhizo_modelsBU1 import *  # Helper class for cylindrical rhizosphere models

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

complicated MPI support (a non-mpi version of richards_cyl is needed, see script dumux3_nompi.sh)
"""

""" 
Parameters  
"""
name = "c12_rhizo_1cm"  # scenario name, to save results

""" soil """
min_b = [-4., -4., -15.]  # cm
max_b = [4., 4., 0.]  # cm
cell_number = [7, 7, 15]  # [8, 8, 15]  # [16, 16, 30]  # [32, 32, 60]  # [8, 8, 15] # [1]
periodic = False
static = True
static_rs = False
if static_rs:
    fname = "../../grids/RootSystem8.rsml"
else:
    fname = "Lupinus_albus_Leitner_2014_nods"
    if rank == 0:
        path = "../../../CPlantBox/modelparameter/structural/rootsystem/"
        
        rs = pb.MappedRootSystem()
        rs.setSeed(0)
        rs.readParameters(path + fname + ".xml")
        rs.setGeometry(pb.SDF_PlantBox( max_b[0]-min_b[0],  max_b[1]-min_b[1], max_b[2]-min_b[2]))  # to not let roots grow out of soil
        rs.initialize()
        rs.simulate(0.1, True)
        rs.write(fname + ".rsml");
    else:
        fname = fname + ".rsml"

domain_volume = np.prod(np.array(max_b) - np.array(min_b))

sand = [0.045, 0.43, 0.15, 3, 1000]
loam = [0.08, 0.43, 0.04, 1.6, 50]
# loam = [0.03, 0.345, 0.01, 2.5, 28.6]  # mai
clay = [0.1, 0.4, 0.01, 1.1, 10]
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
mode = "dumux"  # or "dumux_exact"
NC = 10  # dof+1
logbase = 1.5  # according to Mai et al. (2019)
split_type = 0  # type 0 == volume, type 1 == surface, type 2 == length
cyl_indxs = []

""" simulation time """
sim_time = 0.01  # 0.65  # 0.25  # [day]
dt = 30 / (24 * 3600)  # time step [day], 120 schwankt stark
NT = int(np.ceil(sim_time / dt))  # number of iterations
skip = 1  # for output and results, skip iteration

""" 
Initialize macroscopic soil model (Dumux binding)
"""
s = RichardsWrapper(RichardsSP())
s.initialize()
s.createGrid(min_b, max_b, cell_number, periodic)  # [cm]
s.setHomogeneousIC(initial, True)  # cm pressure head, equilibrium
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
print("len(rs.segments)Before_setRectangularGrid",len(rs.segments))
print(np.array([np.array(ns) for ns in rs.nodes]))
print(np.array([np.array(ns) for ns in rs.segments]))
picker = lambda x, y, z: s.pick([x, y, z])  #  function that return the index of a given position in the soil grid (should work for any grid - needs testing)

rs.setSoilGrid(picker,pb.Vector3d(min_b[0], min_b[1], min_b[2]), 
                        pb.Vector3d(max_b[0], max_b[1], max_b[2]),
                        pb.Vector3d(cell_number[0], cell_number[1], cell_number[2]), True)
print(np.array([np.array(ns) for ns in rs.nodes]))
print(np.array([np.array(ns) for ns in rs.segments]))
print("len(rs.segments)After_setRectangularGrid",len(rs.segments))
print(np.array([min_b[0], min_b[1], min_b[2]]), 
                        np.array([max_b[0], max_b[1], max_b[2]]),
                       np.array([cell_number[0], cell_number[1], cell_number[2]]))
#raise Exception
if rank == 0:
    rs = RhizoMappedSegments(rs, wilting_point, NC, logbase, mode)
else:
    rs = RhizoMappedSegments(fname, wilting_point, NC, logbase, mode)
    
r = XylemFluxPython(rs)  # wrap the xylem model around the MappedSegments
init_conductivities_growth(r, age_dependent, 0.05)  # age_dependent is a boolean, root conductivies are given in the file src/python_modules/root_conductivities.py

#print("rs.setSoilGrid(picker) ",rank)
print("len(rs.segments)Before",len(rs.segments))
#rs.setSoilGrid(picker)  # maps segments, maps root segements and soil grid indices to each other in both directions


print("len(rs.segments)After",len(rs.segments))
print("rs.setSoilGrid(picker) ",rank)
rs.set_xylem_flux(r)
# For debugging
# r.plot_conductivities()
r.test()
""" 
Initialize local soil models (around each root segment) 
"""
start_time = timeit.default_timer()
if rank == 0:
    x = s.getSolutionHead()[:, 0]  # initial condition of soil [cm]
# y the "[:, 0]" ???
x = comm.bcast(x, root = 0)  # Soil part runs parallel
ns = len(rs.segments)
dcyl = int(np.floor(ns / max_rank))

repartition = np.array([dcyl for i in range(max_rank)])
repartition[max_rank -1] = ns - rank * dcyl
assert np.sum(repartition) == ns
if rank == 0:
    rs.initialize(soil_, x, np.array([i for i in range(repartition[rank])]))
    print ("Initialized rank {:g}/{:g} [{:g}-{:g}] in {:g} s".format(rank + 1, max_rank, 0,repartition[rank], timeit.default_timer() - start_time))
else:
    rs.initialize(soil_, x, np.array([i for i in range(repartition[rank-1],repartition[rank])]))
    print ("Initialized rank {:g}/{:g} [{:g}-{:g}] in {:g} s".format(rank + 1, max_rank, repartition[rank-1],repartition[rank], timeit.default_timer() - start_time))

#print("press any key"); input()

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
cci = picker(rs.nodes[0].x, rs.nodes[0].y, rs.nodes[0].z)  # collar cell index
cell_volumes = s.getCellVolumes()  # cm3
cell_volumes = comm.bcast(cell_volumes, root = 0)
net_flux = np.zeros(cell_volumes.shape)

for i in range(0, NT):

    """ 0. plant growth """
    # simulate
    nsOld = ns
    dcylOld = dcyl
    repartitionOld = repartition
    
    ns = len(rs.segments)
    dcyl = int(np.floor(ns / max_rank))
    repartition = np.array([dcyl for i in range(max_rank)])
    repartition[max_rank -1] = ns - rank * dcyl
    toAdd = repartition - repartitionOld
    if toAdd[rank] > 0:
        print(rank, toAdd, toAdd[rank])
        rs.update(soil_, x, range(toAdd[rank-1],toAdd[rank]) + nsOld )
        print ("Initialized rank {:g}/{:g} [{:g}-{:g}] in {:g} s".format(rank + 1, max_rank, toAdd[rank-1] + nsOld, toAdd[rank] + nsOld, timeit.default_timer() - start_time))
        
    
    #print("press any key"); input()
    """ flow """
    
    if rank == 0:
        x = s.getSolutionHead() # initial condition of soil [cm]
    x = comm.bcast(x, root = 0)  # Soil part runs parallel
    
    wall_iteration = timeit.default_timer()

    t = i * dt  # current simulation time

    """ 1. xylem model """
    wall_root_model = timeit.default_timer()
    rsx = rs.get_inner_heads()  # matric potential at the root soil interface, i.e. inner values of the cylindric models [cm]
    print("inner_head",rsx)
    soil_k = np.divide(vg.hydraulic_conductivity(rsx, soil), rs.radii)  # only valid for homogenous soil
    if rank == 0:
        rx = r.solve(rs_age + t, -trans * sinusoidal(t), 0., rsx, False, wilting_point, soil_k)  # [cm]   False means that rsx is given per root segment not per soil cell
        proposed_inner_fluxes = r.segFluxes(rs_age + t, rx, rsx, approx = False, cells = False, soil_k = soil_k)  # [cm3/day]
    else:
        proposed_inner_fluxes = None
        rx = None
    print("xylem",rx,proposed_inner_fluxes)
    if rank == 0:
        #r.collar_flux(rs_age + t, rx, rsx, soil_k, False))
        collar_seg_id = np.where(np.array([ns.x for ns in rs.segments]) == 0)[0][0]
        #collar_flux_ = r.collar_flux(rs_age + t, rx.copy(), rsx.copy(), k_soil = soil_k, cells = False)  # validity checks
        collar_flux_= r.axial_flux(collar_seg_id, rs_age + t, rx.copy(), rsx.copy(), k_soil = soil_k, cells = False)
        #collar_flux_3= r.axial_fluxes(rs_age + t, rx.copy(), rsx.copy(), k_soil = soil_k, cells = False)
        err = np.linalg.norm(np.abs(np.sum(proposed_inner_fluxes) )- np.abs(collar_flux_))
        
        if err > 1.e-6:
            #print(np.abs(rsx1-rsx2))
            print(np.array([np.array(ns) for ns in rs.nodes]))
            print(np.array([np.array(ns) for ns in rs.segments]))
            #print(collar_flux_3,collar_seg_id)
            print("error: summed root surface fluxes and root collar flux differ" , err,  collar_flux_, np.sum(proposed_inner_fluxes))
            
            raise Exception
        err2 = np.linalg.norm(np.abs(-trans * sinusoidal(t)) -np.abs( collar_flux_))
        #print(seg_fluxes)
        print("errors: proposed_inner_fluxes",np.sum(proposed_inner_fluxes),"collar",collar_flux_,"transpiration",-trans * sinusoidal(t),err,err2)
        if r.last == "neumann":
            if err2 > 1.e-6:
                print("error: potential transpiration differs root collar flux in Neumann case" , err2)
                raise Exception
    #print("success: summed root surface fluxes and root collar flux agree" , err, collar_flux_, np.sum(proposed_inner_fluxes))
    #raise Exception()
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

            """ Additional sink plot """
            if i % (60 * 12) == 0:  # every 6h
                ana = pb.SegmentAnalyser(r.rs)
                fluxes = r.segFluxes(rs_age + t, rx, rsx, approx = False, cells = False, soil_k = soil_k)  # cm3/day
                ana.addData("fluxes", fluxes)  # cut off for vizualisation
                ana.addData("fluxes2", realized_inner_fluxes)  # cut off for vizualisation
                flux1d = ana.distribution("fluxes", max_b[2], min_b[2], 15, False)
                flux1d2 = ana.distribution("fluxes2", max_b[2], min_b[2], 15, False)
                sink1d.append(np.array(flux1d))
                sink1d2.append(np.array(flux1d2))
                # realized_inner_fluxes!!!!!!!!!!!!

""" plots and output """
if rank == 0:
    print ("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s")

    # vp.plot_roots_and_soil(r.rs, "pressure head", rsx, s, periodic, min_b, max_b, cell_number, name)  # VTK vizualisation

# #     rsx = rs.get_inner_heads()  # matric potential at the root soil interface, i.e. inner values of the cylindric models [cm]
# #     soil_k = np.divide(vg.hydraulic_conductivity(rsx, soil), rs.radii)  # only valid for homogenous soil
# #     rx = r.solve(rs_age + t, -trans * sinusoidal(t), 0., rsx, False, wilting_point, soil_k)
# #     fluxes = r.segFluxes(rs_age + t, rx, rsx, approx=False, cells=False, soil_k=soil_k)
# #     ana = pb.SegmentAnalyser(r.rs)
# #     ana.addData("rsx", rsx)
# #     ana.addData("rx", rx)
# #     ana.addData("fluxes", fluxes)
# #     vp.plot_roots(ana, "rsx")  # VTK vizualisation
# #     vp.plot_roots(ana, "fluxes")  # VTK vizualisation

    # crit_min_i, crit_max_i, crit_min_o, crit_max_o = rs.plot_cylinders()
    # # print(crit_min_i)
    # rs.plot_cylinder(crit_min_i)
    # print(rs.radii[crit_min_i])

    # plot_transpiration(out_times, water_uptake, collar_flux, lambda t: trans * sinusoidal(t))  # in rhizo_models.py
    # # plot_info(out_times, water_collar_cell, water_cyl, collar_sx, min_sx, min_rx, min_rsx, water_uptake, water_domain)  # in rhizo_models.py

    # np.savetxt(name, np.vstack((out_times, -np.array(collar_flux), -np.array(water_uptake))), delimiter = ';')

    # sink1d = np.array(sink1d)
    # np.save("sink1d_rhizo", sink1d)
    # np.save("sink1d2_rhizo", sink1d2)

    # print(sink1d.shape)
