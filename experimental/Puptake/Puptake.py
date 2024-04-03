import sys; sys.path.append("../../python/modules"); 
sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src")

import os
import plantbox as pb  # CPlantBox
from rosi_richardsnc import RichardsNCSP
from richards import RichardsWrapper  # Python part, macroscopic soil model
from functional.xylem_flux import *  # root system Python hybrid solver
from rhizo_models_Pupt import *  # Helper class for cylindrical rhizosphere models

import visualisation.vtk_plot as vp
import functional.van_genuchten as vg
from functional.root_conductivities import *

import numpy as np
import timeit
import matplotlib.pyplot as plt
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank(); max_rank = comm.Get_size()


def write_file_array(name, data, space =","):
    name2 = './results/Puptake/'+ name+ '.txt'
    with open(name2, 'a') as log:
        log.write(space.join([num for num in map(str, data)])  +'\n')
results_dir="./results/Puptake/"
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
troubleShooting = True

""" soil """
min_b = [-4., -4., -15.]  # cm
max_b = [4., 4., 0.]  # cm
domain_volume = np.prod(np.array(max_b) - np.array(min_b))
cell_number = [7, 7, 15]  # [8, 8, 15]  # [16, 16, 30]  # [32, 32, 60]  # [8, 8, 15] # [1]
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
mode = "dumux_nc"  #
NC = 10  # dof+1
logbase = 1.5  # according to Mai et al. (2019)
split_type = 0  # type 0 == volume, type 1 == surface, type 2 == length

""" simulation time """
sim_time = 10/(24*60)  # 0.65  # 0.25  # [day]
dt =1/ (24 * 60)  #  1/24/60#30 / (24 * 3600)  # time step [day]
NT = int(np.ceil(sim_time / dt))  # number of iterations
skip = 1  # for output and results, skip iteration

""" 
Initialize macroscopic soil model (Dumux binding)
"""
s = RichardsWrapper(RichardsNCSP())
s.initialize()
s.results_dir = results_dir
#### soil shape and resolution
s.createGrid(min_b, max_b, cell_number, periodic)  # [cm]

#### soil biochemical parameters + solver parameter

# ATT! if you change parmeters here, need to do the same changes for the parameters
# of the 1d models (in rhizo_models_Pupt::initialize_dumux_nc_neumann() )
s.setParameter("Component.MolarMass", "3.1e-2") # kg/mol
s.setParameter("Component.LiquidDiffusionCoefficient", "6.e-10")  # m^2 s-1
s.setParameter("Component.freundlichN_","0")# "124.8.")# adsorption
s.setParameter("Component.freundlichK_", "0")#".4")# adsorption
s.setVGParameters([soil_])
s.setParameter("Newton.EnableAbsoluteResidualCriterion", "True")
s.setParameter("Newton.Verbosity","0")
s.setParameter("Newton.EnableChop", "True")# manually reset the results values within preset range 
######

##### IC and BC, no boundary flux (no rain/rvaporation/exchange with the deeper soil)
s.setHomogeneousIC(initial, True)  # cm pressure head, equilibrium
s.setTopBC("noFlux") # for water
s.setBotBC("noFlux") # for water

IC_C = 0.8 # g / g solution <== has to be < 1.
s.setParameter("Soil.IC.C",str(IC_C))  # g / g solution, = g/cm3 water 
s.setParameter("Soil.BC.Top.SType", "2")  
s.setParameter("Soil.BC.Top.CValue", "0.")  # 
s.setParameter("Soil.BC.Bot.SType", "2")  
s.setParameter("Soil.BC.Bot.CValue", "0.")


s.initializeProblem()
s.setCriticalPressure(wilting_point)  # to remain above wilting point
s.ddt = 1.e-5  # [day] initial Dumux time step

""" 
Initialize xylem model 
"""
rs = RhizoMappedSegments("../../grids/RootSystem8.rsml", wilting_point, NC, logbase, mode, s)
rs.setRectangularGrid(pb.Vector3d(min_b[0], min_b[1], min_b[2]), pb.Vector3d(max_b[0], max_b[1], max_b[2]),
                        pb.Vector3d(cell_number[0], cell_number[1], cell_number[2]), True)
# True: root segments are cut  to the soil grid so that each segment is completely within one soil control element, this works only for rectangular grids so far
picker = lambda x, y, z: s.pick([x, y, z])  #  function that return the index of a given position in the soil grid (should work for any grid - needs testing)
rs.setSoilGrid(picker)  # maps segments, maps root segements and soil grid indices to each other in both directions
r = XylemFluxPython(rs)  # wrap the xylem model around the MappedSegments
init_conductivities(r, age_dependent)  # age_dependent is a boolean, root conductivies are given in the file /root_conductivities.py
rs.set_xylem_flux(r)


""" 
Initialize local soil models (around each root segment) 
"""
start_time = timeit.default_timer()
x = s.getSolutionHead()  # initial condition of soil [cm]
x = comm.bcast(x, root = 0)  # Soil part runs parallel
cc_ = s.getSolution(1)# g/g or g/cm water? 

ns = len(rs.segments)
dcyl = int(np.floor(ns / max_rank))
if rank + 1 == max_rank:
    rs.initialize(soil_, x, np.array(range(rank * dcyl, ns)), cc = cc_)
    print ("Initialized rank {:g}/{:g} [{:g}-{:g}] in {:g} s".format(rank + 1, max_rank, rank * dcyl, ns, timeit.default_timer() - start_time))
else:
    rs.initialize(soil_, x, np.array(range(rank * dcyl, (rank + 1) * dcyl)), cc = cc_)
    print ("Initialized rank {:g}/{:g} [{:g}-{:g}] in {:g} s".format(rank + 1, max_rank, rank * dcyl, (rank + 1) * dcyl, timeit.default_timer() - start_time))
# print("press any key"); input()

if troubleShooting:
    rs.checkVolumeBalance()
    #print(rs.allDiff1d3dVol_rel,rs.allDiff1d3dVol_abs )


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
cci = picker(rs.nodes[0].x, rs.nodes[0].y, rs.nodes[0].z)  # collar cell index (correct?)
cell_volumes = s.getCellVolumes()  # cm3
cell_volumes = comm.bcast(cell_volumes, root = 0)
net_flux = np.zeros(cell_volumes.shape) # for water fluxes between the d soil voxels
net_mass_flux = net_flux # for solute fluxes between the d soil voxels



if troubleShooting:
    rs.checkMassOMoleBalance(sourceWat=net_flux, # cm3/day 
                                     sourceSol=net_mass_flux, # g/day
                                     dt=dt,        # day                                         
                                     diff1d3dCW_abs_lim = np.Inf,
                              verbose_ = False,
                              diff1d3dCW_rel_lim =np.full(10, np.Inf))
    print("rs.allDiff1d3dCW_rel",max(rs.allDiff1d3dCW_rel[0]),max(rs.allDiff1d3dCW_rel[1])) # initial relative error. to compare with error during simulation
    

for i in range(0, NT):

    wall_iteration = timeit.default_timer()

    t = i * dt  # current simulation time

    """ 1. xylem model """
    
    wall_root_model = timeit.default_timer()
    rsx = rs.get_inner_heads()  # matric potential at the root soil interface, i.e. inner values of the cylindric models [cm]
    deltaR_cm =  rs.getDeltaR() # 1/2 length of the 1d domains inner cells [cm]
    
    soil_k = np.divide(vg.hydraulic_conductivity(rsx, soil), deltaR_cm)  
    
    if rank == 0:
        rx = r.solve(rs_age + t, -trans * sinusoidal(t), 0., rsx, False, wilting_point, soil_k)  # [cm]   False means that rsx is given per root segment not per soil cell
        proposed_inner_fluxes = np.array(r.segFluxes(rs_age + t, rx, rsx, approx = False, cells = False, soil_k = soil_k))  # [cm3/day]
    else:
        proposed_inner_fluxes = None
        rx = None
    wall_root_model = timeit.default_timer() - wall_root_model

    """ 2. local soil models """
    wall_rhizo_models = timeit.default_timer()
    proposed_outer_fluxes = r.splitSoilFluxes(net_flux / dt, split_type)
    proposed_outer_fluxes_C = r.splitSoilFluxes(net_mass_flux / dt, split_type)
    
    if rs.mode == "dumux_nc":
        # inner water flux from plant water uptake
        proposed_inner_fluxes = comm.bcast(proposed_inner_fluxes, root = 0)
        # no need to send inner_flux_C: computed in dumux according to the
        # active uptake parameters
        
        rs.solve(dt, proposed_inner_fluxes, proposed_outer_fluxes, proposed_outer_fluxes_C)  # mass fluxes
    else:
        print(rs.mode)
        raise Exception("this script is for dumux_nc only")
    
    # flux rate at the inner and outer boundaries
    Q_inW, Q_outW, Q_inS, Q_outS  = rs.getSavedBCs()# cm3/d, g/d
    
    
    realized_inner_fluxes = rs.get_inner_fluxes() # Q_inW#get_inner_fluxes only works if time step is small enough
    realized_inner_fluxes = comm.bcast(realized_inner_fluxes, root = 0)
    
    realized_mass_fluxes =  rs.get_inner_mass_fluxes() #Q_inS# get_inner_mass_fluxes only works if time step is small enough
    realized_mass_fluxes = comm.bcast(realized_mass_fluxes, root = 0)
    wall_rhizo_models = timeit.default_timer() - wall_rhizo_models
    """ 3a. macroscopic soil model """
    wall_soil_model = timeit.default_timer()
    water_content = np.array(s.getWaterContent())
    water_content = comm.bcast(water_content, root = 0)
    soil_water = np.multiply(water_content, cell_volumes)
    solute_conc = np.array(s.getSolution(1))  # solute concentration
    solutes = np.multiply(solute_conc, soil_water)  # water per cell [cm3]
    soil_fluxes = r.sumSegFluxes(realized_inner_fluxes)  # [cm3/day]  per soil cell
    s.setSource(soil_fluxes.copy())  # [cm3/day], in richards.py
    soil_mass_fluxes = r.sumSegFluxes(realized_mass_fluxes)  # [cm3/day]  per soil cell
    s.setSource(soil_mass_fluxes.copy(), 1)  # [g/day], in richards.py
    
    """ 3b. check mass balance """
    
    if troubleShooting:
        rs.checkMassOMoleBalance(sourceWat=soil_fluxes, # cm3/day 
                                         sourceSol=soil_mass_fluxes, # g/day
                                         dt=dt,        # day                                         
                                         diff1d3dCW_abs_lim = np.Inf,
                                  verbose_ = False,
                                  diff1d3dCW_rel_lim =np.full(10, np.Inf))
        print("rs.allDiff1d3dCW_rel",max(rs.allDiff1d3dCW_rel[0]),max(rs.allDiff1d3dCW_rel[1])) #relative error between 1d and 3d
        # shape = (num component, num soil cells)

    """ 3c. solve """
    s.solve(dt)  # in solverbase.py

    """ 3d. calculate net fluxes """
    water_content = np.array(s.getWaterContent())
    water_content = comm.bcast(water_content, root = 0)
    new_soil_water = np.multiply(water_content, cell_volumes)  # calculate net flux
    net_flux = new_soil_water - soil_water  # change in water per cell [cm3]
    # print('net_flux',net_flux,'soil_fluxes',soil_fluxes[0]*dt, sum(realized_inner_fluxes)*dt)
    # print(net_flux - sum(realized_inner_fluxes)*dt)
    # print('avgd', s.base.getAvgDensity())
    for k, root_flux in soil_fluxes.items():
        net_flux[k] -= root_flux * dt

    solute_conc = np.array(s.getSolution(1))  # solute concentration, g/g = g/cm3
    solute_conc = comm.bcast(solute_conc, root = 0)
    new_solutes = np.multiply(solute_conc, new_soil_water)  # water per cell [cm3]
    net_mass_flux = new_solutes - solutes
    
    if (rs.allDiff1d3dCW_rel[1] > 1e-10).any():
        raise Exception
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

