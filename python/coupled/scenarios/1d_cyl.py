""" 
Jan's new scenario

Coupled to cylindrical rhizosphere models using 1d richards equation (DUMUX solver)

complicated MPI support (a non-mpi version of richards_cyl is needed, see script dumux3_nompi.sh)
"""
import sys; sys.path.append("../../../modules/"); sys.path.append("../../../../../CPlantBox/");  sys.path.append("../../../../../CPlantBox/src/python_modules")
sys.path.append("../../../../build-cmake/cpp/python_binding/"); sys.path.append("../../../modules/fv/");

from scenario_setup import *
from rhizo_models import *  # Helper class for cylindrical rhizosphere models

import timeit
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy import sparse
import scipy.sparse.linalg as LA
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank(); max_rank = comm.Get_size()

""" 
Initialize  
"""

name = "small_cyl"  # name to export resutls
sstr = "_dry"

r, rs_age, trans, wilting_point, soil_, s, sra_table_lookup, mapping, sim_time, dt, skip = set_scenario1D(sstr)
min_b, max_b, cell_number = get_domain1D()

""" rhizosphere models """
mode = "dumux"  # or "dumux_exact"
NC = 10  # dof+1
logbase = 1.5  # according to Mai et al. (2019)
split_type = 0  # type 0 == volume, type 1 == surface, type 2 == length

""" simulation time """
dt = 36 / (24 * 3600)  # time step [day], 120 schwankt stark
skip = 10  # for output and results, skip iteration

""" 
Initialize xylem model 
"""
fname = "../../../../grids/RootSystem_verysimple2.rsml"
rs = RhizoMappedSegments(fname, wilting_point, NC, logbase, mode)
rs.setRectangularGrid(pb.Vector3d(min_b[0], min_b[1], min_b[2]), pb.Vector3d(max_b[0], max_b[1], max_b[2]),
                        pb.Vector3d(cell_number[0], cell_number[1], cell_number[2]), True)
r.rs = rs
picker = lambda x, y, z: s.pick([0., 0., z])  #  function that return the index of a given position in the soil grid (should work for any grid - needs testing)
rs.setSoilGrid(picker)  # maps segments, maps root segements and soil grid indices to each other in both directions
rs.set_xylem_flux(r)

""" 
Initialize local soil models (around each root segment) 
"""
start_time = timeit.default_timer()
x = s.getSolutionHead()  # initial condition of soil [cm]
x = comm.bcast(x, root = 0)  # Soil part runs parallel
ns = len(rs.segments)
dcyl = int(np.floor(ns / max_rank))
if rank + 1 == max_rank:
    rs.initialize(soil_, x[:, 0], np.array(range(rank * dcyl, ns)))
    print ("Initialized rank {:g}/{:g} [{:g}-{:g}] in {:g} s".format(rank + 1, max_rank, rank * dcyl, ns, timeit.default_timer() - start_time))
else:
    rs.initialize(soil_, x[:, 0], np.array(range(rank * dcyl, (rank + 1) * dcyl)))
    print ("Initialized rank {:g}/{:g} [{:g}-{:g}] in {:g} s".format(rank + 1, max_rank, rank * dcyl, (rank + 1) * dcyl, timeit.default_timer() - start_time))

""" 
Simulation 

loop
1. xylem model
2. local soil models
3. macroscopic soil model 
"""
print("Starting simulation")
NT = int(np.ceil(sim_time / dt))  # number of iterations

start_time = timeit.default_timer()

# for post processing
water_domain = []  # cm3
cci = picker(rs.nodes[0].x, rs.nodes[0].y, rs.nodes[0].z)  # collar cell index
cell_volumes = s.getCellVolumes()  # cm3
cell_volumes = comm.bcast(cell_volumes, root = 0)
net_flux = np.zeros(cell_volumes.shape)

psi_x_, psi_s_, sink_, psi_s2_ = [], [], [], []  # for xls output
x_, y_ = [], []

for i in range(0, NT):

    wall_iteration = timeit.default_timer()

    t = i * dt  # current simulation time

    """ 1. xylem model """
    wall_root_model = timeit.default_timer()
    rsx = rs.get_inner_heads()  # matric potential at the root soil interface, i.e. inner values of the cylindric models [cm]
    # soil_k = np.divide(vg.hydraulic_conductivity(rsx, soil), rs.radii)  # only valid for homogenous soil
    if i > 0:
        soil_k = rs.get_soil_k(rx)
    else:
        soil_k = []

    if rank == 0:
        rx = r.solve(rs_age + t, -trans * sinusoidal(t), 0., rsx, False, wilting_point, soil_k)  # [cm]   False means that rsx is given per root segment not per soil cell
        proposed_inner_fluxes = r.segFluxes(rs_age + t, rx, rsx, approx = False, cells = False, soil_k = soil_k)  # [cm3/day]
    else:
        proposed_inner_fluxes = None
        rx = None
    rx = comm.bcast(rx, root = 0)  # needed by rs.get_soil_k
    wall_root_model = timeit.default_timer() - wall_root_model

    """ 2. local soil models """
    wall_rhizo_models = timeit.default_timer()
    proposed_outer_fluxes = r.splitSoilFluxes(net_flux / dt, split_type)
    if rs.mode == "dumux":
        proposed_inner_fluxes = comm.bcast(proposed_inner_fluxes, root = 0)
        rs.solve(dt, proposed_inner_fluxes, proposed_outer_fluxes)
    elif rs.mode == "dumux_exact":
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
        min_sx = np.min(s.getSolutionHead())
        dd = np.array(s.getSolutionHead())

        if rank == 0:
            x_.append(t)
            psi_x_.append(rx[1:])
            psi_s_.append(rsx)
            collar_flux = r.collar_flux(rs_age + t, rx, rsx, soil_k, False)

            min_soil_fluxes, max_soil_fluxes, summed_soil_fluxes = 1.e9, -1.e9, 0.
            sink = np.zeros(water_content[:, 0].shape)
            for k, v in soil_fluxes.items():
                sink[k] += v
                summed_soil_fluxes += v
                if max_soil_fluxes < v:
                    max_soil_fluxes = v
                if min_soil_fluxes > v:
                    min_soil_fluxes = v
            sink_.append(sink)
            y_.append(summed_soil_fluxes)  # cm3/day
            psi_s2_.append(dd[:, 0])
            water_domain.append(np.min(soil_water))  # cm3

            min_rsx = np.min(np.array(rsx))
            min_rx = np.min(np.array(rx))
            print("Cylindrical model: minimum root soil interface {:g} cm, soil {:g} cm, root xylem {:g} cm"
                  .format(min_rsx, min_sx, min_rx))
            print("Fluxes: summed local fluxes {:g}, collar flux {:g}, predescribed {:g}"
                  .format(summed_soil_fluxes, collar_flux, -trans * sinusoidal(t)))
            n = round(float(i) / float(NT) * 100.)
            print("[" + ''.join(["*"]) * n + ''.join([" "]) * (100 - n) + "], {:g} days".format(s.simTime))
            print("Iteration {:g} took {:g} seconds [{:g}% root, {:g}% rhizo {:g}% soil ]\n".
                  format(i, wall_iteration, wall_root_model / wall_iteration, wall_rhizo_models / wall_iteration, wall_soil_model / wall_iteration))

if rank == 0:
    print ("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s")
    file1 = 'results/psix_' + name + sstr  # per segment
    np.save(file1, np.array(psi_x_))
    file2 = 'results/psiinterface_' + name + sstr  # per segment
    np.save(file2, np.array(psi_s_))
    file3 = 'results/sink_' + name + sstr
    np.save(file3, np.array(-np.array(sink_)))
    file4 = 'results/transpiration_' + name + sstr
    np.save(file4, np.vstack((x_, -np.array(y_))))
    file5 = 'results/soil_' + name + sstr
    np.save(file5, np.array(psi_s2_))
    print("fin")

