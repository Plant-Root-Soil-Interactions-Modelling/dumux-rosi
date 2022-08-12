"""
functions for the cylindrical coupling approach
"""
import numpy as np
import timeit
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank(); max_rank = comm.Get_size()

from xylem_flux import sinusoidal2


def simulate_const(s, rs, trans, sim_time, dt):
    """     
    simulates the coupled scenario       
        root architecture is not gowing  
        conductivities are not changing over time
        
    s
    rs
    trans
    sim_time    
    dt
    
    TODO recyle factorisation of left hand side ... 
    """
    wilting_point = -15000  # cm
    skip = 3 * 6  # for output and results, skip iteration
    rs_age = 0.  # day
    split_type = 0  # type 0 == volume, type 1 == surface, type 2 == length

    """ 
    Initialize local soil models (around each root segment) 
    """
    start_time = timeit.default_timer()
    sx = s.getSolutionHead()  # initial condition of soil [cm]
    sx = comm.bcast(sx, root = 0)  # Soil part runs parallel
    ns = len(rs.segments)
    dcyl = int(np.floor(ns / max_rank))
    if rank + 1 == max_rank:
        print ("Initialized final rank {:g}/{:g} [{:g}-{:g}] in {:g} s".format(rank + 1, max_rank, rank * dcyl, ns, timeit.default_timer() - start_time))
        rs.initialize(s.soils[0], sx, np.array(range(rank * dcyl, ns)))
    else:
        print ("Initialized rank {:g}/{:g} [{:g}-{:g}] in {:g} s".format(rank + 1, max_rank, rank * dcyl, (rank + 1) * dcyl, timeit.default_timer() - start_time))
        rs.initialize(s.soils[0], sx, np.array(range(rank * dcyl, (rank + 1) * dcyl)))

    r = rs.rs  # rename (XylemFluxPython)

    start_time = timeit.default_timer()

    psi_x_, psi_s_, sink_ , x_, y_, psi_s2_ = [], [], [], [], [], []  # for post processing

    cell_volumes = s.getCellVolumes()  # cm3
    cell_volumes = comm.bcast(cell_volumes, root = 0)
    net_flux = np.zeros(cell_volumes.shape)

    rsx = rs.get_inner_heads()  # matric potential at the root soil interface, i.e. inner values of the cylindric models (not extrapolation to the interface!) [cm]
    if rank == 0:
        rx = r.solve(rs_age, -trans * sinusoidal2(0, dt), 0., rsx, cells = False, wilting_point = wilting_point, soil_k = [])
    else:
        rx = None

    N = int(np.ceil(sim_time / dt))  # number of iterations

    """ simualtion loop """
    for i in range(0, N):

        t = i * dt  # current simulation time

        """ 1. xylem model """
        rsx = rs.get_inner_heads()  # matric potential at the root soil interface, i.e. inner values of the cylindric models (not extrapolation to the interface!) [cm]
        rx = comm.bcast(rx, root = 0)
        soil_k = rs.get_soil_k(rx)
        if rank == 0:
            # soil_k = np.divide(vg.hydraulic_conductivity(rsx, soil), rs.radii)  # only valid for homogenous soil
            rx = r.solve(rs_age, -trans * sinusoidal2(t, dt), 0., rsx, cells = False, wilting_point = wilting_point, soil_k = soil_k)
            proposed_inner_fluxes = r.segFluxes(rs_age, rx.copy(), rsx.copy(), approx = False, cells = False, soil_k = soil_k.copy())  # [cm3/day]
        else:
            proposed_inner_fluxes = None
            rx = None

        # validity check
        if rank == 0:
            collar_flux = r.collar_flux(rs_age, rx.copy(), rsx.copy(), k_soil = soil_k.copy(), cells = False)
            err = np.linalg.norm(np.sum(proposed_inner_fluxes) - collar_flux)
            if err > 1.e-8:
                print("error: summed root surface fluxes and root collar flux differ" , err)

        """ 2. local soil models """
        proposed_inner_fluxes = comm.bcast(proposed_inner_fluxes, root = 0)
        proposed_outer_fluxes = r.splitSoilFluxes(net_flux / dt, split_type)
        rs.solve(dt, proposed_inner_fluxes, proposed_outer_fluxes)  # left and right neumann fluxes
        realized_inner_fluxes = rs.get_inner_fluxes()  # identical for mode = "dumux"
        realized_inner_fluxes = comm.bcast(realized_inner_fluxes, root = 0)

        # validity check
        err = np.linalg.norm(np.array(proposed_inner_fluxes) - np.array(realized_inner_fluxes))
        if err > 1.e-8:
            print("error: summed root surface fluxes and cylindric model fluxes differ" , err)

        """ 3a. macroscopic soil model """
        water_content = np.array(s.getWaterContent())  # theta per cell [1]
        water_content = comm.bcast(water_content, root = 0)
        soil_water = np.multiply(water_content, cell_volumes)  # water per cell [cm3]
        soil_fluxes = r.sumSegFluxes(realized_inner_fluxes)  # [cm3/day]  per soil cell
        s.setSource(soil_fluxes.copy())  # [cm3/day], in moduels/richards.py
        s.solve(dt)  # in modules/solverbase.py

        """ 3b. calculate net fluxes """
        water_content = np.array(s.getWaterContent())
        water_content = comm.bcast(water_content, root = 0)
        new_soil_water = np.multiply(water_content, cell_volumes)  # calculate net flux
        net_flux = new_soil_water - soil_water  # change in water per cell [cm3]
        for k, root_flux in soil_fluxes.items():
            net_flux[k] -= root_flux * dt
            # net_flux[k] = 0.  # same as SRA approach (currently using a no-flux at boundary)
        soil_water = new_soil_water

        """ remember results ... """
        if i % skip == 0:
            sx = s.getSolutionHead()
            if rank == 0:
                sx = sx[:, 0]
                psi_x_.append(rx.copy())  # cm (per root node)
                psi_s_.append(rsx.copy())  # cm (per root segment)
                sink = np.zeros(sx.shape)
                for k, v in soil_fluxes.items():
                    sink[k] += v
                sink_.append(sink)  # cm3/day (per soil cell)
                x_.append(t)  # day
                y_.append(np.sum(sink))  # cm3/day
                psi_s2_.append(sx)  # cm (per soil cell)
                print("{:g}/{:g} iterations".format(i, N))

    if rank == 0:
        print ("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s")

    return psi_x_, psi_s_, sink_, x_, y_, psi_s2_

