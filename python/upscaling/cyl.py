"""
functions for the cylindrical coupling approach
"""
import numpy as np
import timeit
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()

from xylem_flux import sinusoidal2

def simulate_const(s, rs, trans, sim_time, dt):
    """     
    simulates the coupled scenario       
        root architecture is not gowing  
        conductivities are not changing over time
        
    s
    rs
    sra_table_lookup
    trans
    sim_time    
    dt
    
    TODO recyle factorisation of left hand side ... 
    """
    wilting_point = -15000 # cm    
    skip = 3  # for output and results, skip iteration
    rs_age = 0.
    split_type = 0  # type 0 == volume, type 1 == surface, type 2 == length
    
    """ 
    Initialize local soil models (around each root segment) 
    """
    start_time = timeit.default_timer()
    x = s.getSolutionHead()  # initial condition of soil [cm]
    rs.initialize(s.soils[0], x)
    r = rs.rs # rename (XylemFluxPython)
    
    """ 
    Simulation 
    
    loop
    1. xylem model
    2. local soil models
    3. macroscopic soil model 
    """
    print("Starting simulation")
    start_time = timeit.default_timer()
    
    psi_x_, psi_s_, sink_ , x_, y_, psi_s2_ = [], [], [], [], [], []  # for post processing
    
    water_uptake, water_collar_cell, water_cyl, water_domain = [], [], [], []  # cm3
    
    cell_volumes = s.getCellVolumes()  # cm3
    cell_volumes = comm.bcast(cell_volumes, root = 0)
    net_flux = np.zeros(cell_volumes.shape)
    
    N = int(np.ceil(sim_time / dt))  # number of iterations

    rsx = rs.get_inner_heads()  # matric potential at the root soil interface, i.e. inner values of the cylindric models (not extrapolation to the interface!) [cm]
    rx = r.solve(rs_age, -trans * sinusoidal2(0, dt), 0., rsx, cells = False, wilting_point = wilting_point, soil_k = [])

    for i in range(0, N):
    
        t = i * dt  # current simulation time
    
        """ 1. xylem model """
        rsx = rs.get_inner_heads()  # matric potential at the root soil interface, i.e. inner values of the cylindric models (not extrapolation to the interface!) [cm]
        # soil_k = np.divide(vg.hydraulic_conductivity(rsx, soil), rs.radii)  # only valid for homogenous soil
        soil_k = rs.get_soil_k(rx)
        rx = r.solve(rs_age + t, -trans * sinusoidal2(t, dt), 0., rsx, cells = False, wilting_point = wilting_point, soil_k = soil_k)
    
        # validity check
        proposed_inner_fluxes = r.segFluxes(rs_age + t, rx.copy(), rsx.copy(), approx = False, cells = False, soil_k = soil_k.copy())  # [cm3/day]
        collar_flux = r.collar_flux(rs_age + t, rx.copy(), rsx.copy(), k_soil = soil_k.copy(), cells = False)
        err = np.linalg.norm(np.sum(proposed_inner_fluxes) - collar_flux)
        if err > 1.e-10:
            print("error: summed root surface fluxes and root collar flux differ" , err)
    
        """ 2. local soil models """
        proposed_outer_fluxes = r.splitSoilFluxes(net_flux / dt, split_type)
        rs.solve(dt, proposed_inner_fluxes, proposed_outer_fluxes)  # left and right neumann fluxes
        realized_inner_fluxes = rs.get_inner_fluxes()  # identical for mode = "dumux"
    
        # validity check
        err = np.linalg.norm(np.array(proposed_inner_fluxes) - np.array(realized_inner_fluxes))
        if err > 1.e-15:
            print("error: summed root surface fluxes and cylindric model fluxes differ" , err)
    
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
            # net_flux[k] = 0.  # same as SRA approach (currently using a no-flux at boundary)
        soil_water = new_soil_water
    
        """ remember results ... """
        if i % skip == 0:
            psi_x_.append(rx.copy())  # cm
            psi_s_.append(rsx.copy()) # cm
            sink_.append(proposed_inner_fluxes.copy()) # cm3/day
            x_.append(t) # day
            y_.append(np.sum(proposed_inner_fluxes)) # cm3/day 
            psi_s2_.append(s.getSolutionHead()[:, 0])
            print("{:g}/{:g} iterations".format(i, N)) 
            
            
    print ("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s")
    
    return psi_x_, psi_s_, sink_, x_, y_, psi_s2_          


  