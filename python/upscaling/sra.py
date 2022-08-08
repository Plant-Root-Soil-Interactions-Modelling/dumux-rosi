"""
functions for the steady rate approach
"""
import numpy as np
import timeit

from xylem_flux import sinusoidal2
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import fsolve

import van_genuchten as vg

def open_sra_lookup(filename):
    """ opens the look-up table from a file, to quickly find soil root interface potential """
    sra_table = np.load(filename + ".npy")
    x = np.load(filename + "_.npy", allow_pickle = True)
    kx_ = x[0]
    sx_ = x[1]
    inner_ = x[2]
    outer_ = x[3]
    # kr_ = x[4]
    # soil = x[5]
    # print(np.min(sx_), np.max(sx_))
    # print("shape", sra_table.shape)
    # print(sra_table.shape)
    # print(inner_.shape)
    # print(outer_.shape)
    return RegularGridInterpolator((kx_, sx_, inner_, outer_), sra_table)  # method = "nearest" fill_value = None , bounds_error=False


def soil_root_interface(rx, sx, inner_kr, rho, sp):
    """
    finds potential at the soil root interface
    
    rx             xylem matric potential [cm]
    sx             bulk soil matric potential [cm]
    inner_kr       root radius times hydraulic conductivity [cm/day] 
    rho            geometry factor [1]
    sp             soil van Genuchten parameters (type vg.Parameters)
    """
    k_soilfun = lambda hsoil, hint: (vg.fast_mfp[sp](hsoil) - vg.fast_mfp[sp](hint)) / (hsoil - hint)
    # rho = outer_r / inner_r  # Eqn [5]
    rho2 = rho * rho  # rho squared
    # b = 2 * (rho2 - 1) / (1 + 2 * rho2 * (np.log(rho) - 0.5))  # Eqn [4]
    b = 2 * (rho2 - 1) / (1 - 0.53 * 0.53 * rho2 + 2 * rho2 * (np.log(rho) + np.log(0.53)))  # Eqn [7]
    fun = lambda x: (inner_kr * rx + b * sx * k_soilfun(sx, x)) / (b * k_soilfun(sx, x) + inner_kr) - x
    rsx = fsolve(fun, rx)
    return rsx


def soil_root_interface_table(rx, sx, inner_kr_, rho_, f):
    """
    finds potential at the soil root interface
        
    rx             xylem matric potential [cm]
    sx             bulk soil matric potential [cm]
    inner_kr       root radius times hydraulic conductivity [cm/day] 
    rho            geometry factor [1]
    sp             soil van Genuchten parameters (type vg.Parameters)
    """
    rsx = f((rx, sx, inner_kr_ , rho_))
    return rsx

def simulate_const(s, r, sra_table_lookup, trans, sim_time, dt):
    """     
    simulates the coupled scenario       
        root architecture is not gowing  
        conductivities are not changing over time
        
    s                            soil model (RichardsWrapper(RichardsSP()))
    r                            xylem flux model (XylemFluxPython wrapping MappedSegments mapped to soil @param s)
    sra_table_lookup             potentials a root soil interface    
    trans                        daily transpiration
    sim_time                     simulation time
    dt                           time step
    
    TODO recyle factorisation of left hand side ... 
    """
    wilting_point = -15000 # cm
    skip = 1  # for output and results, skip iteration
    rs_age = 0.     
    max_iter = 1000 # maximum for fix point iteration
    
    start_time = timeit.default_timer()
    
    nodes = r.rs.nodes
    segs = r.rs.segments
    ns = len(segs)
    mapping = np.array([r.rs.seg2cell[j] for j in range(0, ns)])
    cell_centers = s.getCellCenters()
    cell_centers_z = np.array([cell_centers[mapping[j]][2] for j in range(0, ns)])
    seg_centers_z = np.array([0.5 * (nodes[s.x].z + nodes[s.y].z) for s in segs])
    
    outer_r = r.rs.segOuterRadii()  
    inner_r = r.rs.radii
    types = r.rs.subTypes
    rho_ = np.divide(outer_r, np.array(inner_r))
       
    psi_x_, psi_s_, sink_ , x_, y_, psi_s2_ = [], [], [], [], [], []  # for post processing

    sx = s.getSolutionHead()  # inital condition, solverbase.py
    hsb = np.array([sx[mapping[j]][0] for j in range(0, ns)])  # soil bulk matric potential per segment
    rsx = hsb.copy()  # initial values for fix point iteration
    
    rx = r.solve(rs_age, -trans * sinusoidal2(0, dt), 0., rsx, False, wilting_point, soil_k = [])
    rx_old = rx.copy()
    
    kr_ = np.array([r.kr_f(rs_age, types[j]) for j in range(0, len(outer_r))])  # here const
    inner_kr_ = np.multiply(inner_r, kr_)  # multiply for table look up; here const
    
    N = int(np.ceil(sim_time / dt))  # number of iterations

    for i in range(0, N):
    
        t = i * dt  # current simulation time
    
        wall_iteration = timeit.default_timer()
        wall_fixpoint = timeit.default_timer()
    
        err = 1.e6  # cm
        c = 0
        while err > 1 and c < max_iter:
    
            """ interpolation """
            wall_interpolation = timeit.default_timer()
            rx_ = rx[1:] - seg_centers_z  # from total matric potential to matric potential
            hsb_ = hsb - cell_centers_z  # from total matric potential to matric potential
            rsx = soil_root_interface_table(rx_ , hsb_, inner_kr_, rho_, sra_table_lookup)  
            # rsx = soil_root_interface(rx[1:] , hsb, inner_kr_, rho_, sra_table_lookup) # put soil to sra_table_lookup
            rsx = rsx + seg_centers_z  # from matric potential to total matric potential
            wall_interpolation = timeit.default_timer() - wall_interpolation
    
            """ xylem matric potential """
            wall_xylem = timeit.default_timer()
            rx = r.solve(rs_age, -trans * sinusoidal2(t, dt), 0., rsx, False, wilting_point, soil_k = [])  # xylem_flux.py, cells = False
            err = np.linalg.norm(rx - rx_old)
            wall_xylem = timeit.default_timer() - wall_xylem
    
            rx_old = rx.copy()
            c += 1
    
        print(i, c, "iterations", wall_interpolation / (wall_interpolation + wall_xylem), wall_xylem / (wall_interpolation + wall_xylem))
    
        wall_fixpoint = timeit.default_timer() - wall_fixpoint
    
        wall_soil = timeit.default_timer()
        fluxes = r.segFluxes(rs_age, rx, rsx, approx = False, cells = False)
        collar_flux = r.collar_flux(rs_age, rx.copy(), rsx.copy(), k_soil = [], cells = False)  # validity checks
        err = np.linalg.norm(np.sum(fluxes) - collar_flux)
        if err > 1.e-10:
            print("error: summed root surface fluxes and root collar flux differ" , err)
            raise
        err2 = np.linalg.norm(-trans * sinusoidal2(t, dt) - collar_flux)
        if r.last == "neumann":
            if err2 > 1.e-10:
                print("error: potential transpiration differs root collar flux in Neumann case" , err2)
                raise
        soil_fluxes = r.sumSegFluxes(fluxes)
        s.setSource(soil_fluxes.copy())  # richards.py
        s.solve(dt)
        sx = s.getSolutionHead()[:, 0]  # richards.py
        hsb = np.array([sx[mapping[j]] for j in range(0, ns)])  # soil bulk matric potential per segment
        wall_soil = timeit.default_timer() - wall_soil
    
        wall_iteration = timeit.default_timer() - wall_iteration
    
        """ remember results ... """
        if i % skip == 0:
            psi_x_.append(rx.copy())  # cm
            psi_s_.append(rsx.copy())  # cm
            sink_.append(fluxes.copy())  # cm3/day
            x_.append(t)  # day
            y_.append(np.sum(fluxes))  # cm3/day
            psi_s2_.append(np.array(sx))  # cm
    
    print ("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s")
    
    return psi_x_, psi_s_, sink_, x_, y_, psi_s2_
    
    
    
