"""
functions for the steady rate approach
"""
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import fsolve
import timeit

from xylem_flux import sinusoidal2

import van_genuchten as vg

from sra import soil_root_interface
from sra import soil_root_interface_table


def layer_sumB(r, b):
    """ sums up the xylem matric potential vector"""
    segs = r.rs.segments
    ns = len(segs)
    mapping = np.array([r.rs.seg2cell[j] for j in range(0, ns)])
    cell_max = np.max(mapping)
    b_ups = np.zeros((cell_max + 1,))
    mapping[0] = -1  # ignrore first entry
    b_ups[0] = b[0]
    for i in range(0, cell_max):
        b_ups[i + 1] = np.sum(b[mapping == i])
    return b_ups


def layer_sumA(r, A):
    """ sums up the system matrix """
    segs = r.rs.segments
    ns = len(segs)
    mapping = np.array([r.rs.seg2cell[j] for j in range(0, ns)])
    rx_ = rx[1:]
    cell_max = np.max(mapping)
    for i in range(0, cell_max):
        rx_small[i] = np.sum(rx_[mapping == i])
    return rx_small


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
    wilting_point = -15000  # cm
    skip = 6  # for output and results, skip iteration
    rs_age = 0.  # day
    max_iter = 1000  # maximum for fix point iteration

    if isinstance(sra_table_lookup, RegularGridInterpolator):
        root_interface = soil_root_interface_table
    else:
        root_interface = soil_root_interface

    start_time = timeit.default_timer()

    nodes = r.rs.nodes
    segs = r.rs.segments
    ns = len(segs)
    mapping = np.array([r.rs.seg2cell[j] for j in range(0, ns)])  # because seg2cell is a map
    cell_centers = s.getCellCenters()
    cell_centers_z = np.array([cell_centers[j][2] for j in mapping])
    seg_centers_z = np.array([0.5 * (nodes[s.x].z + nodes[s.y].z) for s in segs])

    outer_r = r.rs.segOuterRadii()
    inner_r = r.rs.radii
    types = r.rs.subTypes
    rho_ = np.divide(outer_r, np.array(inner_r))
    rho_ = np.minimum(rho_, np.ones(rho_.shape) * 200)  ############################################ (too keep within table)

    psi_x_, psi_s_, sink_ , x_, y_, psi_s2_ = [], [], [], [], [], []  # for post processing

    sx = s.getSolutionHead()  # inital condition, solverbase.py
    hsb = np.array([sx[j][0] for j in mapping])  # soil bulk matric potential per segment
    rsx = hsb.copy()  # initial values for fix point iteration

    # r.init_solve_static(rs_age, rsx, False, wilting_point, soil_k = [])  # speed up & and forever static...

    rx = r.solve(rs_age, -trans * sinusoidal2(0, dt), 0., rsx, False, wilting_point, soil_k = [])
    print(rx.shape)
    rx = layer_sumB(r, rx)
    print(rx.shape)
    rx_old = rx.copy()
    ddd

    kr_ = np.array([r.kr_f(rs_age, types[j]) for j in range(0, len(outer_r))])  # here const
    inner_kr_ = np.multiply(inner_r, kr_)  # multiply for table look up; here const

    N = int(np.ceil(sim_time / dt))  # number of iterations

    """ simulation loop """
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
            rsx = root_interface(rx_ , hsb_, inner_kr_, rho_, sra_table_lookup)
            rsx = rsx + seg_centers_z  # from matric potential to total matric potential
            wall_interpolation = timeit.default_timer() - wall_interpolation

            """ xylem matric potential """
            wall_xylem = timeit.default_timer()
            rx = r.solve(rs_age, -trans * sinusoidal2(t, dt), 0., rsx, False, wilting_point, soil_k = [])  # xylem_flux.py, cells = False
            err = np.linalg.norm(rx - rx_old)
            wall_xylem = timeit.default_timer() - wall_xylem

            rx_old = rx.copy()
            c += 1

        wall_fixpoint = timeit.default_timer() - wall_fixpoint

        wall_soil = timeit.default_timer()
        fluxes = r.segFluxes(rs_age, rx, rsx, approx = False, cells = False)
        collar_flux = r.collar_flux(rs_age, rx.copy(), rsx.copy(), k_soil = [], cells = False)  # validity checks
        err = np.linalg.norm(np.sum(fluxes) - collar_flux)
        if err > 1.e-6:
            print("error: summed root surface fluxes and root collar flux differ" , err)
            raise
        err2 = np.linalg.norm(-trans * sinusoidal2(t, dt) - collar_flux)
        if r.last == "neumann":
            if err2 > 1.e-6:
                print("error: potential transpiration differs root collar flux in Neumann case" , err2)
                raise
        soil_fluxes = r.sumSegFluxes(fluxes)
        s.setSource(soil_fluxes.copy())  # richards.py
        s.solve(dt)
        sx = s.getSolutionHead()[:, 0]  # richards.py
        hsb = np.array([sx[j] for j in mapping])  # soil bulk matric potential per segment
        wall_soil = timeit.default_timer() - wall_soil

        wall_iteration = timeit.default_timer() - wall_iteration

        """ remember results ... """
        if i % skip == 0:
            psi_x_.append(rx.copy())  # cm (per root node)
            psi_s_.append(rsx.copy())  # cm (per root segment)
            sink = np.zeros(sx.shape)
            for k, v in soil_fluxes.items():
                sink[k] += v
            sink_.append(sink)  # cm3/day (per soil cell)
            x_.append(t)  # day
            y_.append(np.sum(sink))  # cm3/day
            psi_s2_.append(sx.copy())  # cm (per soil cell)
            print("{:g}/{:g} {:g} iterations".format(i, N, c), wall_interpolation / (wall_interpolation + wall_xylem), wall_xylem / (wall_interpolation + wall_xylem))

    print ("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s")

    return psi_x_, psi_s_, sink_, x_, y_, psi_s2_
