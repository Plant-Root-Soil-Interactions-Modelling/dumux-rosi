"""
functions for the cylindrical coupling approach
"""
import numpy as np
import timeit
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank(); max_rank = comm.Get_size()

from xylem_flux import sinusoidal2


def simulate_const(s, rs, trans, sim_time, dt, trans_f = None):
    """     
    simulates the coupled scenario       
        root architecture is not gowing  
        conductivities are not changing over time
        
    s
    rs
    trans
    sim_time    
    dt
    """
    wilting_point = -15000  # cm
    skip = 1  # 3 * 6  # for output and results, skip iteration
    rs_age = 0.  # day
    split_type = 0  # type 0 == volume, type 1 == surface, type 2 == length

    """ set defalut potential transpiration """
    if not trans_f:
        trans_f = lambda age, dt:-trans * sinusoidal2(age, dt)

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

    wv_old = rs.get_water_volume()

    rsx = rs.get_inner_heads()  # matric potential at the root soil interface, i.e. inner values of the cylindric models (not extrapolation to the interface!) [cm]
    print("\nINITIAL root-soil interface matric potentials", np.min(rsx), np.max(rsx))
    if rank == 0:
        rx = r.solve(rs_age, trans_f(0., dt), 0., rsx, cells = False, wilting_point = wilting_point, soil_k = [])
    else:
        rx = None

    rs.solve(dt, rsx, np.zeros(rsx.shape))  # left dirichlet, right neumann
    fluxes = r. segFluxes(rs_age, rx, rsx, False, False, [])
    print("INITIAL segment fluxes", np.min(fluxes), np.max(fluxes))
    print("INITIAL root-soil interface matric potentials", np.min(rs.get_inner_heads()), np.max(rs.get_inner_heads()))
    print("INITIAL change", np.min(rs.get_water_volume() - wv_old), np.max(rs.get_water_volume() - wv_old))
    # rs.plot_cylinders()

    # outer_r = r.rs.segOuterRadii()
    inner_surface = 2 * np.pi * np.multiply(rs.radii, rs.seg_length)  # cm2 (for approximation)
    nodes = r.rs.nodes
    segs = r.rs.segments

    types = r.rs.subTypes

    N = int(np.ceil(sim_time / dt))  # number of iterations

    """ simualtion loop """
    for i in range(0, N):

        t = i * dt  # current simulation time

        """ 1. xylem model """
        if rank == 0:
            print("[x", end = "")
        wall_xylem = timeit.default_timer()

        wv = rs.get_water_volume()
        # print(np.min(wv), np.max(wv))
        seg_fluxes = ((wv - wv_old) / dt) / (24 * 3600)  # [cm3 / day]
        # seg_fluxes = np.array(r.segFluxes(rs_age + t, rx, rsx, False, False, []))

        print("based on water volume", np.sum(seg_fluxes), "segment fluxes", seg_fluxes)
        seg_fluxes_ = np.array(r.segFluxes(rs_age + t, rx, rsx, False, False, []))
        print("based on segment flux", np.sum(seg_fluxes_), "segment fluxes", seg_fluxes_)

        print ("\nvolume change in soil", np.min(wv - wv_old), np.max(wv - wv_old))
        print("roots-soil fluxes", np.min(seg_fluxes), np.max(seg_fluxes))

        kr_ = np.array([r.kr_f(rs_age + t, types[j]) for j in range(0, len(inner_surface))])
        kx_ = np.array([r.kx_f(rs_age + t, types[j]) for j in range(0, len(inner_surface))])

        # # exact
        # f = -2. * np.pi * np.multiply(rs.radii, kr_)  # double f = -perimeter*kr; // flux is proportional to f // *rho*g
        # tau = np.sqrt(2. * np.pi * np.multiply(rs.radii, np.divide(kr_, kx_)))  # double tau = std::sqrt(perimeter*kr/kx); // sqrt(c) [cm-1]
        # taul = np.multiply(tau, rs.seg_length)
        # d = np.exp(-taul) - np.exp(taul)  # double d = std::exp(-tau*l)-std::exp(tau*l); /
        # aa = -np.divide(seg_fluxes, 2.*np.ones(seg_fluxes.shape) - np.exp(taul) - np.exp(-taul))
        # bb = np.divide(np.multiply(tau, d), f)
        #
        # rsx = np.array([-0.5 * (np.multiply(aa, bb) - rx[s.x] - rx[s.y]) for s in segs])
        # rsx0 = np.array([0.5 * (rx[s.x] + rx[s.y]) for s in segs])  # for aa == 0

        # # approximation
        flux = np.divide(seg_fluxes, inner_surface)  # [cm3 / cm2 / day]
        rsx = np.divide(flux, kr_) + np.array([0.5 * (rx[s.x] + rx[s.y]) for s in segs])

        # # should be const for const kr and kx ...
        # print("taul", np.min(taul), np.max(taul))
        # print("f", np.min(f), np.max(f))
        # print("d", np.min(d), np.max(d))
        # print("radial conductivity", np.min(kr_), np.max(kr_))
        # print("axial conductivity", np.min(kx_), np.max(kx_))

        print("potential transpiration", trans_f(t - dt, dt), np.sum(seg_fluxes))
        print("xylem matric potentials", np.min(rx), np.max(rx))
        print("root-soil interface matric potentials", np.min(rsx), np.max(rsx))
        # print("root-soil interface matric potentials", np.min(rsx0), np.max(rsx0))
        print()

        wv_old = wv

        rx = comm.bcast(rx, root = 0)
        # soil_k = rs.get_soil_k(rx)
        if rank == 0:
            # soil_k = np.divide(vg.hydraulic_conductivity(rsx, soil), rs.radii)  # only valid for homogenous soil
            rx = r.solve(rs_age + t, trans_f(t, dt), 0., rsx, cells = False, wilting_point = wilting_point)  # soil_k = soil_k

            seg_fluxes_ = np.array(r.segFluxes(rs_age + t, rx, rsx, False, False, []))
            print("\n\nThey should agree.... \n", trans_f(t, dt), np.sum(seg_fluxes_))
            # proposed_inner_fluxes = r.segFluxes(rs_age + t, rx.copy(), rsx.copy(), approx = False, cells = False, soil_k = soil_k.copy())  # [cm3/day]
        else:
            # proposed_inner_fluxes = None
            rx = None
        # # validity check
        # if rank == 0:
        #     collar_flux = r.collar_flux(rs_age + t, rx.copy(), rsx.copy(), k_soil = soil_k.copy(), cells = False)
        #     err = np.linalg.norm(np.sum(proposed_inner_fluxes) - collar_flux)
        #     if err > 1.e-8:
        #         print("error: summed root surface fluxes and root collar flux differ" , err)

        if rank == 0:
            print("]", end = "")
        wall_xylem = timeit.default_timer() - wall_xylem

        """ 2. local soil models """
        if rank == 0:
            print("[l", end = "")
        wall_local = timeit.default_timer()

        # proposed_inner_fluxes = comm.bcast(proposed_inner_fluxes, root = 0)
        proposed_outer_fluxes = r.splitSoilFluxes(net_flux / dt, split_type)
        rs.solve(dt, rsx, proposed_outer_fluxes)  # left dirchlet, right neumann <----

        # rs.plot_cylinders()

        wall_local = timeit.default_timer() - wall_local
        if rank == 0:
            print("]", end = "")

        """ 3a. macroscopic soil model """
        if rank == 0:
            print("[m", end = "")
        wall_macro = timeit.default_timer()

        water_content = np.array(s.getWaterContent())  # theta per cell [1]
        water_content = comm.bcast(water_content, root = 0)
        soil_water = np.multiply(water_content, cell_volumes)  # water per cell [cm3]
        soil_fluxes = r.sumSegFluxes(seg_fluxes)  # [cm3/day]  per soil cell
        s.setSource(soil_fluxes.copy())  # [cm3/day], in moduels/richards.py
        s.solve(dt)  # in modules/solverbase.py

        wall_macro = timeit.default_timer() - wall_macro
        if rank == 0:
            print("]", end = "")

        """ 3b. calculate net fluxes """
        if rank == 0:
            print("[n", end = "")
        wall_netfluxes = timeit.default_timer()

        water_content = np.array(s.getWaterContent())
        water_content = comm.bcast(water_content, root = 0)
        new_soil_water = np.multiply(water_content, cell_volumes)  # calculate net flux
        net_flux = new_soil_water - soil_water  # change in water per cell [cm3]
        for k, root_flux in soil_fluxes.items():
            net_flux[k] -= root_flux * dt
            # net_flux[k] = 0.  # same as SRA approach (currently using a no-flux at boundary)
        soil_water = new_soil_water

        wall_netfluxes = timeit.default_timer() - wall_netfluxes
        if rank == 0:
            print("]", end = "")

        """ remember results ... """
        if rank == 0:
            print("[r", end = "")

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
                print("{:g}/{:g} iterations".format(i, N), "time", t, "wallt times",
                      wall_xylem / (wall_xylem + wall_local + wall_macro + wall_netfluxes),
                      wall_local / (wall_xylem + wall_local + wall_macro + wall_netfluxes),
                      wall_macro / (wall_xylem + wall_local + wall_macro + wall_netfluxes),
                      wall_netfluxes / (wall_xylem + wall_local + wall_macro + wall_netfluxes),
                      "number of segments", rs.getNumberOfSegments())

        if rank == 0:
            print("]")

    if rank == 0:
        print ("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s")

    return psi_x_, psi_s_, sink_, x_, y_, psi_s2_

