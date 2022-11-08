"""
functions for the cylindrical coupling approach
"""
import numpy as np
import timeit
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank(); max_rank = comm.Get_size()

from xylem_flux import sinusoidal2
import vtk_plot as vtk
import plantbox as pb


def simulate_const(s, rs, sim_time, dt, trans_f, rs_age):
    """     
    simulates the coupled scenario       
        root architecture is not gowing  
        conductivities are not changing over time
        
    s
    rs
    trans
    sim_time    
    dt
    trans_f
    rs_age
    """
    wilting_point = -15000  # cm
    skip = 10  # 3 * 6  # for output and results, skip iteration
    split_type = 1  # type 0 == volume, type 1 == surface, type 2 == length

    """ 
    Initialize local soil models (around each root segment) 
    """
    start_time = timeit.default_timer()
    sx = s.getSolutionHead()  # initial condition of soil [cm]
    sx = comm.bcast(sx, root = 0)  # soil part might run parallel
    cc = s.getSolution(1)  # UNITS?
    cc = comm.bcast(cc, root = 0)  # soil part might run parallel
    comm.barrier()
    ns = len(rs.segments)
    dcyl = int(np.floor(ns / max_rank))
    if rank + 1 == max_rank:
        rs.initialize(s.soils[0], sx, np.array(range(rank * dcyl, ns)), cc)
        print ("\nInitialized final rank {:g}/{:g} [{:g}-{:g}] in {:g} s\n".format(rank + 1, max_rank, rank * dcyl, ns, timeit.default_timer() - start_time))
    else:
        rs.initialize(s.soils[0], sx, np.array(range(rank * dcyl, (rank + 1) * dcyl)), cc)
        print ("\nInitialized rank {:g}/{:g} [{:g}-{:g}] in {:g} s\n".format(rank + 1, max_rank, rank * dcyl, (rank + 1) * dcyl, timeit.default_timer() - start_time))

    r = rs.rs  # rename (XylemFluxPython)
    comm.barrier()

    rsx = rs.get_inner_heads()  # matric potential at the root soil interface, i.e. inner values of the cylindric models (not extrapolation to the interface!) [cm]
    rsx = comm.bcast(rsx, root = 0)  # Soil part runs parallel
    print("\nINITIAL root-soil interface matric potentials", np.min(rsx), np.max(rsx))
    if rank == 0:
        rx = r.solve(rs_age, trans_f(rs_age + 0., dt), 0., rsx, cells = False, wilting_point = wilting_point, soil_k = [])
    else:
        rx = None

    psi_x_, psi_s_, sink_ , x_, y_, psi_s2_, soil_c_ = [], [], [], [], [], [], []  # for post processing

    cell_volumes = s.getCellVolumes()  # cm3
    cell_volumes = comm.bcast(cell_volumes, root = 0)
    net_flux = np.zeros(cell_volumes.shape)

    segs = r.rs.segments

    N = int(np.ceil(sim_time / dt))  # number of iterations

    """ simualtion loop """
    for i in range(0, N):

        t = i * dt  # current simulation time

        """ 1. xylem model """
        # if rank == 0:
        #     print("[x", end = "")
        wall_xylem = timeit.default_timer()

        rsx = rs.get_inner_heads(1)  # matric potential at the root soil interface, 2nd node  [cm]
        rsc = rs.get_inner_solutes(1)  # TODO UNITS !
        comm.barrier()
        if rank == 0:
            rx = r.solve(rs_age + t, trans_f(rs_age + t, dt), 0., rsx, cells = False, wilting_point = wilting_point)  # soil_k = soil_k
            # print("*", end = "")
            seg_fluxes = np.array(r.segFluxes(rs_age + t, rx, rsx, False, False, []))
            # print("\n\nShould agree if not in stress \n", trans_f(t, dt), np.sum(seg_fluxes))
            # print("*", end = "")

            seg_sol_fluxes = np.array(r.solute_fluxes(rsc))
        else:
            rx = None
            seg_fluxes = None
        comm.barrier()
        rx = comm.bcast(rx, root = 0)
        # print("*", end = "")
        seg_fluxes = comm.bcast(seg_fluxes, root = 0)
        seg_sol_fluxes = comm.bcast(seg_sol_fluxes, root = 0)

        # if rank == 0:
        #     print("]", end = "")
        wall_xylem = timeit.default_timer() - wall_xylem

        """ 2. local soil models """
        # if rank == 0:
        #     print("[l", end = "")
        wall_local = timeit.default_timer()

        seg_rx = np.array([0.5 * (rx[s.x] + rx[s.y]) for s in segs])
        proposed_outer_fluxes = r.splitSoilFluxes(net_flux / dt, split_type)  # if this fails, a segment is not mapped, i.e. out of soil domain
        rs.solve(dt, seg_rx, proposed_outer_fluxes)  # left dirchlet, right neumann <----
        # TODO mass_net_fluxes

        wall_local = timeit.default_timer() - wall_local
        # if rank == 0:
        #     print("]", end = "")

        """ 3a. macroscopic soil model """
        # if rank == 0:
        #     print("[m", end = "")
        wall_macro = timeit.default_timer()

        water_content = np.array(s.getWaterContent())  # theta per cell [1]
        water_content = comm.bcast(water_content, root = 0)
        soil_water = np.multiply(water_content, cell_volumes)  # water per cell [cm3]
        # TODO mass source

        soil_fluxes = r.sumSegFluxes(seg_fluxes)  # [cm3/day]  per soil cell
        soil_sol_fluxes = r.sumSegFluxes(seg_sol_fluxes)
        s.setSource(soil_fluxes.copy(), eq_idx = 0)  # [cm3/day], in moduels/richards.py
        # s.setSource(soil_sol_fluxes.copy(), eq_idx = 1)
        s.solve(dt)  # in modules/solverbase.py

        wall_macro = timeit.default_timer() - wall_macro
        # if rank == 0:
        #     print("]", end = "")

        """ 3b. calculate net fluxes """
        # if rank == 0:
        #     print("[n", end = "")
        wall_netfluxes = timeit.default_timer()
        water_content = np.array(s.getWaterContent())
        water_content = comm.bcast(water_content, root = 0)
        new_soil_water = np.multiply(water_content, cell_volumes)  # calculate net flux
        net_flux = new_soil_water - soil_water  # change in water per cell [cm3]
        for k, root_flux in soil_fluxes.items():
            net_flux[k] -= root_flux * dt
            # net_flux[k] = 0.  # same as SRA approach (currently using a no-flux at boundary)
        # TODO mass net flux
        soil_water = new_soil_water

        wall_netfluxes = timeit.default_timer() - wall_netfluxes
        # if rank == 0:
        #     print("]", end = "")
        # if rank == 0:
        #     if rs_age + t > 1.5:
        #         print(ns)
        #         print(len(rs.cyls))
        #         min_b = [-19, -2.5, -200.]  # for soybean
        #         max_b = [19, 2.5, 0.]
        #         cell_number = [1, 1, 200]
        #         # vtk.plot_roots_and_soil(rs, "fluxes", seg_fluxes.copy(), s, True, min_b, max_b, cell_number, "nice_plot")
        #         # vtk.plot_roots(pd, p_name:str, win_title:str = "", render:bool = True):
        #         ind0 = s.pick([0, 0, -3.5])
        #         # ind1 = s.pick([0, 0, -15.])
        #         # ind2 = s.pick([0, 0, -25.])
        #         print("cell0", ind0)
        #         # print("cell1", ind1)
        #         # print("cell2", ind2)
        #         cell2seg = r.rs.cell2seg
        #         segs0 = cell2seg[ind0]
        #         # # segs1 = cell2seg[ind1]
        #         # # segs2 = cell2seg[ind2]
        #         print(segs0)
        #         for i in segs0:
        #             rs.plot_cylinder(i)
        #         dd

        """ remember results ... """
        # if rank == 0:
        #     print("[r", end = "")

        sx = s.getSolutionHead()
        if rank == 0:
            sx = sx[:, 0]
            sink = np.zeros(sx.shape)
            for k, v in soil_fluxes.items():
                sink[k] += v
            sink_.append(sink)  # cm3/day (per soil cell)
            x_.append(rs_age + t)  # day
            y_.append(np.sum(sink))  # cm3/day
            psi_s2_.append(sx)  # cm (per soil cell)
            cc = s.getSolution(1)[:, 0]  # UNITS?
            soil_c_.append(cc)

        if i % skip == 0 and rank == 0:
            print("{:g}/{:g} iterations".format(i, N), "time", rs_age + t, "wallt times",
                  wall_xylem / (wall_xylem + wall_local + wall_macro + wall_netfluxes),
                  wall_local / (wall_xylem + wall_local + wall_macro + wall_netfluxes),
                  wall_macro / (wall_xylem + wall_local + wall_macro + wall_netfluxes),
                  wall_netfluxes / (wall_xylem + wall_local + wall_macro + wall_netfluxes),
                  "segments:", rs.getNumberOfSegments(), "root collar:", rx[0], "\n")
            print("rsx", np.min(rsx), np.max(rsx), "trans", trans_f(rs_age + t, dt), "time", rs_age + t)
            psi_x_.append(rx.copy())  # cm (per root node)
            psi_s_.append(rsx.copy())  # cm (per root segment)

        # if rank == 0:
        #     print("]")

    if rank == 0:
        print ("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s")

    return psi_x_, psi_s_, sink_, x_, y_, psi_s2_, soil_c_

