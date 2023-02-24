"""
functions for the cylindrical coupling approach
"""
import numpy as np
import timeit
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank(); max_rank = comm.Get_size()

from xylem_flux import sinusoidal2
import vtk_plot as vtk
import plantbox as pb
import evapotranspiration as evap


def simulate_const(s, rs, sim_time, dt, trans_f, rs_age, type):
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
    sx = s.getSolutionHead_()  # initial condition of soil [cm]
    cc = s.getSolution_(1)  # kg/m3

    segs = rs.rs.rs.segments  # this is not nice (rs RhizoMappedSegments, rs.rs XylemFluxPython, rs.rs.rs MappedRootSystem(MappedSegments)
    seg2cell = rs.rs.rs.seg2cell
    cell2seg = rs.rs.rs.cell2seg
    ns = len(segs)

    dcyl = int(np.floor(ns / max_rank))
    if rank + 1 == max_rank:
        rs.initialize(s.soils[0], sx, np.array(range(rank * dcyl, ns)), cc)
        print ("\nInitialized final rank {:g}/{:g} [{:g}-{:g}] in {:g} s\n".format(rank + 1, max_rank, rank * dcyl, ns, timeit.default_timer() - start_time))
    else:
        rs.initialize(s.soils[0], sx, np.array(range(rank * dcyl, (rank + 1) * dcyl)), cc)
        print ("\nInitialized rank {:g}/{:g} [{:g}-{:g}] in {:g} s\n".format(rank + 1, max_rank, rank * dcyl, (rank + 1) * dcyl, timeit.default_timer() - start_time))

    r = rs.rs  # rename (XylemFluxPython)

    rsx = rs.get_inner_heads()  # matric potential at the root soil interface, i.e. inner values of the cylindric models (not extrapolation to the interface!) [cm]
    if type == 1:
        rsc = [cc[seg2cell[i]] for i in range(0, ns)]  # kg/m3
    else:
        rsc = rs.get_inner_solutes(1)  # kg/m3

    if rank == 0:
        print("\nINITIAL root-soil interface matric potentials", np.min(rsx), np.max(rsx), np.min(rsc), np.max(rsc))
        rx = r.solve(rs_age, trans_f(rs_age + 0., dt), 0., rsx, cells = False, wilting_point = wilting_point, soil_k = [])
    else:
        rx = None

    psi_x_, psi_s_, sink_ , x_, y_, psi_s2_, soil_c_, c_ = [], [], [], [], [], [], [], []  # for post processing

    cell_volumes = s.getCellVolumes_()  # cm3
    net_flux = np.zeros(cell_volumes.shape)

    N = int(np.ceil(sim_time / dt))  # number of iterations

    """ simualtion loop """
    for i in range(0, N):

        t = i * dt  # current simulation time

        """ 1. xylem model """
        # if rank == 0:
        #     print("[x", end = "")
        wall_xylem = timeit.default_timer()

        rsx = rs.get_inner_heads(1)  # matric potential at the root soil interface, 2nd node  [cm] (in rhizo_models.py)
        if type == 1:
            rsc = [cc[seg2cell[i]] for i in range(0, ns)]  # kg/m3
            # print(np.min(rsc), np.max(rsc))
        else:
            rsc = rs.get_inner_solutes(1)  # kg/m3
        comm.barrier()

        if rank == 0:
            rx = r.solve(rs_age + t, trans_f(rs_age + t, dt), 0., rsx, cells = False, wilting_point = wilting_point)  # soil_k = soil_k
            # print("*", end = "")
            seg_fluxes = np.array(r.segFluxes(rs_age + t, rx, rsx, False, False, []))
            # print("\n\nShould agree if not in stress \n", trans_f(t, dt), np.sum(seg_fluxes))
            # print("*", end = "")
            seg_sol_fluxes = np.array(r.solute_fluxes(rsc))  # [g/day]
        else:
            rx = None
            seg_fluxes = None
            seg_sol_fluxes = None

        rx = comm.bcast(rx, root = 0)
        seg_fluxes = comm.bcast(seg_fluxes, root = 0)
        seg_sol_fluxes = comm.bcast(seg_sol_fluxes, root = 0)
        # print("*", end = "")

        if rank == 0:
            collar_flux = r.collar_flux(rs_age + t, rx.copy(), rsx.copy(), k_soil = [], cells = False)  # validity checks
            err = np.linalg.norm(np.sum(seg_fluxes) - collar_flux)
            if err > 1.e-6:
                print("error: summed root surface fluxes and root collar flux differ" , err, r.neumann_ind, collar_flux, np.sum(seg_fluxes))
            err2 = np.linalg.norm(trans_f(rs_age + t, dt) - collar_flux)
            if r.last == "neumann":
                if err2 > 1.e-6:
                    print("error: potential transpiration differs root collar flux in Neumann case" , err2)

        comm.barrier()

        # if rank == 0:
        #     print("]", end = "")
        wall_xylem = timeit.default_timer() - wall_xylem

        """ 2. local soil models """
        # if rank == 0:
        #     print("[l", end = "")
        wall_local = timeit.default_timer()

        if rank == 0:
            for key, value in cell2seg.items():  # check cell2seg
                if key < 0:
                    nodes = rs.rs.rs.nodes
                    print("key is negative", key)
                    print("segments", cell2seg[key])
                    print("coresponding nodes")
                    for s in cell2seg[key]:
                        print(segs[s])
                        print(nodes[segs[s].x], nodes[segs[s].y])
                    ana = pb.SegmentAnalyser(rs.rs.rs.mappedSegments())
                    ana.addCellIds(rs.rs.rs.mappedSegments())
                    vp.plot_roots(ana, "cell_id")
                # for v in value:
                #     if not (v >= 0 and v < ns):
                #         print("segments", ns)
                #         print("value", v)
                #         print("age", rs_age, "rank", rank, "dt", dt)
                #         print("Mapping error in cell", key)
                #         print("mapped to segments", value)
                #         ana = pb.SegmentAnalyser(rs.rs.rs.mappedSegments())
                #         ana.addCellIds(rs.rs.rs.mappedSegments())
                #         vp.plot_roots(ana, "cell_id")

            proposed_outer_fluxes = r.splitSoilFluxes(net_flux / dt, split_type)  # if this fails, a segment is not mapped, i.e. out of soil domain
        else:
            proposed_outer_fluxes = None
        proposed_outer_fluxes = comm.bcast(proposed_outer_fluxes, root = 0)

        seg_rx = np.array([0.5 * (rx[seg.x] + rx[seg.y]) for seg in segs])
        rs.solve(dt, seg_rx, proposed_outer_fluxes)  # left dirchlet, right neumann <----
        # TODO mass_net_fluxes

        wall_local = timeit.default_timer() - wall_local
        # if rank == 0:
        #     print("]", end = "")

        """ 3a. macroscopic soil model """
        # if rank == 0:
        #     print("[m", end = "")
        wall_macro = timeit.default_timer()

        water_content = np.array(s.getWaterContent_())  # theta per cell [1]
        soil_water = np.multiply(water_content, cell_volumes)  # water per cell [cm3]
        # TODO for nitrate

        soil_fluxes = r.sumSegFluxes(seg_fluxes)  # [cm3/day]  per soil cell
        soil_sol_fluxes = r.sumSegFluxes(seg_sol_fluxes)  # [g/day]
        evap.add_nitrificatin_source(s, soil_sol_fluxes, nit_flux = 0.)  # 1.e-5
        s.setSource(soil_fluxes.copy(), eq_idx = 0)  # [cm3/day], in moduels/richards.py
        s.setSource(soil_sol_fluxes.copy(), eq_idx = 1)  # [g/day], in moduels/richards.py
        s.solve(dt)  # in modules/solverbase.py

        wall_macro = timeit.default_timer() - wall_macro
        # if rank == 0:
        #     print("]", end = "")

        """ 3b. calculate net fluxes """
        # if rank == 0:
        #     print("[n", end = "")
        wall_netfluxes = timeit.default_timer()
        water_content = np.array(s.getWaterContent_())
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

        sx = s.getSolutionHead_()
        cc = s.getSolution_(1)  # [kg/m3]
        if rank == 0:
            sink = np.zeros(sx.shape)
            for k, v in soil_fluxes.items():
                sink[k] += v
            sink_.append(sink)  # cm3/day (per soil cell)
            x_.append(rs_age + t)  # day
            y_.append(np.sum(sink))  # cm3/day
            c_.append(-np.sum(seg_sol_fluxes))  # [cm3/day]
            psi_s2_.append(sx)  # cm (per soil cell)
            soil_c_.append(cc)  # [kg/m3]

        if i % skip == 0 and rank == 0:
            print("{:g}/{:g} iterations".format(i, N), "time", rs_age + t, "wallt times",
                  wall_xylem / (wall_xylem + wall_local + wall_macro + wall_netfluxes),
                  wall_local / (wall_xylem + wall_local + wall_macro + wall_netfluxes),
                  wall_macro / (wall_xylem + wall_local + wall_macro + wall_netfluxes),
                  wall_netfluxes / (wall_xylem + wall_local + wall_macro + wall_netfluxes),
                  "segments:", rs.getNumberOfMappedSegments(), "root collar:", rx[0], "\n")
            print("time", rs_age + t, "rsx", np.min(rsx), np.max(rsx), "ccx", np.min(cc), np.max(cc), "trans", trans_f(rs_age + t, dt))
            psi_x_.append(rx.copy())  # cm (per root node)
            psi_s_.append(rsx.copy())  # cm (per root segment)

        # if rank == 0:
        #     print("]")

    if rank == 0:
        print ("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s")

    return psi_x_, psi_s_, sink_, x_, y_, psi_s2_, soil_c_, c_

