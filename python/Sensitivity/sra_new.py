"""
    Water uptake simulation taking nonlinear rhizosphere resistance into account
"""
import numpy as np
import timeit

import plantbox as pb
import visualisation.vtk_plot as vtk
import functional.van_genuchten as vg
from functional.Perirhizal import PerirhizalPython
import evapotranspiration as evap


def simulate_dynamic(s, r, lookuptable_name, sim_time, dt, trans_f, initial_age = 1.):
    """     
    simulates the coupled scenario       
        
    s                            soil model (RichardsWrapper(RichardsSP()))
    r                            PlantHydraulicModel 
    lookuptable_name             matric potentials at the root soil interface (precomputed in a 4D table)    
    sim_time                     simulation time
    dt                           time step
    trans_f                      potential transpiration function 
    initial_age                       initial root system age  
    
    return: 
    psi_x                        root xylem potential (per node) 
    psi_rsx_                     matric potential at the soil root interface (per segment)
    sink_                        water uptake (per cell)   
    times_                       times
    act_trans_                   actual transpiration 
    psi_s_                       soil matric potential (per cell) 
    vol_                         root system volume 
    surf_                        root system surface [cm2]
    krs_                         root system hydraulic conductivity [cm2/day]
    depth_                       root system depth [cm] 
    """

    wilting_point = -15000  # cm
    skip = 10  # for output and results, skip iteration (TODO)
    matimes_iter = 10  # maximum for fix point iteration

    print("simulate_dynamic starting" , flush = True)
    peri = PerirhizalPython(r.ms)
    peri.open_lookup("data/" + lookuptable_name)
    print("opening look up table:", "data/" + lookuptable_name, flush = True)

    start_time = timeit.default_timer()

    psi_x, psi_rsx_, sink_ , times_, act_trans_, psi_s_ = [], [], [], [], [], []  # for post processing
    vol_ = [[], [], [], [], [], []]
    surf_ = [[], [], [], [], [], []]
    net_change_, krs_, depth_, pot_trans_ = [], [], [], []
    collar_pot_ = []

    rs = r.ms
    nodes = rs.nodes
    segs = rs.segments

    sx = s.getSolutionHead_()  # richards.py
    hsb_ = np.array(rs.getHs(sx))  # matric potential per segment
    rsx = hsb_.copy()  # initial values for fix point iteration
    rx = r.solve(initial_age, trans_f(0, dt), rsx, cells = False)
    rx_ = rx.copy()

    N = int(np.ceil(sim_time / dt))  # number of iterations

    print("Starting simulation loop", flush = True)

    """ simulation loop """
    for i in range(0, N):

        t = i * dt  # current simulation time

        """ grow root system and update everything"""
        rs.simulate(dt, False)

        cell2seg = rs.cell2seg  # for debugging
        mapping = rs.getSegmentMapper()  # because seg2cell is a dict
        sx = s.getSolutionHead_()  # richards.py
        hsb_ = np.array(rs.getHs(sx))

        hsb_ = np.maximum(hsb_, np.ones(hsb_.shape) * -15000.)  ############################################ (too keep within table)
        hsb_ = np.minimum(hsb_, np.zeros(hsb_.shape))  ############################################ (too keep within table)

        rsx = np.hstack((rsx, hsb_[rsx.shape[0]:]))  # initial values for fix point iteration

        outer_r = peri.get_outer_radii("length")
        n = len(rs.radii)
        inner_r = np.array([rs.getEffectvieRadius(i) for i in range(0, n)])  # rs.radii
        types = rs.subTypes
        rho_ = np.divide(outer_r, np.array(inner_r))

        rho_ = np.minimum(rho_, np.ones(rho_.shape) * 200)  ############################################ (too keep within table)

        kr_ = r.params.getKr(initial_age + t)
        inner_kr_ = np.multiply(inner_r, kr_)  # multiply for table look up; here const

        inner_kr_ = np.maximum(inner_kr_, np.ones(inner_kr_.shape) * 1.e-7)  ############################################ (too keep within table)
        inner_kr_ = np.minimum(inner_kr_, np.ones(inner_kr_.shape) * 1.e-4)  ############################################ (too keep within table)

        wall_iteration = timeit.default_timer()

        wall_fixpoint = timeit.default_timer()

        err_ = 1.e6  # cm
        c = 0

        rx = r.solve(initial_age + t, trans_f(t, dt), rsx, cells = False)
        rx_ = rx.copy()

        while err_ > 1 and c < matimes_iter:

            """ interpolation """
            wall_interpolation = timeit.default_timer()
            # print("rx", len(rx[1:]))
            # print("hsb_", len(hsb_))
            # print("inner_kr", len(inner_kr_))
            # print("rho", len(rho_))

            rx = np.maximum(rx, np.ones(rx.shape) * -15000.)
            rx = np.minimum(rx, np.zeros(rx.shape))
            rsx = peri.soil_root_interface_potentials(rx[1:], hsb_, inner_kr_, rho_)
            wall_interpolation = timeit.default_timer() - wall_interpolation

            """ xylem matric potential """
            wall_xylem = timeit.default_timer()
            rx = r.solve_again(initial_age + t, trans_f(t, dt), rsx, cells = False)
            err_ = np.linalg.norm(rx - rx_) / np.sqrt(rx.shape[0])  # rmse
            wall_xylem = timeit.default_timer() - wall_xylem

            rx_ = rx.copy()
            c += 1

        wall_fixpoint = timeit.default_timer() - wall_fixpoint

        wall_soil = timeit.default_timer()
        fluxes = r.radial_fluxes(initial_age + t, rx, rsx)
        collar_flux = r.get_transpiration(initial_age + t, rx.copy(), rsx.copy())  # validity checks
        err = np.linalg.norm(np.sum(fluxes) - collar_flux)
        if err > 1.e-6:
            print("error: summed root surface fluxes and root collar flux differ" , err, r.neumann_ind, collar_flux, np.sum(fluxes))
        if rx[0] > -14999:
            err2 = np.linalg.norm(trans_f(t, dt) - collar_flux)
            if err2 > 1.e-6:
                print("")
                print("simulate_dynamic() potential transpiration differs from root collar flux in Neumann case" , err2, rx[0])
                print("fluxes: ", trans_f(t, dt), collar_flux)
                print("")

        soil_fluxes = r.sumSegFluxes(fluxes.copy())
        s.setSource(soil_fluxes.copy())  # richards.py

        water_ = s.getWaterVolume()
        s.solve(dt, doMPIsolve_ = False)
        net_change_.append(s.getWaterVolume() - water_)

        # for key, value in cell2seg.items():  # check cell2seg
        #     if key < 0:
        #         nodes = r.rs.nodes
        #         print("key is negative", key)
        #         print("segments", cell2seg[key])
        #         print("coresponding nodes")
        #         segs = r.rs.segments
        #         for s in cell2seg[key]:
        #             print(segs[s])
        #             print(nodes[segs[s].x], nodes[segs[s].y])
        #         ana = pb.SegmentAnalyser(r.rs.mappedSegments())
        #         ana.addCellIds(r.rs.mappedSegments())
        #         vtk.plot_roots(ana, "cell_id")

        wall_soil = timeit.default_timer() - wall_soil

        wall_iteration = timeit.default_timer() - wall_iteration

        # if initial_age + t > 24.5:
        #     pass
        #     min_b = [-19, -2.5, -200.]  # for soybean
        #     matimes_b = [19, 2.5, 0.]
        #     cell_number = [1, 1, 200]
        #     vtk.plot_roots_and_soil(rs, "fluxes", fluxes.copy()[1:], s, True, min_b, matimes_b, cell_number, "nice_plot")
        #     # vtk.plot_roots(pd, p_name:str, win_title:str = "", render:bool = True):
        #     # ind0 = s.pick([0, 0, -3.5])
        #     # ind1 = s.pick([0, 0, -15.])
        #     # ind2 = s.pick([0, 0, -25.])
        #     # print("cell0", ind0)
        #     # print("cell1", ind1)
        #     # print("cell2", ind2)
        #     # cell2seg = r.rs.cell2seg
        #     # segs0 = cell2seg[ind0]
        #     # # segs1 = cell2seg[ind1]
        #     # # segs2 = cell2seg[ind2]
        #     # for i in segs:
        #     #     rs.plot_cylinder(i)
        #     dd

        """ remember results ... """
        sink = np.zeros(sx.shape)
        for k, v in soil_fluxes.items():
            sink[k] += v
        times_.append(initial_age + t)  # day
        # act_trans_.append(np.sum(sink))  # cm3/day
        act_trans_.append(np.sum(fluxes))  # cm3/day
        collar_pot_.append(rx[0])

        pot_trans_.append(trans_f(t, dt))  # cm3/day

        if i % skip == 0:

            print("age {:g}".format(initial_age + t), "{:g}/{:g} {:g} iterations, rmse {:g}".format(i, N, c, err_), "; wall times {:g} {:g}".format(
                wall_interpolation / (wall_interpolation + wall_xylem), wall_xylem / (wall_interpolation + wall_xylem)),
                  "number of segments", rs.getNumberOfSegments(), "collar potential {:g}".format(rx[0]), flush = True)

            sink_.append(sink)  # cm3/day (per soil cell)
            psi_s_.append(sx.copy())  # cm (per soil cell)

            ana = pb.SegmentAnalyser(r.ms.mappedSegments())  # VOLUME and SURFACE
            for j in range(0, 6):  # root types
                anac = pb.SegmentAnalyser(ana)
                anac.filter("subType", j)
                vol_[j].append(anac.getSummed("volume"))
                surf_[j].append(anac.getSummed("surface"))
            krs, _ = r.get_krs(initial_age + t)
            krs_.append(krs)  # KRS
            depth_.append(ana.getMinBounds().z)

            """ direct vtp output """
            # psi_x.append(rx.copy())  # cm (per root node)
            # psi_rsx_.append(rsx.copy())  # cm (per root segment)

            # rs = r.ms
            # n = len(rs.radii)
            # radii = np.array([rs.getEffectvieRadius(i) for i in range(0, n)])
            # old_radii = rs.radii
            # rs.radii = radii
            # ana = pb.SegmentAnalyser(rs.mappedSegments())  # VOLUME and SURFACE
            # rs.radii = old_radii
            # ana.addData("rx", rx[1:])
            # ana.addData("rsx", rsx)
            # ana.addAge(initial_age + t)  # "age"
            # ana.addHydraulicConductivities(r.params, initial_age + t)  # "kr", "kx"
            # ana.addFluxes(r, rx[1:], rsx, initial_age + t)  # "axial_flux", "radial_flux"
            # ana.write("results/rs{0:05d}.vtp".format(int(i / skip)), ["radius", "subType", "creationTime", "organType", "rx", "rsx", "age", "kr", "kx", "axial_flux", "radial_flux"])

    # rs = r.ms
    # n = len(rs.radii)
    # radii = np.array([rs.getEffectvieRadius(i) for i in range(0, n)])
    # rs.radii = radii
    # ana = pb.SegmentAnalyser(rs.mappedSegments())  # VOLUME and SURFACE
    # for j in range(0, 6):  # root types
    #     anac = pb.SegmentAnalyser(ana)
    #     anac.filter("subType", j)
    #     vol_[j].append(anac.getSummed("volume"))
    #     surf_[j].append(anac.getSummed("surface"))
    # krs, _ = r.get_krs(initial_age + t)
    # krs_.append(krs)  # KRS
    # depth_.append(ana.getMinBounds().z)
    # ana.addData("rx", rx[1:])
    # ana.addData("rsx", rsx)
    # ana.addHydraulicConductivities(r.params, sim_time)
    # ana.addFluxes(r, rx[1:], rsx, initial_age + t)  # "axial_flux", "radial_flux"
    # vtk.plot_roots(ana, "radial_flux")

    print ("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s")

    # psi_x = []  # what do I do with nested arrays?
    # psi_rsx = []

    return pot_trans_, psi_x, psi_rsx_, sink_, times_, act_trans_, psi_s_, vol_, surf_, krs_, depth_, collar_pot_
    """ 
    pot_trans_                   potential transpiration
    psi_x                        root xylem potential (per node) 
    psi_rsx_                     soil matric potential (per segment)
    sink_                        water uptake (per cell)  
    times_                       times
    act_trans_                   actual transpiration
    psi_s_                       soil matric potentials (per cell) 
    vol_                         root system volume 
    surf_                        root system surface [cm2]
    krs_                         root system hydraulic conductivity [cm2/day]
    depth_                       root system depth [cm]
    collar_pot_                  root collar potential [cm]  
    """
