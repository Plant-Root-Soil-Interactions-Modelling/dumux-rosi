"""
functions for the steady rate approach
"""
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import fsolve
import timeit

import plantbox as pb
import visualisation.vtk_plot as vtk
import functional.van_genuchten as vg
from functional.Perirhizal import PerirhizalPython
import evapotranspiration as evap


def simulate_dynamic(s, r, lookuptable_name, sim_time, dt, trans_f, initial_age = 1., type_ = 1):
    """     
    simulates the coupled scenario       
        root architecture is not gowing  
        conductivities are not changing over time
        
    s                            soil model (RichardsWrapper(RichardsSP()))
    r                            xylem flux model (XylemFluxPython wrapping MappedSegments mapped to soil @param s)
    lookuptable_name             potentials a root soil interface    
    sim_time                     simulation time
    dt                           time step
    trans_f                      potential transpiration function 
    initial_age                       initial root system age  
    type_                        1 = water only, 2 = water and nitrate
    """

    wilting_point = -15000  # cm
    skip = 10  # for output and results, skip iteration
    max_iter = 10  # maximum for fix point iteration

    peri = PerirhizalPython(r.ms)
    peri.open_lookup(lookuptable_name)

    start_time = timeit.default_timer()

    psi_x_, psi_s_, sink_ , x_, y_, psi_s2_ = [], [], [], [], [], []  # for post processing
    soil_c_, c_ = [], []
    vol_ = [[], [], [], [], [], []]
    surf_ = [[], [], [], [], [], []]
    krs_ = []
    depth_ = []

    rs = r.ms
    nodes = rs.nodes
    segs = rs.segments
    ns = len(segs)
    mapping = rs.getSegmentMapper()  # because seg2cell is a dict

    for i in range(0, len(segs)):
        if segs[i].x == 0:
            collar_ind = i  # segment index of root collar
            break

    sx = s.getSolutionHead_()  # richards.py
    hsb = np.array([sx[j] for j in mapping])  # soil bulk matric potential per segment
    rsx = hsb.copy()  # initial values for fix point iteration
    rx = r.solve(initial_age, trans_f(0, dt), rsx, cells = False)
    rx_old = rx.copy()

    N = int(np.ceil(sim_time / dt))  # number of iterations

    print("Starting simulation loop")

    """ simulation loop """
    for i in range(0, N):

        t = i * dt  # current simulation time

        """ grow root system and update everything"""
        rs.simulate(dt, False)

        cell2seg = rs.cell2seg  # for debugging
        mapping = rs.getSegmentMapper()
        sx = s.getSolutionHead_()  # richards.py
        hsb = np.array([sx[j] for j in mapping])  # soil bulk matric potential per segment
        rsx = hsb.copy()  # initial values for fix point iteration

        cell_centers = s.getCellCenters_()
        cell_centers_z = np.array([cell_centers[j][2] for j in mapping])
        seg_centers_z = rs.getSegmentZ()

        outer_r = peri.get_outer_radii("length")
        inner_r = r.ms.radii
        types = r.ms.subTypes
        rho_ = np.divide(outer_r, np.array(inner_r))
        rho_ = np.minimum(rho_, np.ones(rho_.shape) * 200)  ############################################ (too keep within table)

        kr_ = r.params.getKr(initial_age + t)
        inner_kr_ = np.multiply(inner_r, kr_)  # multiply for table look up; here const
        inner_kr_ = np.maximum(inner_kr_, np.ones(inner_kr_.shape) * 1.e-7)  ############################################ (too keep within table)
        inner_kr_ = np.minimum(inner_kr_, np.ones(inner_kr_.shape) * 1.e-4)  ############################################ (too keep within table)

        wall_iteration = timeit.default_timer()
        wall_fixpoint = timeit.default_timer()

        err = 1.e6  # cm
        c = 0

        rx = r.solve(initial_age + t, trans_f(t, dt), rsx, False)
        rx_old = rx.copy()

        hsb_ = hsb - cell_centers_z  # from total matric potential to matric potential
        hsb_ = np.maximum(hsb_, np.ones(hsb_.shape) * -15000.)  ############################################ (too keep within table)
        hsb_ = np.minimum(hsb_, np.zeros(hsb_.shape))  ############################################ (too keep within table)

        while err > 1 and c < max_iter:

            """ interpolation """
            wall_interpolation = timeit.default_timer()
            rx_ = rx[1:] - seg_centers_z  # from total matric potential to matric potential
            rx_ = np.maximum(rx_, np.ones(rx_.shape) * -15000.)  ############################################ (too keep within table)

            # rsx = root_interface(rx_ , hsb_, inner_kr_, rho_, sra_table_lookup)
            rsx = peri.soil_root_interface_potentials(rx_ , hsb_, inner_kr_, rho_)

            rsx = rsx + seg_centers_z  # from matric potential to total matric potential
            wall_interpolation = timeit.default_timer() - wall_interpolation

            """ xylem matric potential """
            wall_xylem = timeit.default_timer()
            # print("Segment size from Python ", len(r.ms.segments), ns)
            rx = r.solve(initial_age + t, trans_f(initial_age + t, dt), rsx, False)  # xylem_flux.py, cells = False
            err = np.linalg.norm(rx - rx_old)
            wall_xylem = timeit.default_timer() - wall_xylem

            rx_old = rx.copy()
            c += 1

        wall_fixpoint = timeit.default_timer() - wall_fixpoint

        if type_ == 2:
            cc = s.getSolution_(1)  # kg/m3
            rsc = np.array([cc[i] for i in mapping])  # kg/m3
            seg_sol_fluxes = np.array(r.solute_fluxes(rsc))  # [g/day]
            soil_sol_fluxes = r.sumSegFluxes(seg_sol_fluxes)  # [g/day]
            # evap.add_nitrificatin_source(s, soil_sol_fluxes, nit_flux = 0.)  # = 1.14e-4 g/day
            evap.add_nitrificatin_source(s, soil_sol_fluxes, nit_flux = 0.)  # = 1.14e-4 g/day 1.e-7 * (75 * 16 * 1)
            s.setSource(soil_sol_fluxes.copy(), eq_idx = 1)  # [g/day], in moduels/richards.py

        wall_soil = timeit.default_timer()
        fluxes = r.segFluxes(initial_age + t, rx, rsx, approx = False, cells = False)
        collar_flux = r.collar_flux(initial_age + t, rx.copy(), rsx.copy(), k_soil = [], cells = False)  # validity checks
        err = np.linalg.norm(np.sum(fluxes) - collar_flux)
        if err > 1.e-6:
            print("error: summed root surface fluxes and root collar flux differ" , err, r.neumann_ind, collar_flux, np.sum(fluxes))
        err2 = np.linalg.norm(trans_f(initial_age + t, dt) - collar_flux)
        if r.last == "neumann":
            if err2 > 1.e-6:
                print("error: potential transpiration differs root collar flux in Neumann case" , err2)
        soil_fluxes = r.sumSegFluxes(fluxes)
        s.setSource(soil_fluxes.copy())  # richards.py
        s.solve(dt)

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
        #     max_b = [19, 2.5, 0.]
        #     cell_number = [1, 1, 200]
        #     vtk.plot_roots_and_soil(rs, "fluxes", fluxes.copy()[1:], s, True, min_b, max_b, cell_number, "nice_plot")
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
        x_.append(initial_age + t)  # day
        y_.append(np.sum(sink))  # cm3/day
        if type_ == 2:
            c_.append(-np.sum(seg_sol_fluxes))  # [g/day]

        if i % skip == 0:

            # if i % (24 * skip) == 0:
            print("time", initial_age + t, "{:g}/{:g} {:g} iterations".format(i, N, c), "wall times",
                  wall_interpolation / (wall_interpolation + wall_xylem), wall_xylem / (wall_interpolation + wall_xylem),
                  "number of segments", rs.getNumberOfSegments(), "root collar", rx[0])

            sink_.append(sink)  # cm3/day (per soil cell)

            psi_s2_.append(sx.copy())  # cm (per soil cell)

            if type_ == 2:
                soil_c_.append(cc)  # [kg/m3]

            ana = pb.SegmentAnalyser(r.ms.mappedSegments())  # VOLUME and SURFACE
            for j in range(0, 6):  # root types
                anac = pb.SegmentAnalyser(ana)
                anac.filter("subType", j)
                vol_[j].append(anac.getSummed("volume"))
                surf_[j].append(anac.getSummed("surface"))
            krs, _ = r.get_krs(initial_age + t, [collar_ind])
            krs_.append(krs)  # KRS
            depth_.append(ana.getMinBounds().z)

            """ direct vtp output """
            # psi_x_.append(rx.copy())  # cm (per root node)
            # psi_s_.append(rsx.copy())  # cm (per root segment)
            # ana.addData("rx", rx[1:])
            # ana.addData("rsx", rsx)
            # ana.addAge(initial_age + t)  # "age"
            # ana.addConductivities(r, initial_age + t)  # "kr", "kx"
            # ana.addFluxes(r, rx, rsx, initial_age + t)  # "axial_flux", "radial_flux"
            # ana.write("results/rs{0:05d}.vtp".format(int(i / skip)), ["radius", "subType", "creationTime", "organType", "rx", "rsx", "age", "kr", "kx", "axial_flux", "radial_flux"])

    print ("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s")

    return psi_x_, psi_s_, sink_, x_, y_, psi_s2_, vol_, surf_, krs_, depth_, soil_c_, c_

