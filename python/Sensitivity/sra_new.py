"""
    Dynamic:

    Dynamic water uptake simulation taking nonlinear rhizosphere resistance into account 
    (called by run_sra.py)    
    
    Daniel Leitner, 2025    
"""
import numpy as np
import timeit

import plantbox as pb
import visualisation.vtk_plot as vtk
import functional.van_genuchten as vg
from functional.Perirhizal import PerirhizalPython

import evapotranspiration as evap
import carbon_cost

def simulate_dynamic(s, r, lookuptable_name, sim_time, dt, trans_f, output_times, initial_age = 1., cc_data = None):
    """     
    simulates the coupled scenario       
        
    s                            soil model (RichardsWrapper(RichardsSP()))
    r                            PlantHydraulicModel 
    lookuptable_name             matric potentials at the root soil interface (precomputed in a 4D table)    
    sim_time                     simulation time
    dt                           time step
    trans_f                      potential transpiration function 
    output_times                 time points for additional vtp outputs
    initial_age                  initial root system age  
    cc_data                      (optional) dict with "radii" and "anatomy" for carbon model 2 and 3   
    
    return:

    List 1 (r1):
    times_                       times (simulation time plus initial age)
    pot_trans_                   potential transpiration [cm3/day]
    act_trans_                   actual transpiration [cm3/day]
    collar_pot_                  root collar potential [cm]  

    List 2 (r2): for each 10th time step
    times_lr_                    times for the following arrays  
    sink_                        water uptake (per cell)  
    psi_s_                       soil matric potentials (per cell)   
    net_change_                  net domain water change (including plant uptake)
    len_                         root system length [cm] per subType
    surf_                        root system surface [cm2] per subType
    vol_                         root system volume [cm3] per subType
    depth_                       root system depth [cm] 
    krs_                         root system hydraulic conductivity [cm2/day]
    carbon                       carbon cost [3 models]

    List 3 (r3): of defined output times
    out_times_                   output times including start and final simulation time (plus initial age)
    ana_                         list of segment analysers at the out_times 
    """

    wilting_point = -15000  # cm
    skip = 10  # for output and results, skip iteration (TODO)
    matimes_iter = 10  # maximum for fix point iteration

    print("simulate_dynamic starting" , flush = True)
    peri = PerirhizalPython(r.ms)
    peri.open_lookup("data/" + lookuptable_name)
    print("opening look up table:", "data/" + lookuptable_name, flush = True)

    start_time = timeit.default_timer()

    times_, pot_trans_, act_trans_, collar_pot_, = [], [], [], [] # for all time steps
    times_lr_, sink_, psi_s_ = [],[], []
    net_change_, depth_, krs_ , carbon_  = [], [], [], []    
    len_ = [[], [], [], [], [], []]
    surf_ = [[], [], [], [], [], []]    
    vol_ = [[], [], [], [], [], []]
    ana_ = []
   
    rs = r.ms
    nodes = rs.nodes
    segs = rs.segments

    sx = s.getSolutionHead_()  # richards.py
    hsb_ = np.array(rs.getHs(sx))  # matric potential per segment
    rsx = hsb_.copy()  # initial values for fix point iteration
    rx = r.solve(initial_age, trans_f(0, dt), rsx, cells = False)
    rx_ = rx.copy()

    N = int(np.ceil(sim_time / dt))  # number of iterations
    output_time_indices = np.round(np.array(output_times)/dt).astype(int)
    output_time_indices = np.append(output_time_indices, N-1)
    output_time_indices = np.insert(output_time_indices, 0, 0)
    out_times_ = output_times.copy()
    out_times_.append(sim_time)
    out_times_.insert(0,0.)
    print("Starting simulation loop",N, "iterations", flush = True)
    
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
        inner_r = np.array([rs.getEffectiveRadius(i) for i in range(0, n)])  # rs.radii
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
            rx = np.maximum(rx, np.ones(rx.shape) * -15000.)
            rx = np.minimum(rx, np.zeros(rx.shape))
            rsx = peri.soil_root_interface_potentials(rx[1:], hsb_, inner_kr_, rho_)
            wall_interpolation = timeit.default_timer() - wall_interpolation

            """ xylem matric potential """
            wall_xylem = timeit.default_timer()
            rx = r.solve_again(initial_age + t, trans_f(t, dt), rsx, cells = False)
            res_ = rx - rx_
            err_ =  np.linalg.norm(res_)/ np.sqrt(rx.shape[0])  # rmse            
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
        net_change_.append((s.getWaterVolume() - water_)/dt) # cm3/day 

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

        """ remember results ... """
        times_.append(initial_age + t)  # day
        pot_trans_.append(trans_f(t, dt))  # cm3/day
        act_trans_.append(np.sum(fluxes))  # cm3/day
        collar_pot_.append(rx[0])

        if i % skip == 0:
            
            # print("age {:g}".format(initial_age + t), "{:g}/{:g} {:g} iterations, rmse {:g}".format(i, N, c, err_), "; wall times {:g} {:g}".format(
            #     wall_interpolation / (wall_interpolation + wall_xylem), wall_xylem / (wall_interpolation + wall_xylem)),
            #       "number of segments", rs.getNumberOfSegments(), "collar potential {:g}".format(rx[0]), flush = True)
            # print("res_", np.argmax(np.abs(res_)), np.max(np.abs(res_)), err_)
            
            times_lr_.append(initial_age + t)
            
            sink = np.zeros(sx.shape)
            for k, v in soil_fluxes.items():
                sink[k] += v
            sink_.append(sink)  # cm3/day (per soil cell)
            psi_s_.append(sx.copy())  # cm (per soil cell)           
            
            ana = pb.SegmentAnalyser(rs.mappedSegments())  # VOLUME and SURFACE
            for j in range(0, 6):  # root types
                anac = pb.SegmentAnalyser(ana)
                anac.filter("subType", j)
                len_[j].append(anac.getSummed("length"))
                surf_[j].append(anac.getSummed("surface"))
                vol_[j].append(anac.getSummed("volume"))
            
            vol = ana.getParameter("volume")
            c1 = carbon_cost.carbon_cost_volume(vol)
            if cc_data: # untested
                subType = ana.getParameter("subType")
                segLeng = ana.getParameter("length") 
                r = cc_data["radii"]
                a = cc_data["anatomy"]
                c2 = carbon_cost.carbon_cost_simple(vol, subType, a[0][0], a[1][0], a[2][0])
                c3 = carbon_cost.carbon_cost_anatomy(vol, subType, segLen, a[0], a[1], a[2], r[0], r[1], r[2])
                carbon_.append([c1,c2,c3])
            else:
                carbon_.append([c1])    
                       
            depth_.append(rs.getMinBounds().z)                
            krs, _ = r.get_krs(initial_age + t)
            krs_.append(krs)  # KRS

        """ direct vtp output """
        if i in output_time_indices:
            n = len(rs.radii)
            radii = np.array([rs.getEffectiveRadius(i) for i in range(0, n)])
            old_radii = rs.radii
            rs.radii = radii
            ana = pb.SegmentAnalyser(rs.mappedSegments())  # VOLUME and SURFACE
            rs.radii = old_radii
            ana.addData("rx", rx[1:])
            ana.addData("rsx", rsx)
            ana.addAge(initial_age + t)  # "age"
            ana.addHydraulicConductivities(r.params, initial_age + t)  # "kr", "kx"
            ana.addFluxes(r, rx[1:], rsx, initial_age + t)  # "axial_flux", "radial_flux"
            ana_.append(ana)

    print ("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s")
    
    out_times_ = np.array(out_times_)+np.ones(output_time_indices.shape)*initial_age

    r1 = [times_, pot_trans_, act_trans_, collar_pot_] 
    r2 = [times_lr_, sink_, psi_s_, net_change_, len_, surf_, vol_,  depth_, krs_, carbon_]
    r3 = [out_times_, ana_]
    
    return r1, r2, r3
    """ 
    List 1 (r1):
    times_                       times (simulation time plus initial age)
    pot_trans_                   potential transpiration [cm3/day]
    act_trans_                   actual transpiration [cm3/day]
    collar_pot_                  root collar potential [cm]  

    List 2 (r2): for each 10th time step
    times_lr_                    times for the following arrays  
    sink_                        water uptake (per cell)  
    psi_s_                       soil matric potentials (per cell)   
    net_change_                  net domain water change (including plant uptake)
    len_                         root system length [cm] per subType
    surf_                        root system surface [cm2] per subType
    vol_                         root system volume [cm3] per subType
    depth_                       root system depth [cm] 
    krs_                         root system hydraulic conductivity [cm2/day]
    carbon                       carbon cost [3 models]

    List 3 (r3): of defined output times
    out_times_                   output times including start and final simulation time (plus initial age)
    ana_                         list of segment analysers at the out_times 
    """
