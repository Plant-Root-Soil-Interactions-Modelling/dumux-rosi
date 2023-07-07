"""
functions for the cylindrical coupling approach
"""
import numpy as np
import timeit
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank(); max_rank = comm.Get_size()
import sys
from xylem_flux import *
import vtk_plot as vtk
import plantbox as pb
import evapotranspiration as evap
import timeit
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import fsolve
import visualisation.vtk_plot as vp
from pyevtk.hl import gridToVTK
import scenario_setup as scenario
import os 


def open_sri_lookup(filename):
    """ opens the look-up table from a file, to quickly find soil root interface potential """
    sri_table = np.load(filename + ".npy")
    x = np.load(filename + "_.npy", allow_pickle = True)
    rx_ = x[0]
    sx_ = x[1]
    inner_ = x[2]
    outer_ = x[3]
    return RegularGridInterpolator((rx_, sx_, inner_, outer_), sri_table)  # method = "nearest" fill_value = None , bounds_error=False


def soil_root_interface_table(rx, sx, inner_kr_, rho_, f):
    """
    finds potential at the soil root interface
        
    rx             xylem matric potential [cm]
    sx             bulk soil matric potential [cm]
    inner_kr       root radius times hydraulic conductivity [cm/day] 
    rho            geometry factor [1]
    f              function to look up the potentials
    """
    try:
        rsx = f((rx, sx, inner_kr_ , rho_))
    except:
        print("rx", np.min(rx), np.max(rx))  # 0, -16000
        print("sx", np.min(sx), np.max(sx))  # 0, -16000
        print("inner_kr", np.min(inner_kr_), np.max(inner_kr_))  # 1.e-7 - 1.e-4
        print("rho", np.min(rho_), np.max(rho_))  # 1. - 200.
    return rsx

def simulate_const(fname, s, rs, sri_table_lookup, sim_time, dt, trans_f, comp, rs_age, min_b, max_b, type):
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

    
    """ tabularized values for finding the zeros """
    if isinstance(sri_table_lookup, RegularGridInterpolator):
        root_interface = soil_root_interface_table  # function defined above
    else:
        raise
        root_interface = soil_root_interface  # function defined above
        
    
    
    wilting_point = -15000  # cm
    skip = 10  # 3 * 6  # for output and results, skip iteration
    split_type = 1  # type 0 == volume, type 1 == surface, type 2 == length
    max_iter = 1000

    """ 
    Initialize local soil models (around each root segment) 
    """
    start_time = timeit.default_timer()
    sx = s.getSolutionHead_()  # initial condition of soil [cm]
    cc = s.getSolution_(1)  # initial solute concentration [g/cm3]
    nodes = rs.rs.rs.nodes
    segs = rs.rs.rs.segments  # this is not nice (rs RhizoMappedSegments, rs.rs XylemFluxPython, rs.rs.rs MappedRootSystem(MappedSegments)
    
    seg2cell = rs.rs.rs.seg2cell
    cell2seg = rs.rs.rs.cell2seg
    mapping = rs.rs.rs.getSegmentMapper()
    hsb = np.array([sx[j] for j in mapping])  #soil bulk matric potential per segment
    cell_centers = s.getCellCenters_()
    cell_centers_z = np.array([cell_centers[j][2] for j in mapping])
    seg_centers_z = rs.getSegmentZ() 
    
    ns = len(segs)

    for i in range(0, len(segs)):
        if segs[i].x == 0:
            collar_ind = i  # segment index of root collar
            break

    dcyl = int(np.floor(ns / max_rank))
    if rank + 1 == max_rank:
        rs.initialize(s.soils[0], sx, np.array(range(rank * dcyl, ns)), cc)
        print ("\nInitialized final rank {:g}/{:g} [{:g}-{:g}] in {:g} s\n".format(rank + 1, max_rank, rank * dcyl, ns, timeit.default_timer() - start_time))
    else:
        rs.initialize(s.soils[0], sx, np.array(range(rank * dcyl, (rank + 1) * dcyl)), cc)
        print ("\nInitialized rank {:g}/{:g} [{:g}-{:g}] in {:g} s\n".format(rank + 1, max_rank, rank * dcyl, (rank + 1) * dcyl, timeit.default_timer() - start_time))

    r = rs.rs  # rename (XylemFluxPython)
    outer_r = r.rs.segOuterRadii()
    inner_r = r.rs.radii
    types = r.rs.subTypes
    rho_ = np.divide(outer_r, np.array(inner_r))
    rho_ = np.minimum(rho_, np.ones(rho_.shape) * 200)  

    psi_x_, psi_s_, sink_ , x_, y_, psi_s2_, soil_c_, c_, mass_soil_c = [], [], [], [], [], [], [], [], []
    vol_ = [[], [], [], [], [], []]
    surf_ = [[], [], [], [], [], []]
    krs_ = []
    depth_ = []
    # for post processing

    
    cell_volumes = s.getCellVolumes_()  # cm3
    net_flux = np.zeros(cell_volumes.shape)
    net_sol_flux = np.zeros(cell_volumes.shape)

    N = int(np.ceil(sim_time / dt))  # number of iterations

    """ simulation loop """
    for i in range(0, N):

        t = i * dt  # current simulation time
        
        sx = s.getSolutionHead_()  # [cm] richards.py
        hsb = np.array([sx[j] for j in mapping])  # soil bulk matric potential per segment
        rsx = hsb.copy()  # initial values for fix point iteration
        
        kr_ = r.getKr(rs_age + t)
        inner_kr_ = np.multiply(inner_r, kr_)  # multiply for table look up; const
        inner_kr_ = np.maximum(inner_kr_, np.ones(inner_kr_.shape) * 1.e-7)  
        inner_kr_ = np.minimum(inner_kr_, np.ones(inner_kr_.shape) * 1.e-4)  
        
        err = 1.e6  # cm
        c = 0
        
        rx = r.solve(rs_age, trans_f(rs_age + 0., dt), 0., rsx, cells = False, wilting_point = wilting_point, soil_k = [])
        rx_old = rx.copy()
        
        hsb_ = hsb - cell_centers_z  # from total matric potential to matric potential
        hsb_ = np.maximum(hsb_, np.ones(hsb_.shape) * -15000.)
        hsb_ = np.minimum(hsb_, np.zeros(hsb_.shape))  

        while err > 1 and c < max_iter:

            """ interpolation """
            wall_interpolation = timeit.default_timer()
            rx_ = rx[1:] - seg_centers_z  # from total matric potential to matric potential
            rx_ = np.maximum(rx_, np.ones(rx_.shape) * -15000.)  ############################################ (too keep within table)
            rsx = root_interface(rx_ , hsb_, inner_kr_, rho_, sri_table_lookup)
            rsx = rsx + seg_centers_z  # from matric potential to total matric potential

            """ xylem matric potential """
            wall_xylem = timeit.default_timer()
            # print("Segment size from Python ", len(r.rs.segments), ns)
            rx = r.solve(rs_age, trans_f(rs_age + t, dt), 0., rsx, cells = False, wilting_point = wilting_point, soil_k = [])
            err = np.linalg.norm(rx - rx_old)
            wall_xylem = timeit.default_timer() - wall_xylem

            rx_old = rx.copy()
            c += 1
            print('number of iterations', c, err) 

        #get current exudation rate (g/day/root apex) for comp
        kexu = scenario.exudation_rates(rs_age+t, comp) #[kg/(m2 day)]
        sol_result = r.exudate_fluxes(rs_age+t, kexu) #[g/day/seg], [g/cm2/day/seg]
        
        if rank == 0:
            seg_fluxes = np.array(r.segFluxes(rs_age + t, rx, rsx, False, False, []))
            seg_sol_fluxes = np.array(sol_result[0])  # [g/day/seg]
        else:
            seg_fluxes = None
            seg_sol_fluxes = None

        seg_fluxes = comm.bcast(seg_fluxes, root = 0)
        seg_sol_fluxes = comm.bcast(seg_sol_fluxes, root = 0)

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
        wall_xylem = timeit.default_timer() - wall_xylem

        """ 2. local soil models """
        wall_local = timeit.default_timer()
        if rank == 0:
            for key, value in cell2seg.items():  # check cell2seg
                if key < 0: #part of the root system is above ground and out of the soild omain 
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


            proposed_outer_fluxes = r.splitSoilFluxes(net_flux / dt, split_type) 
            proposed_outer_sol_fluxes = r.splitSoilFluxes(net_sol_flux / dt, split_type)
            # if this fails, a segment is not mapped, i.e. out of soil domain
        else:
            proposed_outer_fluxes = None
            proposed_outer_sol_fluxes = None
        proposed_outer_fluxes = comm.bcast(proposed_outer_fluxes, root = 0)
        proposed_outer_sol_fluxes = comm.bcast(proposed_outer_sol_fluxes, root = 0)
        
        rs.solve(dt, rsx, proposed_outer_fluxes,sol_result[1],proposed_outer_sol_fluxes) #or 0?
        # water: left dirchlet, right neumann; solute: left and right Neumann<----

        wall_local = timeit.default_timer() - wall_local
        
        """ 3a. macroscopic soil model """
        wall_macro = timeit.default_timer()

        water_content = np.array(s.getWaterContent_())  # theta per cell [1]
        soil_water = np.multiply(water_content, cell_volumes)  # water per cell [cm3]
        solute_conc = np.array(s.getSolution_(1)) #g/cm3
        soil_solute = np.multiply(solute_conc, soil_water) # [g] 

        soil_fluxes = r.sumSegFluxes(seg_fluxes)  # [cm3/day]  per soil cell
        soil_sol_fluxes = r.sumSegFluxes(seg_sol_fluxes)  # [g/day]
        soil_sol_fluxes = evap.decay(soil_sol_fluxes, dt, s.decay)  #[g/day]
        s.setSource(soil_fluxes.copy(), eq_idx = 0)  # [cm3/day], in modules/richards.py
        s.setSource(soil_sol_fluxes.copy(), eq_idx = 1)  # [g/day], in modules/richards.py
        s.solve(dt)  # in modules/solverbase.py

        wall_macro = timeit.default_timer() - wall_macro
        # if rank == 0:
        #     print("]", end = "")

        """ 3b. calculate net fluxes """
        wall_netfluxes = timeit.default_timer()
        water_content = np.array(s.getWaterContent_())
        new_soil_water = np.multiply(water_content, cell_volumes)  # [cm3], calculate net flux
        net_flux = new_soil_water - soil_water  # change in water per cell [cm3]
        for k, root_flux in soil_fluxes.items():
            net_flux[k] -= root_flux * dt

        """ 3c. calculate mass net fluxes """
        solute_conc = np.array(s.getSolution_(1)) #g/cm3 
        new_soil_solute = np.multiply(solute_conc, new_soil_water) #[g] 

        net_sol_flux = new_soil_solute - soil_solute  # change in solute per cell [g]
        for k, root_sol_flux in soil_sol_fluxes.items():
            net_sol_flux[k] += root_sol_flux * dt
        
        soil_water = new_soil_water #[cm3]
        soil_solute = new_soil_solute #[g]

        wall_netfluxes = timeit.default_timer() - wall_netfluxes

        """ remember results ... """
        sx = s.getSolutionHead_() #[cm]
        cc = s.getSolution_(1)  # [g/cm3]
        if rank == 0:
            sink = np.zeros(sx.shape)
            for k, v in soil_fluxes.items():
                sink[k] += v
            sink_.append(sink)  # cm3/day (per soil cell)
            x_.append(rs_age + t)  # day
            y_.append(np.sum(sink))  # cm3/day
            c_.append(-np.sum(seg_sol_fluxes))  # [g/day]
            psi_s2_.append(sx)  # cm (per soil cell)
            soil_c_.append(cc)  # g/cm3
            mass_soil_c.append(np.sum(new_soil_solute)) #[g]

            ana = pb.SegmentAnalyser(r.rs.mappedSegments())  # VOLUME and SURFACE
            for j in range(0, 6):  # root types
                anac = pb.SegmentAnalyser(ana)
                anac.filter("subType", j)
                vol_[j].append(anac.getSummed("volume"))
                surf_[j].append(anac.getSummed("surface"))
            krs, _ = r.get_krs(rs_age + t, [collar_ind])
            krs_.append(krs)  # KRS
            depth_.append(ana.getMinBounds().z)

            
        if i % skip == 0 and rank == 0:
            print("{:g}/{:g} iterations".format(i, N), "time", rs_age + t, "wall times",
                  wall_xylem / (wall_xylem + wall_local + wall_macro + wall_netfluxes),
                  wall_local / (wall_xylem + wall_local + wall_macro + wall_netfluxes),
                  wall_macro / (wall_xylem + wall_local + wall_macro + wall_netfluxes),
                  wall_netfluxes / (wall_xylem + wall_local + wall_macro + wall_netfluxes),
                  "segments:", rs.getNumberOfMappedSegments(), "root collar:", rx[0], "\n")
            print("time", rs_age + t, "rsx", np.min(rsx), np.max(rsx), "ccx", np.min(cc), np.max(cc), "trans", trans_f(rs_age + t, dt))
            psi_x_.append(rx.copy())  # cm (per root node)
            psi_s_.append(rsx.copy())  # cm (per root segment)


        #print('current time is',i/240)
        #print('current root age', rs_age+i/240) 
        if (rs_age % 10 == 0) and (i/240 == 0.5):  # every 2.5 days i % (20 * 12) == 0
            # map concentration cylinders to 1mm grid in all directions 
            xx_ = np.linspace(min_b[0], max_b[0], int(10*(max_b[0]-min_b[0]))) 
            yy_ = np.linspace(min_b[1], max_b[1], int(10*(max_b[1]-min_b[1])))
            zz_ = np.linspace(min_b[2],max_b[2], int(10*(max_b[2]-min_b[2])))
            XX, YY, ZZ = np.meshgrid(xx_, yy_, zz_, indexing='ij')
            C  = rs.map_cylinders_solute(XX,YY,ZZ,'cyl_exu') #[g/cm3]

            if not os.path.exists('results/vts_'+fname):
                os.makedirs('results/vts_'+fname)
                os.makedirs('results/concentration_'+fname)
            
            gridToVTK("results/vts_"+fname+"/./exudate_day"+('{:03d}'.format(int(rs_age))), XX, YY, ZZ, pointData = {"Exudate concentration (g/cm3)":C}) 
            np.save("results/concentration_"+fname+"/day"+("{:03d}".format(int(rs_age))), C) 

    if rank == 0:
        print ("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s")

    return psi_x_, psi_s_, sink_, x_, y_, psi_s2_, vol_, surf_, krs_, depth_, soil_c_, c_, mass_soil_c

