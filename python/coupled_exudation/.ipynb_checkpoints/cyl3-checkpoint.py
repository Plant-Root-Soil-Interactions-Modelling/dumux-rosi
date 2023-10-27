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
import visualisation.vtk_plot as vp
from scenario_setup import weather


def exudate_fluxes(exud ):
    """ input [mol/d, returns [g/day]
    
    """
    #if not isinstance(kex,int):
    #    ages1 = self.get_ages(rs_age)

    #types = np.asarray(self.rs.subTypes, int)
    #a = self.rs.radii
    #l = self.rs.segLength()
    #sf = np.zeros(len(l),)
    #sf = np.array([2 * np.pi * a[i] * l[i] * 1.e-4 * kex for i in range(len(l))])
    #notRoots = np.where(organtypes != 2)[0]
    #sf[notRoots] = 0
    
    #take out seed node to have it per segments
    return np.array(exud)[1:] * 342.3# g/mol  # kg/day -> g/day

def simulate_const(s, rs, sim_time, dt, kexu, rs_age,repartition, type, Q_Exud):
    """     
    simulates the coupled scenario       
        root architecture is not growing  
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
    weatherX = weather(sim_time+rs_age) 
    """ 
    Initialize local soil models (around each root segment) 
    """
    start_time = timeit.default_timer()
    sx = s.getSolutionHead_()  # initial condition of soil [cm]
    cc = s.getSolution_(1)  # solute concentration [kg/m3].
    
    
    cell2segVals = np.concatenate((list(rs.cell2seg.values()))).flatten()
    if len(cell2segVals) != len(set(cell2segVals)):
        print(rs.seg2cell)
        print(rs.cell2seg)
        print(cell2segVals)
        print(len(cell2segVals), len(set(cell2segVals)))
        raise Exception
    
    
    repartitionOld = repartition
    nsOld = sum(repartitionOld)
    ns = len(rs.segments)
    dcyl = int(np.floor(ns / max_rank))
    repartition = np.array([dcyl for i in range(max_rank)])
    repartition[max_rank -1] = ns - rank * dcyl
    toAdd = repartition - repartitionOld
    
    if toAdd[rank] > 0:
        if len(sx.shape)!=1:
            raise Exception
        if rank == 0:
            rs.update( sx, np.array([i for i in range(toAdd[rank])]) + nsOld, cc )
            print ("Initialized rank {:g}/{:g} [{:g}-{:g}] in {:g} s".format(rank + 1, max_rank, nsOld, toAdd[rank] + nsOld, timeit.default_timer() - start_time))
        else:
            rs.update( sx,  np.array([i for i in range(toAdd[rank-1],toAdd[rank])]) + nsOld, cc )
            print ("Initialized rank {:g}/{:g} [{:g}-{:g}] in {:g} s".format(rank + 1, max_rank, toAdd[rank-1] + nsOld, toAdd[rank] + nsOld, timeit.default_timer() - start_time))
        
        
    

    segs = rs.segments  
    Nt = len(rs.nodes)
    print(len(segs), Nt)
    assert len(segs) == (Nt -1)
    seg2cell = rs.seg2cell
    cell2seg = rs.cell2seg
    
    ns = len(segs)
    r = rs.rs  # rename (XylemFluxPython)
    for i in range(0, len(segs)):#sure about that, even in case of artificial shoot?
        if segs[i].x == 0:
            collar_ind = i  # segment index of root collar
            break

    psi_x_, psi_s_, sink_ , x_, y_, psi_s2_, soil_c_, c_ , c_All= [], [], [], [], [], [], [], [] ,[]
    vol_ = [[], [], [], [], [], []]
    surf_ = [[], [], [], [], [], []]
    krs_ = []
    depth_ = []
    # for post processing

    
    cell_volumes = s.getCellVolumes_()  # cm3
    net_flux = np.zeros(cell_volumes.shape)
    net_sol_flux = np.zeros(cell_volumes.shape)

    N = 1#int(np.ceil(sim_time / dt))  # number of iterations

    """ simualtion loop """
    for i in range(0, N):

        t = i * dt  # current simulation time
        
        weatherX = weather(t)

        """ 1. xylem model """
        # if rank == 0:
        #     print("[x", end = "")
        wall_xylem = timeit.default_timer()

    
        ###
        rsx = rs.get_inner_heads(weather=weatherX)  # matric potential at the root soil interface, i.e. inner values of the cylindric models (not extrapolation to the interface!) [cm]
        if type == "dumux_dirichlet":
            rsc = [cc[seg2cell[i]] for i in range(0, ns)]  # kg/m3
        else:
            rsc = rs.get_inner_solutes(1)  # kg/m3
        comm.barrier()
        
        subtypes = np.asarray(rs.subTypes, int)
        organTypes = np.asarray(rs.organTypes, int)
        a = np.array(rs.radii)
        l = np.array(rs.segLength())
        
        if rank == 0:
            print("need to add fixed point iteration for water and carbon fluxes between plant-rhizosphere")
            print("\nINITIAL root-soil interface matric potentials", np.nanmin(rsx), np.nanmax(rsx), np.nanmin(rsc), np.nanmax(rsc))
            
            
            r.minLoop = 1000
            r.maxLoop = 5000
            r.Qlight = weatherX["Qlight"]
            r.solve_photosynthesis(sim_time_ = sim_time+rs_age, 
                        sxx_=rsx, #will not use the output of the air bc. so it s just a stand-in for now
                        cells_ = False,
                        ea_ = weatherX["ea"],#not used
                        es_=weatherX["es"],#not used
                        verbose_ = False, doLog_ = False,
                        TairC_= weatherX["TairC"],#not used
                        outputDir_= "./results/rhizoplantExud")
            rx = np.array(r.psiXyl)
            seg_fluxes = np.array(r.outputFlux)# [cm3/day]
            errLeuning = sum(seg_fluxes)
            
            seg_sol_fluxes = exudate_fluxes(Q_Exud)#g/day for segments
            #r.exudate_fluxes(rs_age+1, kexu))  # [g/day]
        else:
            rx = None
            seg_fluxes = None
            seg_sol_fluxes = None
        rx = comm.bcast(rx, root = 0)
        seg_fluxes = comm.bcast(seg_fluxes, root = 0)
        seg_sol_fluxes = comm.bcast(seg_sol_fluxes, root = 0)
        #print("rx",rx)
        leavesSegs = np.where(np.array(rs.organTypes )==4)
        
        """ some checks """
        
        assert np.abs(np.sum(seg_fluxes)) < 1e-10 #no change in water storage
        fluxes_leaves = seg_fluxes[leavesSegs]
        assert (np.abs(sum(fluxes_leaves))- np.abs(sum(np.array(r.Ev)))) < 1e-10
        if (min(r.Ev) < 0) or (min(r.Jw) < 0) or (min(fluxes_leaves)<0):
            print("leaf looses water", min(r.Ev),min(r.Jw), min(fluxes_leaves))
            raise Exception
            
        #does not use rsx for the leaves, but the pg (guard cells wat. pot.) values
        #[cm3/day]
        proposed_inner_fluxes = r.segFluxes(t, rx.copy(), rsx.copy(), approx = False, cells = False, soil_k = [])  # [cm3/day]
        assert (np.array(proposed_inner_fluxes) == seg_fluxes).all()
        ###


        """ 2. local soil models """
        wall_local = timeit.default_timer()
        if rank == 0:
                        
            proposed_outer_fluxes = r.splitSoilFluxes(net_flux / dt, split_type) 
            proposed_outer_sol_fluxes = r.splitSoilFluxes(net_sol_flux / dt, split_type)
            #print("proposed_outer_fluxes",sum(proposed_outer_fluxes),sum(net_flux / dt))
            #print("proposed_outer_sol_fluxes",sum(proposed_outer_sol_fluxes),sum(net_sol_flux / dt))
            # if this fails, a segment is not mapped, i.e. out of soil domain
        else:
            proposed_outer_fluxes = None
            proposed_outer_sol_fluxes = None
        proposed_outer_fluxes = comm.bcast(proposed_outer_fluxes, root = 0)
        proposed_outer_sol_fluxes = comm.bcast(proposed_outer_sol_fluxes, root = 0)
        seg_rx = np.array([0.5 * (rx[seg.x] + rx[seg.y]) for seg in segs])
        
        #ages = XylemFluxPython.get_ages(r,rs_age+1)
        subTypes = rs.subTypes
        
        
        kex = np.where(np.array(rs.organTypes)==2, 1e-2,0)
                
        rs.solve(dt, seg_rx, proposed_outer_fluxes,kex,proposed_outer_sol_fluxes) #or 0?
        # water: left dirchlet, right neumann; solute: left and right Neumann<----

        wall_local = timeit.default_timer() - wall_local
        
        """ 3a. macroscopic soil model """
        wall_macro = timeit.default_timer()

        water_content = np.array(s.getWaterContent_())  # theta per cell [1]
        soil_water = np.multiply(water_content, cell_volumes)  # water per cell [cm3]
        solute_conc = np.array(s.getSolution_(1))
        soil_solute = np.multiply(solute_conc, soil_water)

        soil_fluxes = r.sumSegFluxes(seg_fluxes)  # [cm3/day]  per soil cell
        soil_sol_fluxes = r.sumSegFluxes(seg_sol_fluxes)  # [g/day]
        print("seg_sol_fluxes1",seg_sol_fluxes,soil_sol_fluxes)
        soil_sol_fluxes = evap.decay(soil_sol_fluxes, dt, s.decay)  #[g/day]
        print("seg_sol_fluxes2",soil_sol_fluxes)
        s.setSource(soil_fluxes.copy(), eq_idx = 0)  # [cm3/day], in modules/richards.py
        s.setSource(soil_sol_fluxes.copy(), eq_idx = 1)  # [g/day], in modules/richards.py
        s.solve(dt)  # in modules/solverbase.py

        wall_macro = timeit.default_timer() - wall_macro
        # if rank == 0:
        #     print("]", end = "")

        """ 3b. calculate net fluxes """
        wall_netfluxes = timeit.default_timer()
        water_content = np.array(s.getWaterContent_())
        new_soil_water = np.multiply(water_content, cell_volumes)  # calculate net flux
        net_flux = new_soil_water - soil_water  # change in water per cell [cm3]
        #print('net_flux', net_flux) 
        for k, root_flux in soil_fluxes.items():
            net_flux[k] -= root_flux * dt

        """ 3c. calculate mass net fluxes """
        solute_conc = np.array(s.getSolution_(1))
        print("solute_conc",solute_conc)
        new_soil_solute = np.multiply(solute_conc, soil_water)
        print("new_soil_solute",new_soil_solute)
        #raise Exception

        net_sol_flux = new_soil_solute - soil_solute  # change in water per cell [cm3]
        #print('net_sol_flux', net_sol_flux) 
        for k, root_sol_flux in soil_sol_fluxes.items():
            net_sol_flux[k] += root_sol_flux * dt
        
        soil_water = new_soil_water
        soil_solute = new_soil_solute

        wall_netfluxes = timeit.default_timer() - wall_netfluxes

        """ remember results ... """
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
            c_All.append(seg_sol_fluxes)  # [cm3/day]
            psi_s2_.append(sx)  # cm (per soil cell)
            
            cutoff = 1e-15 #is get value too small, makes paraview crash
            cc_p = np.array(cc)
            cc_p[abs(cc_p) < cutoff] = 0
            
            soil_c_.append(cc_p)  # [kg/m3]

            ana = pb.SegmentAnalyser(r.rs.mappedSegments())  # VOLUME and SURFACE
            for j in range(0, 6):  # root types
                anac = pb.SegmentAnalyser(ana)
                anac.filter("subType", j)
                vol_[j].append(anac.getSummed("volume"))
                surf_[j].append(anac.getSummed("surface"))
            #krs, _ = r.get_krs(rs_age + t, [collar_ind])
            #krs_.append(krs)  # KRS
            depth_.append(ana.getMinBounds().z)

            
        if i % skip == 0 and rank == 0:
            print("{:g}/{:g} iterations".format(i, N), "time", rs_age + t, "wall times",
                  wall_xylem / (wall_xylem + wall_local + wall_macro + wall_netfluxes),
                  wall_local / (wall_xylem + wall_local + wall_macro + wall_netfluxes),
                  wall_macro / (wall_xylem + wall_local + wall_macro + wall_netfluxes),
                  wall_netfluxes / (wall_xylem + wall_local + wall_macro + wall_netfluxes),
                  "segments:", rs.getNumberOfMappedSegments(), "root collar:", rx[0], "\n")
            print("time", rs_age + t, "rsx", np.min(rsx), np.max(rsx), "ccx", np.min(cc), np.max(cc), "trans",sum(np.array(r.Ev)))
            psi_x_.append(rx.copy())  # cm (per root node)
            psi_s_.append(rsx.copy())  # cm (per root segment)

        # if rank == 0:
        #     print("]")

    if rank == 0:
        print ("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s")

    return psi_x_, psi_s_, sink_, x_, y_, psi_s2_, vol_, surf_,  depth_, soil_c_, c_, repartition, c_All #krs_,

