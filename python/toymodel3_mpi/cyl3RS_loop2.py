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
from decimal import *

from scenario_setup import write_file_array

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

def simulate_const(s, rs, sim_time, dt, kexu, rs_age,repartition, type, Q_Exud,
                    trans_maize,
                    plantType = "RS", r = [], wilting_point =-15000, 
                    net_sol_flux=[], net_flux = []):
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
    #wilting_point = -15000  # cm
    skip = 10  # 3 * 6  # for output and results, skip iteration
    split_type = 1  # type 0 == volume, type 1 == surface, type 2 == length
    weatherX = weather(rs_age) 
    """ 
    Initialize local soil models (around each root segment) 
    """
    start_time = timeit.default_timer()
    sx = s.getSolutionHead_()  # initial condition of soil [cm]
    cc = s.getSolution_(1)  # solute concentration [kg/m3].
    
    # for i in range(0, len(segs)):
        # if segs[i].x == 0:
            # collar_ind = i  # segment index of root collar
            # break


    
    cell2segVals = np.concatenate((list(rs.rs.cell2seg.values()))).flatten()
    if len(cell2segVals) != len(set(cell2segVals)):#make sure all values are unique
        print(rs.rs.seg2cell)
        print(rs.rs.cell2seg)
        print(cell2segVals)
        print(len(cell2segVals), len(set(cell2segVals)))
        raise Exception
    
    
    repartitionOld = repartition
    nsOld = sum(repartitionOld)
    ns = len(rs.get_segments())
    dcyl = int(np.floor(ns / max_rank))
    repartition = np.array([dcyl for i in range(max_rank)])
    repartition[max_rank -1] = ns - rank * dcyl
    toAdd = repartition - repartitionOld
    
    if toAdd[rank] > 0: # that could be put in main file (after RS growth)
        if len(sx.shape)!=1:
            raise Exception
        if rank == 0:
            r.update( sx, np.array([i for i in range(toAdd[rank])]) + nsOld)#, cc )
            print ("Initialized rank {:g}/{:g} [{:g}-{:g}] in {:g} s".format(rank + 1, max_rank, nsOld, toAdd[rank] + nsOld, timeit.default_timer() - start_time))
        else:
            r.update( sx,  np.array([i for i in range(toAdd[rank-1],toAdd[rank])]) + nsOld)#, cc )
            print ("Initialized rank {:g}/{:g} [{:g}-{:g}] in {:g} s".format(rank + 1, max_rank, toAdd[rank-1] + nsOld, toAdd[rank] + nsOld, timeit.default_timer() - start_time))
    # print("repartitionOld", repartitionOld, len(repartitionOld))
    if sum ( repartitionOld)==0: 
        print("first check")    
        r.checkMassOMoleBalance()

    segs = rs.rs.segments
    Nt = len(rs.rs.nodes)
    # print(len(segs), Nt)
    assert len(segs) == (Nt -1)
    seg2cell = rs.rs.seg2cell
    cell2seg = rs.rs.cell2seg
    
    ns = len(segs)
    #r = rs.rs  # rename (XylemFluxPython)

    psi_x_, psi_s_, sink_ , x_, y_, psi_s2_, soil_c_, c_ , c_All= [], [], [], [], [], [], [], [] ,[]
    vol_ = [[], [], [], [], [], []]
    surf_ = [[], [], [], [], [], []]
    krs_ = []
    depth_ = []
    # for post processing

    
    cell_volumes = s.getCellVolumes_()  # cm3
    #net_flux = np.zeros(cell_volumes.shape)
    #net_sol_flux = np.zeros(cell_volumes.shape)

    N = int(np.ceil(sim_time / dt))  # number of iterations
    """ simualtion loop """
    for i in range(0, N):

        t = i * dt  # current simulation time
        
        weatherX = weather(t)

        """ 1. xylem model """
        # if rank == 0:
        #     print("[x", end = "")
        wall_xylem = timeit.default_timer()

    
        ###
        rsx = r.get_inner_heads(weather=weatherX)  # matric potential at the root soil interface, i.e. inner values of the cylindric models (not extrapolation to the interface!) [cm]
        if type.startswith("dumux"):#that s not used
            rsc = [cc[seg2cell[i]] for i in range(0, ns)]  # kg/m3
        else:
            rsc = r.get_inner_solutes(1)  # kg/m3
        comm.barrier()
        
        subtypes = np.asarray(rs.rs.subTypes, int)
        organTypes = np.asarray(rs.rs.organTypes, int)
        a = np.array(rs.rs.radii)
        l = np.array(rs.rs.segLength())
        
        if rank == 0:
            # print("need to add fixed point iteration for water and carbon fluxes between plant-rhizosphere")
            # print("\nINITIAL root-soil interface matric potentials", np.nanmin(rsx), np.nanmax(rsx), np.nanmin(rsc), np.nanmax(rsc))
            #print("transpiration tot:", trans_maize(rs_age + 0., dt)*10)
            rx = rs.solve(rs_age, 0, 0., rsx, cells = False, 
                                wilting_point = wilting_point, soil_k = [])
            #trans_maize(rs_age + 0., dt)*10
            rx = np.array(rx)
            seg_fluxes =np.full(len( np.array(rs.segFluxes(rs_age + t, rx, rsx, False, False, []))),0.)# [cm3/day]
            #seg_fluxes = np.array(rs.segFluxes(rs_age + t, rx, rsx, False, False, []))# [cm3/day]
            # print("rx, rsx", rx, rsx)
            
            seg_sol_fluxes = Q_Exud # exudate_fluxes(Q_Exud)#g/day for segments
            #r.exudate_fluxes(rs_age+1, kexu))  # [g/day]
            
        else:
            rx = None
            seg_fluxes = None
            seg_sol_fluxes = None
        rx = comm.bcast(rx, root = 0)
        seg_fluxes = comm.bcast(seg_fluxes, root = 0)
        seg_sol_fluxes = comm.bcast(seg_sol_fluxes, root = 0)
        #print("rx",rx)
        
        
        """ 2. local soil models """
        wall_local = timeit.default_timer()
        if rank == 0:
            if len(net_flux) > 0:            
                proposed_outer_fluxes = rs.splitSoilFluxes(net_flux / dt, split_type) 
            else:
                proposed_outer_fluxes = np.full(len(r.cyls), 0.)
            proposed_outer_sol_fluxes = np.full(len(proposed_outer_fluxes), 0.)#rs.splitSoilFluxes(net_sol_flux / dt, split_type)
            # print("proposed_outer_fluxes",proposed_outer_fluxes,net_flux / dt)
            # print("proposed_outer_sol_fluxes",proposed_outer_sol_fluxes)#,net_sol_flux / dt)
            # if this fails, a segment is not mapped, i.e. out of soil domain
        else:
            proposed_outer_fluxes = None
            proposed_outer_sol_fluxes = None
        proposed_outer_fluxes = comm.bcast(proposed_outer_fluxes, root = 0)
        proposed_outer_sol_fluxes = comm.bcast(proposed_outer_sol_fluxes, root = 0)
        seg_rx = np.array([0.5 * (rx[seg.x] + rx[seg.y]) for seg in segs])
        
        soil_solute = r.getCC(len(cell_volumes),konz = False)#np.multiply(solute_conc, soil_water)
        
        if "dirichlet" in type:        
            r.solve(dt, seg_rx, proposed_outer_fluxes, seg_sol_fluxes, proposed_outer_sol_fluxes) #
        else:
            # solute fluxes are in mol/cm2 scv/d
            r.solve(dt, seg_fluxes, proposed_outer_fluxes, seg_sol_fluxes,proposed_outer_sol_fluxes) # 
        # water: left dirchlet, right neumann; solute: left and right Neumann<----
        # cyl = r.cyls[1]
        # write_file_array("pressureHead",np.array(cyl.getSolutionHead()).flatten())
        # write_file_array("coord", cyl.getDofCoordinates().flatten())
        # for i in range(r.numFluidComp):
            # write_file_array("solute_conc"+str(i+1), np.array(cyl.getSolution_(i+1)).flatten()* r.molarDensityWat ) 
        # for i in range(r.numFluidComp, r.numComp):
            # write_file_array("solute_conc"+str(i+1), np.array(cyl.getSolution_(i+1)).flatten()* r.bulkDensity_m3 /1e6 ) 

                    
        # r.checkMassOMoleBalance(doSolute = False)#I think that fails because the sucrose then flows out?
        wall_local = timeit.default_timer() - wall_local
        

        """ 3c. calculate mass net fluxes """
        solute_conc = r.getCC(len(cell_volumes))# np.array(s.getSolution_(1))/1e6 #mol/cm3
        # raise Exception
        try:
            assert min(solute_conc) >=0
        except:
            print("soil_sol_fluxes", soil_sol_fluxes, soil_fluxes)
            print("min(solute_conc)",min(solute_conc), min(solute_concOld))
            raise Exception
            
        new_soil_solute = r.getCC(len(solute_conc),konz = False)#np.multiply(solute_conc, soil_water)
        try:
            assert min(new_soil_solute) >=0
        except:
            print("min(new_soil_solute), min(soil_water)",min(new_soil_solute))
            raise Exception
        # net_sol_flux = np.multiply((new_soil_solute - soil_solute)/dt ,np.array(s.getCellVolumes()).flatten()) / 1e6 # change in water per cell [cm3]
        net_sol_flux = (new_soil_solute - soil_solute)/dt #* 1e6#/ (np.array(s.getCellVolumes()).flatten()/
        # print('net_sol_flux', net_sol_flux, (new_soil_solute - soil_solute).shape, (np.array(s.getCellVolumes())).shape) 
        #for k, root_sol_flux in soil_sol_fluxes.items():
        #    net_sol_flux[k] -= root_sol_flux * dt
            # print(k, dt, root_sol_flux,(new_soil_solute - soil_solute)[k], soil_sol_fluxes[k], net_sol_flux[k])
            # print((new_soil_solute - soil_solute )[k] -(soil_sol_fluxes[k]*dt + net_sol_flux[k]))
            # print((new_soil_solute - soil_solute )[k] -(root_sol_flux*dt + net_sol_flux[k]))
            # a=  Decimal((new_soil_solute - soil_solute )[k])
            # b = Decimal((root_sol_flux*dt))
            # temp = a - b
            # print(a -(b + temp))
            # net_sol_flux[k] = temp
            
        # print("soluteflows", (new_soil_solute - soil_solute)[cellIds],sum(soil_sol_fluxes.values())*dt,net_sol_flux[cellIds] )
        # print(cellIds,(new_soil_solute - soil_solute )[cellIds] ,soil_sol_fluxes,np.array(list(soil_sol_fluxes.values() ))*dt,( net_sol_flux[cellIds]))
        # print(((new_soil_solute - soil_solute)[cellIds] 
                    # -( np.array(list(soil_sol_fluxes.values()))*dt + net_sol_flux[cellIds]))
                    # /(( new_soil_solute - soil_solute)[cellIds])*100)
        #assert(((new_soil_solute - soil_solute)==0).all() )
        # assert((abs( ((new_soil_solute - soil_solute)[cellIds] 
                    # -( np.array(list(soil_sol_fluxes.values()))*dt + net_sol_flux[cellIds]))
                    # /(( new_soil_solute - soil_solute)[cellIds])*100)<0.1).all())
        #assert (abs((new_soil_solute - soil_solute )[cellIds] - np.array(list(soil_sol_fluxes.values()))*dt + net_sol_flux[cellIds]) < 1e-30).all()
        
        """ 3a. macroscopic soil model """
        wall_macro = timeit.default_timer()

        water_content = np.array(s.getWaterContent_())  # theta per cell [1]
        soil_water = np.multiply(water_content, cell_volumes)  # water per cell [cm3]
        solute_conc = np.array(s.getSolution_(1))/1e6 #mol/cm3
        # soil_solute = np.multiply(solute_conc, soil_water) #mol
        
        soil_fluxes = rs.sumSegFluxes(seg_fluxes)  # [cm3/day]  per soil cell
        # soil_sol_fluxes = rs.sumSegFluxes(seg_sol_fluxes)  # [mol/day]
        
        s.setSource(soil_fluxes.copy(), eq_idx = 0)  # [cm3/day], in modules/richards.py
        
        if ((len(net_sol_flux) >0) and max(abs(net_sol_flux)) != 0.):
            test_values = list(net_sol_flux.copy())
            test_keys = np.array([i for i in range(len(test_values))])
            res = {}
            for key in test_keys:
                for value in test_values:
                    if value > 0.:
                        res[key] = value
                    test_values.remove(value)
                    break
            # print("source", res, net_sol_flux[337], dt)
            # print(new_soil_solute[337], soil_solute[337] )
            s.setSource(res.copy(), eq_idx = 1)  # soil_sol_fluxes.copy(), eq_idx = 1)  # [mol/day], in modules/richards.py
            
        solute_concOld = solute_conc
        s.solve(dt)  # in modules/solverbase.py
        assert (s.getSolution_(1) == r.soilModel.getSolution_(1)).all()
        r.checkMassOMoleBalance(doWater = False)#I think that fails because the sucrose then flows out?
        
        test1 = np.array(s.getContent(1, True))
        test2 = np.array(r.getCC(len(solute_conc), konz = False))
        print(test1[np.where(test1!= 0.)], test2[np.where(test2!= 0.)])
        print(np.where(test1!= 0.), np.where(test2!= 0.))
        print((test1 == test2).all())
        raise Exception
        # print(" np.array(s.getCellVolumes()).flatten() ",  np.array(s.getCellVolumes()).flatten())
        wall_macro = timeit.default_timer() - wall_macro
        
        # if rank == 0:
        #     print("]", end = "")

        """ 3b. calculate net fluxes """
        wall_netfluxes = timeit.default_timer()
        water_content = np.array(s.getWaterContent_())
        soil_solute_content = np.array(s.getSolution_(1))
        # print(water_content)
        new_soil_water = np.multiply(water_content, cell_volumes)  # calculate net flux
        net_flux = new_soil_water - soil_water  # change in water per cell [cm3]
        #print('net_flux', net_flux) 
        for k, root_flux in soil_fluxes.items():
            net_flux[k] -= root_flux * dt #all the water lost at the outer boundary
        cellIds = np.fromiter(cell2seg.keys(), dtype=int)
        # print("waterflows",sum((new_soil_water - soil_water)[cellIds]),sum(soil_fluxes.values())*dt,sum(net_flux[cellIds]*dt),
                # sum(seg_fluxes))
        # print((new_soil_water - soil_water )[cellIds] ,np.array(list(soil_fluxes.values() )),( net_flux[cellIds]),dt)
        # print((new_soil_water - soil_water )[cellIds] -( np.array(list(soil_fluxes.values()))*dt + net_flux[cellIds]*dt))
        assert ((new_soil_water - soil_water )[cellIds] == np.array(list(soil_fluxes.values()))*dt + net_flux[cellIds]*dt).all()
        
        """ 3d. backup """
        # soil_water = new_soil_water
        # soil_solute = new_soil_solute

        wall_netfluxes = timeit.default_timer() - wall_netfluxes
        
        """ some checks """
        
        #does not use rsx for the leaves, but the pg (guard cells wat. pot.) values
        #[cm3/day]
        # proposed_inner_fluxes = rs.segFluxes(t, rx.copy(), rsx.copy(), approx = False, cells = False, soil_k = [])  # [cm3/day]
        # assert (np.array(proposed_inner_fluxes) == seg_fluxes).all()
        ###
        r.checkMassOMoleBalance()#I think that fails because the sucrose then flows out?
        
        try:
            assert min(water_content) >=0
            assert min(soil_solute_content) >=0
        except:
            print("min(new_soil_solute), min(soil_water)",min(water_content),min(soil_solute_content))
            raise Exception
            
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

            ana = pb.SegmentAnalyser(rs.rs.mappedSegments())  # VOLUME and SURFACE
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
                  "segments:", rs.rs.getNumberOfMappedSegments(), "root collar:", rx[0], "\n")
            print("time", rs_age + t, "rsx", np.min(rsx), np.max(rsx), "ccx", np.min(cc), np.max(cc))#, "trans",sum(np.array(r.Ev)))
            psi_x_.append(rx.copy())  # cm (per root node)
            psi_s_.append(rsx.copy())  # cm (per root segment)

        # if rank == 0:
        #     print("]")

    if rank == 0:
        print ("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s")

    return psi_x_, psi_s_, sink_, x_, y_, psi_s2_, vol_, surf_,  depth_, soil_c_, c_, repartition, c_All, net_sol_flux, net_flux #krs_,
