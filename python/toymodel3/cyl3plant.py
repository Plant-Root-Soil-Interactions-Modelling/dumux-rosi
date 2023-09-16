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


def simulate_const(s, rs, sim_time, dt, kexu, rs_age,repartition, type,Q_plant,
                    trans_maize,
                    plantType = "RS", r = [], wilting_point =-15000, 
                    outer_R_bc_sol=[], #mol
                    outer_R_bc_wat = []):#m3
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
    Q_Exud = Q_plant[0]; Q_mucil = Q_plant[1] #mol/day
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
    
    
    # repartitionOld = repartition
    # nsOld = sum(repartitionOld)
    # ns = len(rs.get_segments())
    # dcyl = int(np.floor(ns / max_rank))
    # repartition = np.array([dcyl for i in range(max_rank)])
    # repartition[max_rank -1] = ns - rank * dcyl
    # toAdd = repartition - repartitionOld
    
    # if toAdd[rank] > 0: # that could be put in main file (after RS growth)
        # if len(sx.shape)!=1:
            # raise Exception
        # if rank == 0:
            # r.update( sx, np.array([i for i in range(toAdd[rank])]) + nsOld)#, cc )
            # print ("Initialized rank {:g}/{:g} [{:g}-{:g}] in {:g} s".format(rank + 1, max_rank, nsOld, toAdd[rank] + nsOld, timeit.default_timer() - start_time))
        # else:
            # r.update( sx,  np.array([i for i in range(toAdd[rank-1],toAdd[rank])]) + nsOld)#, cc )
            # print ("Initialized rank {:g}/{:g} [{:g}-{:g}] in {:g} s".format(rank + 1, max_rank, toAdd[rank-1] + nsOld, toAdd[rank] + nsOld, timeit.default_timer() - start_time))
    # print("repartitionOld", repartitionOld, len(repartitionOld))
    # if sum ( repartitionOld)==0: 
        # print("first check")    
        # r.checkMassOMoleBalance()

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
            isPlant = True
            if isPlant:            
                # rs.minLoop = 1000
                # rs.maxLoop = 5000
                # rs.Qlight = weatherX["Qlight"]
                # rs.solve_photosynthesis(sim_time_ = rs_age + t, 
                            # sxx_=rsx, #will not use the output of the air bc. so it s just a stand-in for now
                            # cells_ = False,
                            # ea_ = weatherX["ea"],#not used
                            # es_=weatherX["es"],#not used
                            # verbose_ = False, doLog_ = False,
                            # TairC_= weatherX["TairC"],#not used
                            # outputDir_= "./results/rhizoplantExud")
                rx = np.array(rs.psiXyl)
                seg_fluxes = np.array(rs.outputFlux)# [cm3/day]
                errLeuning = sum(seg_fluxes)            
            else:
                rx = rs.solve(rs_age, -2, 0., rsx, cells = False, 
                                    wilting_point = wilting_point, soil_k = [])
                #trans_maize(rs_age + 0., dt)*10
                rx = np.array(rx)
                #seg_fluxes =np.full(len( np.array(rs.segFluxes(rs_age + t, rx, rsx, False, False, []))),0.)# [cm3/day]
                seg_fluxes = np.array(rs.segFluxes(rs_age + t, rx, rsx, False, False, []))# [cm3/day]
                # print("rx, rsx", rx, rsx)
            
            seg_sol_fluxes = Q_Exud # mol/day for segments
            seg_mucil_fluxes = Q_mucil
            #r.exudate_fluxes(rs_age+1, kexu))  # [g/day]
            
        else:
            rx = None
            seg_fluxes = None
            seg_sol_fluxes = None
            seg_mucil_fluxes = None
        rx = comm.bcast(rx, root = 0)
        seg_fluxes = comm.bcast(seg_fluxes, root = 0)
        seg_sol_fluxes = comm.bcast(seg_sol_fluxes, root = 0)
        seg_mucil_fluxes = comm.bcast(seg_mucil_fluxes, root = 0)
        #print("rx",rx)
        
        """ 2. local soil models """
        wall_local = timeit.default_timer()
        if rank == 0:
            if len(outer_R_bc_wat) > 0:            
                proposed_outer_fluxes = rs.splitSoilFluxes(outer_R_bc_wat / dt, split_type) #cm3/day
            else:
                proposed_outer_fluxes = np.full(len(r.cyls), 0.)
            
            if len(outer_R_bc_sol[0]) > 0:            
                proposed_outer_sol_fluxes = rs.splitSoilFluxes(outer_R_bc_sol[0] / dt, split_type)#mol/day
                proposed_outer_mucil_fluxes = rs.splitSoilFluxes(outer_R_bc_sol[1] / dt, split_type)
            else:
                proposed_outer_sol_fluxes = np.full(len(r.cyls), 0.)
                proposed_outer_mucil_fluxes = np.full(len(r.cyls), 0.)
            # if this fails, a segment is not mapped, i.e. out of soil domain
        else:
            proposed_outer_fluxes = None
            proposed_outer_sol_fluxes = None
        proposed_outer_fluxes = comm.bcast(proposed_outer_fluxes, root = 0)
        proposed_outer_sol_fluxes = comm.bcast(proposed_outer_sol_fluxes, root = 0)
        proposed_outer_mucil_fluxes = comm.bcast(proposed_outer_mucil_fluxes, root = 0)
        seg_rx = np.array([0.5 * (rx[seg.x] + rx[seg.y]) for seg in segs])
        
        #mol
        soil_solute = np.array([np.array(r.getCC(len(cell_volumes), idComp = idc + 1, konz = False)) for idc in range(r.numFluidComp)])
        
    
        if False:       
            i = 24
            cyl = r.cyls[i]
            l = r.seg_length[i]            
            watVol = sum(cyl.getWaterVolumesCyl(l))
            a_in = r.radii[i]
            a_out = r.outer_radii[i]
            
            print("cyl id",i,"water volume",watVol, "radius and length",a_in, a_out,l )
            times = [0., 1.0]  # days
            for jj, dt in enumerate(np.diff(times)):
            
                watVolBU = watVol
                QflowOut = proposed_outer_fluxes[i] #qflowOut * (2 * np.pi * a_out * l)
                
                qflowOut = QflowOut/(2 * np.pi * a_out * l)# 0.26 * (i+1)/4
                #-0.009046391194856907
                QflowIn = seg_fluxes[i]#qflow * (2 * np.pi * a_in * l)
                qflow = QflowIn/(2 * np.pi * a_in * l)#-0.009046389614038853 #- 0.26 * (i+1)
                cyl.setInnerFluxCyl(qflow)  # [cm3/day] -> [cm /day]
                cyl.setOuterFluxCyl(qflowOut)  # [cm3/day] -> [cm /day] 
             
                cyl.solve(dt, maxDt = 2500/(24*3600))                
                watVol = sum(cyl.getWaterVolumesCyl(l))
                print("water volume after",watVol,"dt",dt)
                print(  "QflowIn",QflowIn, "Qflow_dt", QflowIn * dt, "QflowOut", QflowOut)
                print("diff",watVolBU + (QflowIn + QflowOut) * dt)
                
        
                print("water conservation ",i, watVol , watVolBU )
                print("QflowIn", QflowIn ,"QflowOut", QflowOut,"dt", dt)
                print("qflowout",QflowOut , "qflowin",qflow)
                print( "l",l,"a_in", r.radii[i], a_in ,"a_out",r.outer_radii[i], a_out )
                print("diff",( watVolBU + (QflowIn + QflowOut) * dt))
                
            raise Exception
        
        if "dirichlet" in type:        
            r.solve(dt, seg_rx, proposed_outer_fluxes, seg_sol_fluxes, proposed_outer_sol_fluxes, seg_mucil_fluxes, proposed_outer_mucil_fluxes) #
        else:
            # solute fluxes are in mol/cm2 scv/d
            #print(seg_sol_fluxes)
            #print(proposed_outer_sol_fluxes)
            r.solve(dt, seg_fluxes, proposed_outer_fluxes, seg_sol_fluxes,proposed_outer_sol_fluxes, 
                    seg_mucil_fluxes, proposed_outer_mucil_fluxes) # cm3/day or mol/day
        
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
        # solute_conc = [r.getCC(len(cell_volumes), idComp = idc + 1, konz = False) for idc in range(r.numFluidComp)]# np.array(s.getSolution_(1))/1e6 #mol/cm3
        # # raise Exception
        # try:
            # assert min(solute_conc[0]) >=0
            # assert min(solute_conc[1]) >=0
        # except:
            # print("soil_sol_fluxes", soil_sol_fluxes, soil_fluxes)
            # print("min(solute_conc)",min(solute_conc[0]),min(solute_conc[1]))
            # raise Exception
        
        # mol
        new_soil_solute =np.array( [np.array(r.getCC(len(cell_volumes), idComp = idc + 1, konz = False)) for idc in range(r.numFluidComp)])
        try:
            assert min(new_soil_solute[0]) >=0
            assert min(new_soil_solute[1]) >=0
        except:
            print("min(new_soil_solute), min(soil_water)",min(new_soil_solute[0]),min(new_soil_solute[1]))
            raise Exception
        
        if(len(outer_R_bc_sol[0]) == 0):
            outer_R_bc_sol = np.full(new_soil_solute.shape, 0.)
        
        #mol/day
        soil_source_sol = np.array([np.array((new_soil_solute[i] - soil_solute[i] - outer_R_bc_sol[i])/dt )for i in range(2)])
        try:
            assert len(soil_source_sol[0].shape) == 1
        except:
            print(soil_source_sol[0].shape)
            print(new_soil_solute[0].shape)
            print(soil_solute[0].shape)
            print(outer_R_bc_sol[0].shape)
            print(soil_solute[0])
            raise Exception
        
        soil_fluxes = rs.sumSegFluxes(seg_fluxes)  # [cm3/day]  per soil cell
        # soil_sol_fluxes = rs.sumSegFluxes(seg_sol_fluxes)  # [mol/day]
        
        """ 3.d  some checks """
        # print("new_soil_solute[0][437] - soil_solute[0][437]  - outer_R_bc_sol[0][437]",
                # new_soil_solute[0][437] , soil_solute[0][437]  , outer_R_bc_sol[0][437] )
        # print("new_soil_solute[1][437] - soil_solute[1][437]  - outer_R_bc_sol[1][437]",
                # new_soil_solute[1][437] , soil_solute[1][437]  , outer_R_bc_sol[1][437] )
        # print((r.soilModel.getWaterVolumes()/1e6)[437] )
        r.checkMassOMoleBalance2(soil_fluxes, soil_source_sol, dt,seg_fluxes =seg_fluxes, doSolid = False)
        r.setSoilData(soil_fluxes, soil_source_sol, dt)
        """ 2.0  global soil models """
        wall_macro = timeit.default_timer()

        water_content = np.array(s.getWaterContent_())  # theta per cell [1]
        soil_water = np.multiply(water_content, cell_volumes)  # water per cell [cm3]
        
        soil_solute_content = np.array([np.array(s.getContent(i+1, isDissolved = True)) for i in range(2)])
        #mol
        
        # solute_conc = np.array(s.getSolution_(1))/1e6 #mol/cm3
        # soil_solute = np.multiply(solute_conc, soil_water) #mol
        
        
        s.setSource(soil_fluxes.copy(), eq_idx = 0)  # [cm3/day], in modules/richards.py
        
        for idComp in range(2):#mol/day
            if ((len(soil_source_sol[idComp]) >0) and max(abs(soil_source_sol[idComp])) != 0.):
                test_values = list(soil_source_sol[idComp].copy())
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
                s.setSource(res.copy(), eq_idx = idComp+1)  # soil_sol_fluxes.copy(), eq_idx = 1)  # [mol/day], in modules/richards.py
                
        
        s.solve(dt)  # in modules/solverbase.py
        assert (s.getSolution_(0) == r.soilModel.getSolution_(0)).all()
        # print(" np.array(s.getCellVolumes()).flatten() ",  np.array(s.getCellVolumes()).flatten())
        wall_macro = timeit.default_timer() - wall_macro
        
        # if rank == 0:
        #     print("]", end = "")

        """ 3b. calculate net fluxes """
        wall_netfluxes = timeit.default_timer()
        water_content = np.array(s.getWaterContent_())
        # print(water_content)
        new_soil_water = np.multiply(water_content, cell_volumes)  # calculate net flux
        outer_R_bc_wat = new_soil_water - soil_water  # change in water per cell [cm3]
        #print('net_flux', net_flux) 
        for k, root_flux in soil_fluxes.items():
            outer_R_bc_wat[k] -= root_flux * dt #all the water lost at the outer boundary
        cellIds = np.fromiter(cell2seg.keys(), dtype=int)
        cellIds =  np.array([x for x in cellIds if x >= 0])#take out air segments
        # print("waterflows",sum((new_soil_water - soil_water)[cellIds]),sum(soil_fluxes.values())*dt,sum(net_flux[cellIds]*dt),
                # sum(seg_fluxes))
        # print((new_soil_water - soil_water )[cellIds] -( np.array(list(soil_fluxes.values()))*dt + net_flux[cellIds]*dt))
        try:
            assert (abs(((new_soil_water - soil_water )[cellIds] - (np.array(list(soil_fluxes.values()))*dt + outer_R_bc_wat[cellIds]))/((new_soil_water - soil_water )[cellIds])) < 0.1).all()
        except:
            print(new_soil_water[cellIds] - (soil_water[cellIds] + np.array(list(soil_fluxes.values()))* dt + outer_R_bc_wat[cellIds]))
            raise Exception
        
        soil_solute_content_new = np.array([np.array(s.getContent(i+1, isDissolved = True)) for i in range(2)])
        outer_R_bc_sol = np.array([soil_solute_content_new[i] - soil_solute_content[i] - soil_source_sol[i]*dt for i in range(2)])
        try:
            assert len(outer_R_bc_sol[0].shape) == 1
        except:
            print(outer_R_bc_sol[0].shape)
            print(soil_solute_content_new[0].shape)
            print(soil_solute_content[0].shape)
            print(soil_source_sol[0].shape)
            raise Exception
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
        # r.checkMassOMoleBalance(doSolute = False)#I think that fails because the sucrose then flows out?
        
        try:
            assert min(water_content) >=0
            assert min(soil_solute_content_new[0]) >=0
            assert min(soil_solute_content_new[1]) >=0
        except:
            print("min(new_soil_solute), min(soil_water)",min(water_content),min(soil_solute_content[0]),min(soil_solute_content[1]))
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

    return psi_x_, psi_s_, sink_, x_, y_, psi_s2_, vol_, surf_,  depth_, soil_c_, c_, repartition, c_All, outer_R_bc_sol, outer_R_bc_wat 

