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
    split_type = 0  # type 0 == volume, type 1 == surface, type 2 == length
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
    
    
    segs = rs.rs.segments
    Nt = len(rs.rs.nodes)
    # print(len(segs), Nt)
    assert len(segs) == (Nt -1)
    seg2cell = rs.rs.seg2cell
    cell2seg = rs.rs.cell2seg
    cellIds = np.fromiter(cell2seg.keys(), dtype=int)
    cellIds =  np.array([x for x in cellIds if x >= 0])#take out air segments
    
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
            
            seg_sol_fluxes = Q_Exud /dt# mol/day for segments
            seg_mucil_fluxes = Q_mucil/dt
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
                waterContent = np.array([sum( r.cyls[sid].getWaterVolumesCyl(r.seg_length[sid]) ) for sid in range(len(r.cyls))]) 
                proposed_outer_fluxes = r.splitSoilVals(outer_R_bc_wat / dt, waterContent) #cm3/day
            else:
                proposed_outer_fluxes = np.full(len(r.cyls), 0.)
            
            if len(outer_R_bc_sol[0]) > 0:  
                comp1content = np.array([sum( r.cyls[sid].getContentCyl(1, True, r.seg_length[sid]) ) for sid in range(len(r.cyls))])
                comp2content = np.array([sum( r.cyls[sid].getContentCyl(2, True, r.seg_length[sid]) ) for sid in range(len(r.cyls))])
                proposed_outer_sol_fluxes = r.splitSoilVals(outer_R_bc_sol[0] / dt, comp1content)#mol/day
                proposed_outer_mucil_fluxes = r.splitSoilVals(outer_R_bc_sol[1] / dt, comp2content)
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
        soil_solute = np.array([np.array(r.getC_rhizo(len(cell_volumes), idComp = idc + 1, konz = False)) for idc in range(r.numComp)])
        
        if len(outer_R_bc_sol[0]) > 0:  #making sure that the mass balance is ok at the voxel level
            seg_sol_fluxes_voxel = rs.sumSegFluxes(seg_sol_fluxes*dt) 
            seg_mucil_fluxes_voxel = rs.sumSegFluxes(seg_mucil_fluxes*dt)
            for cid in cellIds: 
                try:
                    assert soil_solute[0][cid] + outer_R_bc_sol[0][cid] + seg_sol_fluxes_voxel[cid] >= 0
                    assert soil_solute[1][cid] + outer_R_bc_sol[1][cid] + seg_mucil_fluxes_voxel[cid] >= 0
                except:
                    print("soil_solute[0]",soil_solute[0][cid] ," outer_R_bc_sol[0] ", outer_R_bc_sol[0][cid] , "seg_sol_fluxes_voxel",seg_sol_fluxes_voxel[cid])
                    print("soil_solute[1]",soil_solute[1][cid] , "outer_R_bc_sol[1]", outer_R_bc_sol[1][cid] , "seg_mucil_fluxes_voxel",seg_mucil_fluxes_voxel[cid])
                    print(cid)
                    raise Exception
        assert len(seg_fluxes) == len(organTypes)
        assert len(proposed_outer_fluxes) == len(organTypes)
        assert len(seg_sol_fluxes) == len(organTypes)
        assert len(proposed_outer_sol_fluxes) == len(organTypes)
        assert len(seg_mucil_fluxes) == len(organTypes)
        assert len(proposed_outer_mucil_fluxes) == len(organTypes)
        
        if len(outer_R_bc_sol[0]) > 0:  #making sure that the mass balance is ok at the rhizosphere level
            for sid in seg2cell: 
                if ((seg2cell[sid] >= 0) and (organTypes[sid] == 2)): # root segment in the soil
                    try:
                        assert sum( r.cyls[sid].getContentCyl(1, True, r.seg_length[sid]) )+ proposed_outer_sol_fluxes[sid]*dt+ seg_sol_fluxes[sid] *dt>= 0
                        assert sum(r.cyls[sid].getContentCyl(2, True, r.seg_length[sid]))  + proposed_outer_mucil_fluxes[sid] *dt+ seg_mucil_fluxes[sid]*dt>= 0
                    except:
                        cid = seg2cell[sid]
                        sids = cell2seg[cid]
                        totCont = 0
                        totouterBC = 0
                        totinnerBC = 0
                        for ssid in sids:     
                            totCont += sum(r.cyls[ssid].getContentCyl(1, True, r.seg_length[ssid]))
                            totouterBC +=  proposed_outer_sol_fluxes[ssid] 
                            totinnerBC += seg_sol_fluxes[ssid]
                            print("seg id ",ssid,organTypes[sid], r.seg_length[ssid])        
                            print("soil_solute[0]",sum(r.cyls[ssid].getContentCyl(1, True, r.seg_length[ssid])) ,
                                    " outer_R_bc_sol[0] ", proposed_outer_sol_fluxes[ssid]  ,
                                       "innerBC_sol",seg_sol_fluxes[ssid])
                            print("soil_solute[1]",sum(r.cyls[ssid].getContentCyl(2, True, r.seg_length[ssid]) ),
                                    "outer_R_bc_sol[1]", proposed_outer_mucil_fluxes[ssid] , 
                                    "innerBC_mucil",seg_mucil_fluxes[ssid])
                        
                        print("the one with the issue")                        
                        print("soil_solute[0]",sum(r.cyls[sid].getContentCyl(1, True, r.seg_length[sid])) ,
                                " outer_R_bc_sol[0] ", proposed_outer_sol_fluxes[sid]  ,
                                   "innerBC_sol",seg_sol_fluxes[sid])
                        print("soil_solute[1]",sum(r.cyls[sid].getContentCyl(2, True, r.seg_length[sid]) ),
                                "outer_R_bc_sol[1]", proposed_outer_mucil_fluxes[sid] , 
                                "innerBC_mucil",seg_mucil_fluxes[sid])
                        print(sid, r.seg_length[sid])
                        print(totCont,totouterBC*dt, totinnerBC*dt)
                        raise Exception
                    
        
        print("solve 1d soil")
        if "dirichlet" in type:        
            r.solve(dt, seg_rx, proposed_outer_fluxes, seg_sol_fluxes, proposed_outer_sol_fluxes, seg_mucil_fluxes, proposed_outer_mucil_fluxes) #
        else:
            # solute fluxes are in mol/cm2 scv/d
            r.solve(dt, seg_fluxes, proposed_outer_fluxes, seg_sol_fluxes,proposed_outer_sol_fluxes, 
                    seg_mucil_fluxes, proposed_outer_mucil_fluxes) # cm3/day or mol/day
        print("done")
        # cyl = r.cyls[1]
        # write_file_array("pressureHead",np.array(cyl.getSolutionHead()).flatten())
        # write_file_array("coord", cyl.getDofCoordinates().flatten())
        # for i in range(r.numFluidComp):
            # write_file_array("solute_conc"+str(i+1), np.array(cyl.getSolution_(i+1)).flatten()* r.molarDensityWat ) 
        # for i in range(r.numFluidComp, r.numComp):
            # write_file_array("solute_conc"+str(i+1), np.array(cyl.getSolution_(i+1)).flatten()* r.bulkDensity_m3 /1e6 ) 

                    
        wall_local = timeit.default_timer() - wall_local
        

        """ 3c. calculate mass net fluxes """
        
        
        # mol
        new_soil_solute =np.array( [np.array(r.getC_rhizo(len(cell_volumes), idComp = idc + 1, konz = False)) for idc in range(r.numComp)])
        try:
            assert min(new_soil_solute.flatten()) >=0
        except:
            print("min(new_soil_solute)",min(new_soil_solute.flatten()),[min(nss) for nss in new_soil_solute])
            raise Exception
        
        if(len(outer_R_bc_sol[0]) == 0):
            outer_R_bc_sol = np.full(new_soil_solute.shape, 0.)
        
        #mol/day
        soil_source_sol = np.full(new_soil_solute.shape,0. )
        for nc in range(r.numComp):
            soil_source_sol[nc][cellIds] = np.array(new_soil_solute[nc][cellIds] - soil_solute[nc][cellIds] - outer_R_bc_sol[nc][cellIds])/dt
            #assert:
            
        try:
            assert soil_source_sol.shape == (r.numComp, len(cell_volumes))
        except:
            print(soil_source_sol[0].shape)
            print(new_soil_solute[0].shape)
            print(soil_solute[0].shape)
            print(outer_R_bc_sol[0].shape)
            print(soil_solute[0])
            raise Exception
            
        soil_fluxes = rs.sumSegFluxes(seg_fluxes)  # [cm3/day]  per soil cell
        
        """ 3.d  some checks """
        r.checkMassOMoleBalance2(soil_fluxes, soil_source_sol, dt,seg_fluxes =seg_fluxes, doSolid = True)
        r.setSoilData(soil_fluxes, soil_source_sol, dt)
        """ 2.0  global soil models """
        wall_macro = timeit.default_timer()

        water_content = np.array(s.getWaterContent_())  # theta per cell [1]
        soil_water = np.multiply(water_content, cell_volumes)  # water per cell [cm3]
        
        soil_solute_content = np.array([np.array(s.getContent(i+1, isDissolved = (i < r.numFluidComp))) for i in range(r.numComp)])
        #mol
        
        
        s.setSource(soil_fluxes.copy(), eq_idx = 0)  # [cm3/day], in modules/richards.py
        
        for idComp in range(s.numComp):#mol/day
            if ((len(soil_source_sol[idComp]) >0) and max(abs(soil_source_sol[idComp])) != 0.):
                try:
                    assert min(soil_solute_content[idComp] + soil_source_sol[idComp]*dt) >= 0.
                except:
                    print(idComp, min(soil_solute_content[idComp] + soil_source_sol[idComp]*dt), dt)
                    outputPredicted = np.array(soil_solute_content[idComp] + soil_source_sol[idComp]*dt)
                    
                    print(soil_solute_content[idComp] , soil_source_sol[idComp]*dt)
                    print(np.where(outputPredicted == min(outputPredicted)))
                    raise Exception
                test_values = list(soil_source_sol[idComp].copy())
                test_keys = np.array([i for i in range(len(test_values))])
                res = {}
                for key in test_keys:
                    for value in test_values:
                        res[key] = value
                        test_values.remove(value)
                        break
                        
                try:
                    assert len(res) == len(test_keys)
                except:
                    print(res,len(res),len(test_keys))
                    print("soil_source_sol[idComp]",soil_source_sol[idComp])
                    raise Exception
                    
                s.setSource(res.copy(), eq_idx = idComp+1)  # [mol/day], in modules/richards.py
                
        
        buTotCBefore =sum(s.getTotCContent())
        print("solve 3d soil")
        s.solve(dt)  # in modules/solverbase.py
        print("done")
        buTotCAfter = sum(s.getTotCContent())    
        
        try:
            # maybe 0.1% of error is too large
            assert abs((buTotCAfter - ( buTotCBefore + sum(Q_Exud) + sum(Q_mucil)))/buTotCAfter*100) < 0.1
            # assert bellow can fail when Q_exud and Q_mucil low
            #assert abs((buTotCAfter - ( buTotCBefore + sum(Q_Exud) + sum(Q_mucil) ))) < 0.1*( sum(Q_Exud) + sum(Q_mucil))
        except:
            print("buTotCAfter ,  buTotCBefore", buTotCAfter ,  buTotCBefore)
            print( "sum(Q_Exud) , sum(Q_mucil)", sum(Q_Exud) , sum(Q_mucil))
            print(abs((buTotCAfter - ( buTotCBefore + sum(Q_Exud) + sum(Q_mucil) ))))
            raise Exception

        assert (s.getSolution_(0) == r.soilModel.getSolution_(0)).all()
        wall_macro = timeit.default_timer() - wall_macro
        

        """ 3b. calculate net fluxes """
        wall_netfluxes = timeit.default_timer()
        water_content = np.array(s.getWaterContent_())
        # print(water_content)
        new_soil_water = np.multiply(water_content, cell_volumes)  # calculate net flux
        try:
            # use soil_fluxes and not seg_fluxes. seg_fluxes includes air segments. sum(seg_fluxes) ~ 0.
            # maybe 0.1% of error is too large
            assert abs(  (sum(new_soil_water) - (sum(soil_water) + sum(soil_fluxes.values())*dt) ) /sum(new_soil_water) )*100 < 0.1
        except:
            print(new_soil_water , soil_water , seg_fluxes*dt)
            print(sum(new_soil_water) , sum(soil_water) , sum(soil_fluxes.values())*dt)
            print(abs(sum(new_soil_water) - (sum(soil_water) + sum(soil_fluxes.values())*dt)) )
            print( abs(  (sum(new_soil_water) - (sum(soil_water) + sum(soil_fluxes.values())*dt) ) /sum(new_soil_water) )*100  )
            print(new_soil_water)
            raise Exception
        
        
        outer_R_bc_wat = new_soil_water - soil_water  # change in water per cell [cm3]
        
        for k, root_flux in soil_fluxes.items():
            outer_R_bc_wat[k] -= root_flux * dt #all the water lost at the outer boundary
        
        try:
            assert (abs(((new_soil_water - soil_water )[cellIds] - (np.array(list(soil_fluxes.values()))*dt + outer_R_bc_wat[cellIds]))/((new_soil_water - soil_water )[cellIds])) < 0.1).all()
        except:
            print(new_soil_water[cellIds] - (soil_water[cellIds] + np.array(list(soil_fluxes.values()))* dt + outer_R_bc_wat[cellIds]))
            raise Exception
        
        soil_solute_content_new = np.array([np.array(s.getContent(i+1, isDissolved = (i < r.numFluidComp))) for i in range(r.numComp)])# mol
        outer_R_bc_sol = np.array([soil_solute_content_new[i] - soil_solute_content[i] - soil_source_sol[i]*dt for i in range(r.numComp)])# mol
        try:
            assert outer_R_bc_sol.shape == (r.numComp, len(cell_volumes))
        except:
            print(outer_R_bc_sol)
            print(outer_R_bc_sol.shape)
            print(r.numComp, len(cell_volumes))
            print(outer_R_bc_sol[0].shape)
            print(soil_solute_content_new[0].shape)
            print(soil_solute_content[0].shape)
            print(soil_source_sol[0].shape)
            raise Exception
        for nc in range(r.numComp, r.numFluidComp):
            # all changes in cellIds for 3D soil is cause by biochemical reactions computed in 1D models.
            # thus, for elements which do not flow (range(r.numComp, r.numFluidComp)), there are no changes
            # except those prescribed by the 1d model.
            assert (outer_R_bc_sol[nc][cellIds] == 0.).all()
            # outer_R_bc_sol[nc][cellIds] == flow in 3D model
            # outer_R_bc_sol[nc][not cellIds] == flow and reaction in 3D model
        
        """ 3d. backup """
        wall_netfluxes = timeit.default_timer() - wall_netfluxes
        
        """ some checks """
        
        try:
            assert min(water_content) >=0
            assert min(soil_solute_content_new.flatten()) >=0
        except:
            print("min(new_soil_solute), min(soil_water)",min(water_content),[min(ssc) for ssc in soil_solute_content])
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


    if rank == 0:
        print ("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s")

    return psi_x_, psi_s_, sink_, x_, y_, psi_s2_, vol_, surf_,  depth_, soil_c_, c_, repartition, c_All, outer_R_bc_sol, outer_R_bc_wat 

