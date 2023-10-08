"""
functions for the cylindrical coupling approach
"""
import numpy as np
import timeit
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank(); max_rank = comm.Get_size()
import sys
from xylem_flux import *
from air_modelsPlant import AirSegment
import vtk_plot as vtk
import plantbox as pb
import evapotranspiration as evap
import timeit
import visualisation.vtk_plot as vp
from scenario_setup import weather
from decimal import *

from scenario_setup import write_file_array


def simulate_const(s, rs, sim_time, dt, rs_age, demoType,Q_plant,
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
    subtypes = np.asarray(rs.rs.subTypes, int)
    organTypes = np.asarray(rs.rs.organTypes, int)
    
    airSegsId = np.array(list(set(np.concatenate((r.cell2seg[-1],np.where(np.array(r.organTypes) != 2)[0])) )))#aboveground
    #rootAirSegsId  = np.array([(isinstance(cc, AirSegment) and (organTypes[i] == 2)) for i, cc in enumerate(r.cyls)])
    rhizoSegsId = np.array([i for i in range(len(organTypes)) if i not in airSegsId])
    Q_Exud = Q_plant[0]; Q_mucil = Q_plant[1] #mol/day
    try:
        assert (Q_Exud[airSegsId] == 0).all()
        assert (Q_mucil[airSegsId] == 0).all()
    except:
        print(Q_Exud[airSegsId] )
        print(Q_mucil[airSegsId])
        #print(Q_Exud[rootAirSegsId] )
        #print(Q_mucil[rootAirSegsId])
        print(organTypes[airSegsId])
        raise Exception
        
    #wilting_point = -15000  # cm
    skip = 10  # 3 * 6  # for output and results, skip iteration
    
    
        
    assert len(rs.rs.segments) == (len(rs.rs.nodes) -1)
    seg2cell = rs.rs.seg2cell
    cell2seg = rs.rs.cell2seg
    cellIds = np.fromiter(cell2seg.keys(), dtype=int)
    cellIds =  np.array([x for x in cellIds if x >= 0])#take out air segments
    
    #r = rs.rs  # rename (XylemFluxPython)

    psi_x_, psi_s_, sink_ , x_, y_, psi_s2_, soil_c_, c_ , c_All, c_All1= [], [], [], [], [], [], [], [] ,[], []
    vol_ = [[], [], [], [], [], []]
    surf_ = [[], [], [], [], [], []]
    krs_ = []
    depth_ = []
    # for post processing

    
    cell_volumes = comm.bcast(s.getCellVolumes() , root = 0) #cm3
    #net_flux = np.zeros(cell_volumes.shape)
    #net_sol_flux = np.zeros(cell_volumes.shape)

    N = int(np.ceil(sim_time / dt))  # number of iterations
    
    """ simualtion loop """
    for i in range(N):

        t = rs_age + i * dt  # current simulation time
        
        weatherX = weather(t)

        """ 1. xylem model """
        # a = np.array(rs.rs.radii)
        # l = np.array(rs.rs.segLength())
    
        ###
        
        comm.barrier()
        

            
        rx = np.array(rs.psiXyl)
        seg_fluxes = np.array(rs.outputFlux)# [cm3/day] 
        seg_sol_fluxes = Q_Exud /dt# mol/day for segments
        seg_mucil_fluxes = Q_mucil/dt
        
        assert min(seg_sol_fluxes) >= 0
        assert min(seg_mucil_fluxes) >= 0
        
        """ 2. local soil models """
        waterContent = r.getWaterVolumesCyl(doSum = False, reOrder = True)#np.array([sum( r.cyls[sid].getWaterVolumesCyl(r.seg_length[sid]) ) for sid in range(len(r.cyls))]) 
        comp1content = r.getContentCyl(idComp=1, doSum = False, reOrder = True)#np.array([sum( r.cyls[sid].getContentCyl(1, True, r.seg_length[sid]) ) for sid in range(len(r.cyls))])
        comp2content = r.getContentCyl(idComp=2, doSum = False, reOrder = True)#np.array([sum( r.cyls[sid].getContentCyl(2, True, r.seg_length[sid]) ) for sid in range(len(r.cyls))])
        
        try:
            assert waterContent.shape == (len(organTypes), )
        except:
            print(waterContent,waterContent.shape)
            raise Exception
        
        if rank == 0:
            if len(outer_R_bc_wat) > 0:           
                proposed_outer_fluxes = r.splitSoilVals(outer_R_bc_wat / dt, waterContent) #cm3/day
            else:
                proposed_outer_fluxes = np.full(len(organTypes), 0.)
            
            if len(outer_R_bc_sol[0]) > 0:  
                proposed_outer_sol_fluxes = r.splitSoilVals(outer_R_bc_sol[0] / dt, comp1content)#mol/day
                proposed_outer_mucil_fluxes = r.splitSoilVals(outer_R_bc_sol[1] / dt, comp2content)
            else:
                proposed_outer_sol_fluxes = np.full(len(organTypes), 0.)
                proposed_outer_mucil_fluxes = np.full(len(organTypes), 0.)
            # if this fails, a segment is not mapped, i.e. out of soil domain
        else:
            proposed_outer_fluxes = None
            proposed_outer_sol_fluxes = None
            proposed_outer_mucil_fluxes = None
        proposed_outer_fluxes = comm.bcast(proposed_outer_fluxes, root = 0)
        proposed_outer_sol_fluxes = comm.bcast(proposed_outer_sol_fluxes, root = 0)
        proposed_outer_mucil_fluxes = comm.bcast(proposed_outer_mucil_fluxes, root = 0)
        
        
        #mol
        soil_solute = np.array([np.array(r.getC_rhizo(len(cell_volumes), idComp = idc + 1, konz = False)) for idc in range(r.numComp)])
        if (rank == 0) and ( len(outer_R_bc_sol[0]) > 0):  #making sure that the mass balance is ok at the voxel level
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
        
        if False:# see how to handle that in case of MPI. len(outer_R_bc_sol[0]) > 0:  #making sure that the mass balance is ok at the rhizosphere level
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
                        
                        print("the one with the issue", sid, cid)                        
                        print("soil_solute[0]",sum(r.cyls[sid].getContentCyl(1, True, r.seg_length[sid])) ,
                                " outer_R_bc_sol[0] ", proposed_outer_sol_fluxes[sid]  ,
                                   "innerBC_sol",seg_sol_fluxes[sid])
                        print("soil_solute[1]",sum(r.cyls[sid].getContentCyl(2, True, r.seg_length[sid]) ),
                                "outer_R_bc_sol[1]", proposed_outer_mucil_fluxes[sid] , 
                                "innerBC_mucil",seg_mucil_fluxes[sid])
                        print(sid, r.seg_length[sid])
                        print(totCont,totouterBC,dt, "totinnerBC",totinnerBC,dt, totCont+totouterBC*dt+ totinnerBC*dt)
                        
                        r.splitSoilVals(outer_R_bc_sol[0] / dt, comp1content, troubleShootId = cid)
                        raise Exception
                    
        #rhizoWBefore__ = comm.bcast(r._flat0(comm.gather(np.array([cc.getWaterVolumesCyl(r.seg_length_[i]) for i, cc in enumerate(r.cyls)], dtype=object), root = 0)),root = 0)
        rhizoWBefore_ = r.getWaterVolumesCyl(doSum = False, reOrder = True)#np.array([ sum(cylW) for cylW in rhizoWBefore__])
        rhizoWBefore = sum(rhizoWBefore_) 
        #rhizoTotCBefore__ = comm.bcast(r._flat0(comm.gather(np.array([ r.getTotCContent(i, cc,r.seg_length_[i]) for i, cc in enumerate(r.cyls)], dtype=object), root = 0)),root = 0) 
        rhizoTotCBefore_ = r.getTotCContentAll(doSum = False, reOrder = True)#np.array([ sum(cylC) for cylC in rhizoTotCBefore__]) 
        rhizoTotCBefore = sum(rhizoTotCBefore_) 
        
        start_time_rhizo = timeit.default_timer()
        print("solve 1d soil", rank)
        if "dirichlet" in demoType:        
            seg_rx = np.array([0.5 * (rx[seg.x] + rx[seg.y]) for seg in rs.rs.segments])
            r.solve(dt, seg_rx, proposed_outer_fluxes, seg_sol_fluxes, proposed_outer_sol_fluxes, seg_mucil_fluxes, proposed_outer_mucil_fluxes) #
        else:
            # solute fluxes are in mol/cm2 scv/d
            r.solve(dt, seg_fluxes, proposed_outer_fluxes, seg_sol_fluxes,proposed_outer_sol_fluxes, 
                    seg_mucil_fluxes, proposed_outer_mucil_fluxes) # cm3/day or mol/day
        print("solve 1d soil", rank)
        
        rs.time_rhizo_i = (timeit.default_timer() - start_time_rhizo)
    
        
        #rhizoWAfter__ = comm.bcast(r._flat0(comm.gather(np.array([cc.getWaterVolumesCyl(r.seg_length_[i]) for i, cc in enumerate(r.cyls)], dtype=object), root = 0)),root = 0)
        rhizoWAfter_ = r.getWaterVolumesCyl(doSum = False, reOrder = True)#np.array([ sum(cylW) for cylW in rhizoWAfter__])
        rhizoWAfter = sum(rhizoWAfter_) 
        #rhizoTotCAfter__= comm.bcast(r._flat0(comm.gather(np.array([ r.getTotCContent(i, cc,r.seg_length_[i]) for i, cc in enumerate(r.cyls)], dtype=object), root = 0)),root = 0) 
        rhizoTotCAfter_= r.getTotCContentAll(doSum = False, reOrder = True)#np.array([ sum(cylC) for cylC in rhizoTotCAfter__]) 
        rhizoTotCAfter = sum(rhizoTotCAfter_)  
        
        r.rhizoMassCError_abs = abs(rhizoTotCAfter - ( rhizoTotCBefore + sum(Q_Exud) + sum(Q_mucil) + sum(proposed_outer_sol_fluxes) *dt+ sum(proposed_outer_mucil_fluxes)*dt))
        r.rhizoMassCError_rel = abs(r.rhizoMassCError_abs/rhizoTotCAfter*100)
        
        
        errorsEachW = rhizoWAfter_ - ( rhizoWBefore_ + (seg_fluxes + proposed_outer_fluxes)*dt)
        errorsEachC = rhizoTotCAfter_ - ( rhizoTotCBefore_ + (seg_sol_fluxes+ proposed_outer_sol_fluxes+ seg_mucil_fluxes+ proposed_outer_mucil_fluxes)*dt)
        try:
            #assert (errorsEachW[airSegsId] == 0.).all()
            assert (errorsEachC[airSegsId] == 0.).all()
        except:
            print("errors in airSegsId",airSegsId)
            print(errorsEachC[airSegsId] )
            print("fluxes water", seg_fluxes[airSegsId], proposed_outer_fluxes[airSegsId]) 
            print(seg_sol_fluxes[airSegsId], proposed_outer_sol_fluxes[airSegsId],seg_mucil_fluxes[airSegsId], proposed_outer_mucil_fluxes[airSegsId])
            raise Exception
            
        r.rhizoMassWError_abs = sum(abs(errorsEachW[rhizoSegsId]))
        r.rhizoMassWError_rel = abs(r.rhizoMassWError_abs/sum(rhizoWAfter_[np.where(organTypes == 2)])*100)
        print("rhizoMassCError_abs, rel", r.rhizoMassCError_abs,r.rhizoMassCError_rel  )
        if (r.rhizoMassCError_abs > 1e-9) and (r.rhizoMassCError_rel > 1e-9):
            maxdiffrhizoC = np.where(errorsEachC == max(abs(errorsEachC)))
            print("rhizoMassCError_abs, rel", r.rhizoMassCError_abs,r.rhizoMassCError_rel  )
            print(errorsEachC, sum(abs(errorsEachC[np.where(organTypes == 2)])))
            print(maxdiffrhizoC, dt)
            print(rhizoTotCAfter__[maxdiffrhizoC],rhizoTotCBefore__[maxdiffrhizoC], seg_sol_fluxes, 
                proposed_outer_sol_fluxes[maxdiffrhizoC], seg_mucil_fluxes[maxdiffrhizoC], proposed_outer_mucil_fluxes[maxdiffrhizoC])
            
            
            raise Exception
        # cyl = r.cyls[1]
        # write_file_array("pressureHead",np.array(cyl.getSolutionHead()).flatten())
        # write_file_array("coord", cyl.getDofCoordinates().flatten())
        # for i in range(r.numFluidComp):
            # write_file_array("solute_conc"+str(i+1), np.array(cyl.getSolution_(i+1)).flatten()* r.molarDensityWat ) 
        # for i in range(r.numFluidComp, r.numComp):
            # write_file_array("solute_conc"+str(i+1), np.array(cyl.getSolution_(i+1)).flatten()* r.bulkDensity_m3 /1e6 ) 

                    

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
        
        if rank == 0:        
            soil_fluxes = rs.sumSegFluxes(seg_fluxes)  # [cm3/day]  per soil cell
        else:
            soil_fluxes = None
        
        
        soil_fluxes = comm.bcast(soil_fluxes, root=0)
        
        """ 3.d  some checks """
        print('time checkMassOMoleBalance2',t)
        r.checkMassOMoleBalance2(soil_fluxes, soil_source_sol, dt,seg_fluxes =seg_fluxes, doSolid = True)
        r.setSoilData(soil_fluxes, soil_source_sol, dt)
        """ 2.0  global soil models """
        
        water_content =comm.bcast( np.array(s.getWaterContent()),root= 0)  # theta per cell [1]
        soil_water = np.multiply(water_content, cell_volumes)  # water per cell [cm3]
        
        soil_solute_content = comm.bcast(np.array([np.array(s.getContent(i+1, isDissolved = (i < r.numFluidComp))) for i in range(r.numComp)]),root = 0)
        #mol
        
        #if rank == 0:#or? no I think should be for all theads as we have "if (sourceMap.count(gIdx)>0)"
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
                
        
        buTotCBefore = comm.bcast(sum(s.getTotCContent()), root = 0) 
        
        start_time_3ds = timeit.default_timer()
        print("solve 3d soil")
        s.solve(dt)  # in modules/solverbase.py
        print("done")
        
        rs.time_3ds_i = (timeit.default_timer() - start_time_3ds)
    
        
        buTotCAfter =comm.bcast( sum(s.getTotCContent()) , root = 0)    
        
        s.bulkMassErrorPlant_abs = abs(buTotCAfter - ( buTotCBefore + sum(Q_Exud) + sum(Q_mucil)))
        s.bulkMassErrorPlant_rel = abs(s.bulkMassErrorPlant_abs/buTotCAfter*100)
        s.bulkMassError1ds_abs = abs(buTotCAfter - ( buTotCBefore + sum(soil_source_sol.flatten())*dt))
        s.bulkMassError1ds_rel = abs(s.bulkMassError1ds_abs/buTotCAfter*100)
        
        print("errorCbalance soil 3d?",buTotCAfter , buTotCBefore , sum(Q_Exud) , sum(Q_mucil), 
                sum(soil_source_sol.flatten())*dt,s.bulkMassErrorPlant_abs, s.bulkMassError1ds_abs )
        #raise Exception
        assert (s.getSolution(0) == r.soilModel.getSolution(0)).all()
        assert (s.getSolution(8) > 0. ).all()# idComp 8 is the respiration. if a voxel has [C_8] == 0, is means there were no biochemical reactions.

        """ 3b. calculate net fluxes """
        water_content =comm.bcast( np.array(s.getWaterContent()), root = 0) 
        # print(water_content)
        new_soil_water = np.multiply(water_content, cell_volumes)  # calculate net flux
        s.bulkMassErrorWater_abs = abs(sum(new_soil_water) - (sum(soil_water) + sum(soil_fluxes.values())*dt))
        s.bulkMassErrorWater_rel = abs(s.bulkMassErrorWater_abs /sum(new_soil_water) )*100
        try:
            # use soil_fluxes and not seg_fluxes. seg_fluxes includes air segments. sum(seg_fluxes) ~ 0.
            # maybe 0.1% of error is too large
            assert s.bulkMassErrorWater_rel < 0.1
        except:
            print(new_soil_water , soil_water , seg_fluxes*dt)
            print(sum(new_soil_water) , sum(soil_water) , sum(soil_fluxes.values())*dt)
            print(s.bulkMassErrorWater_abs)
            print( s.bulkMassErrorWater_rel )
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
        
        soil_solute_content_new = comm.bcast(np.array([np.array(s.getContent(i+1, isDissolved = (i < r.numFluidComp))) for i in range(r.numComp)]), root = 0) # mol
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
        for nc in range(r.numFluidComp, r.numComp):
            # all changes in cellIds for 3D soil is cause by biochemical reactions computed in 1D models.
            # thus, for elements which do not flow (range(r.numComp, r.numFluidComp)), there are no changes
            # except those prescribed by the 1d model.
            try:
                assert (abs(outer_R_bc_sol[nc][cellIds]) < 1e-16).all()
            except:
                print("outer_R_bc_sol[nc][cellIds] != 0.", nc+1, cellIds)
                print(outer_R_bc_sol[nc][cellIds])
                print(soil_solute_content_new[nc][cellIds] , soil_solute_content[nc][cellIds] , soil_source_sol[nc][cellIds]*dt)
                raise Exception
            # outer_R_bc_sol[nc][cellIds] == flow in 3D model
            # outer_R_bc_sol[nc][not cellIds] == flow and reaction in 3D model
        
        r.checkSoilBC(outer_R_bc_wat, outer_R_bc_sol)# check oldSoil + rhizoSource = newSoil - (bulkFlow + bulkReaction )
        # then (bulkFlow + bulkReaction ) sent as 1ds outer BC at next time step.
        if rank == 0:
            write_file_array("sumdiffSoilData_abs_A", r.sumdiffSoilData_abs, directory_ =r.results_dir)# cumulative (?)
            write_file_array("maxdiffSoilData_abs_A", r.maxdiffSoilData_abs, directory_ =r.results_dir)# cumulative (?)
            write_file_array("sumdiffSoilData_rel_A", r.sumdiffSoilData_rel, directory_ =r.results_dir)# cumulative (?)
            write_file_array("maxdiffSoilData_rel_A", r.maxdiffSoilData_rel, directory_ =r.results_dir)# cumulative (?)
            write_file_array("diffSoilData_abs_A", r.diffSoilData_abs, directory_ =r.results_dir)# cumulative (?)
            write_file_array("diffSoilData_rel_A", r.diffSoilData_rel, directory_ =r.results_dir)# cumulative (?)
            
        
        try:
            assert s.bulkMassErrorPlant_rel < 1.e-5
            assert (r.sumdiffSoilData_rel < 1.).all()#e-10
        except:
            print(soil_source_sol)
            print("buTotCAfter ,  buTotCBefore", buTotCAfter ,  buTotCBefore)
            print( "sum(Q_Exud) , sum(Q_mucil)", sum(Q_Exud) , sum(Q_mucil))
            print("3ds source",sum(soil_source_sol.flatten()))
            print(s.bulkMassErrorPlant_abs ,s.bulkMassErrorPlant_rel)
            print(s.bulkMassError1ds_abs ,s.bulkMassError1ds_rel)
            print(r.rhizoMassCError_abs, r.rhizoMassCError_rel)
            print(r.sumdiffSoilData_abs, r.sumdiffSoilData_rel)
            print(r.maxdiffSoilData_abs, r.maxdiffSoilData_rel)
            print("rhizoErrors:")
            maxdiffrhizoC = np.where(abs(errorsEachC) > 1e-9)
            print("rhizoMassCError_abs, rel", r.rhizoMassCError_abs,r.rhizoMassCError_rel  )
            print("errorsEachC",errorsEachC, sum(abs(errorsEachC[np.where(organTypes == 2)])), organTypes[maxdiffrhizoC])
            cyls = np.array(r.cyls)
            print("maxdiffrhizoC",maxdiffrhizoC, dt, cyls[maxdiffrhizoC])
            print(rhizoTotCAfter__[maxdiffrhizoC],rhizoTotCBefore__[maxdiffrhizoC], seg_sol_fluxes, 
                proposed_outer_sol_fluxes[maxdiffrhizoC], seg_mucil_fluxes[maxdiffrhizoC], proposed_outer_mucil_fluxes[maxdiffrhizoC])
            print("cyls",[[dd.getSolution(i + 1) for i in range(r.numComp)] for dd in cyls[maxdiffrhizoC]])
            raise Exception
            
        
        emptySoilVoxels = []
        if rank == 0:
            emptySoilVoxels = np.array([elCid for elCid in range(len(cell_volumes)) if ((elCid >= 0) and (elCid not in cellIds))])
            assert len(emptySoilVoxels) + len(cellIds) == len(cell_volumes)
        r.setEmptySoilVoxel(emptySoilVoxels)# now oldSoil = newSoil where we have no rhizo
        
        totWatAdded = outer_R_bc_wat[emptySoilVoxels] 
        totSolAdded = outer_R_bc_sol[:,emptySoilVoxels] 
        outer_R_bc_wat[emptySoilVoxels] = 0. #changes were put in old soil and will not be used as BC
        outer_R_bc_sol[:,emptySoilVoxels] = 0. #changes were put in old soil and will not be used as BC
        
        r.checkSoilBC(outer_R_bc_wat, outer_R_bc_sol)# check oldSoil + rhizoSource = newSoil - (bulkFlow + bulkReaction )
        
        if rank == 0:
            try:
                assert (r.sumdiffSoilData_rel < 1).all()#e-10
            except:
                print("buTotCAfter ,  buTotCBefore", buTotCAfter ,  buTotCBefore)
                print( "sum(Q_Exud) , sum(Q_mucil)", sum(Q_Exud) , sum(Q_mucil))
                print(s.bulkMassErrorPlant_abs ,s.bulkMassErrorPlant_rel)
                print(r.rhizoMassCError_abs, r.rhizoMassCError_rel)
                print(r.sumdiffSoilData_abs, r.sumdiffSoilData_rel)
                print(r.maxdiffSoilData_abs, r.maxdiffSoilData_rel)
                raise Exception
            
            for nc in range(r.numFluidComp, r.numComp):
                try:
                    assert (abs(outer_R_bc_sol[nc]) < 1e-16).all()
                except:
                    print("outer_R_bc_sol[nc][cellIds] != 0.",rank, nc+1, cellIds)
                    print(outer_R_bc_sol[nc])
                    print(soil_solute_content_new[nc] , soil_solute_content[nc] , soil_source_sol[nc]*dt)
                    raise Exception
        """ 3d. backup """
        
        """ some checks """
        
        try:
            assert min(water_content) >=0
            assert min(soil_solute_content_new.flatten()) >=0
        except:
            print("min(new_soil_solute), min(soil_water)",min(water_content),[min(ssc) for ssc in soil_solute_content])
            raise Exception
            
        """ remember results ... """
        sx = s.getSolutionHead()
        cc = s.getSolution(1)  # [kg/m3]
        if rank == 0:
            sink = np.zeros(sx.shape)
            for k, v in soil_fluxes.items():
                sink[k] += v
            sink_.append(sink)  # cm3/day (per soil cell)
            x_.append(rs_age + t)  # day
            y_.append(np.sum(sink))  # cm3/day
            c_.append(-np.sum(seg_sol_fluxes))  # [cm3/day]
            c_All.append(seg_sol_fluxes)  # [cm3/day]
            c_All1.append(seg_mucil_fluxes)
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
            print("time", rs_age + t, "sx", np.min(sx), np.max(sx), "ccx", np.min(cc), np.max(cc))#, "trans",sum(np.array(r.Ev)))
            psi_x_.append(rx.copy())  # cm (per root node)
            psi_s_.append(sx.copy())  # cm (per root segment)


    #if rank == 0:
    #    print ("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s")

    return psi_x_, psi_s_, sink_, x_, y_, psi_s2_, vol_, surf_,  depth_, soil_c_, c_, c_All,c_All1, outer_R_bc_sol, outer_R_bc_wat 

