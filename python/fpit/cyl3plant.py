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

from scenario_setup import write_file_array, write_file_float


def simulate_const(s, rs, sim_time, dt, rs_age, Q_plant,
                    r = [], wilting_point =-15000, 
                    outer_R_bc_sol=[], #mol
                    outer_R_bc_wat = [], 
                    results_dir = './results/'):#m3 , directory_ =results_dir
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
    cell_volumes = comm.bcast(s.getCellVolumes() , root = 0) #cm3
    
    airSegsId = np.array(list(set(np.concatenate((r.cell2seg[-1],np.where(np.array(r.organTypes) != 2)[0])) )))#aboveground
    rhizoSegsId = np.array([i for i in range(len(organTypes)) if i not in airSegsId])
    Q_Exud = Q_plant[0]; Q_mucil = Q_plant[1] #mol/day
    if len(Q_Exud) > 0:
        try:
            assert (Q_Exud[airSegsId] == 0).all()
            assert (Q_mucil[airSegsId] == 0).all()
        except:
            print(Q_Exud[airSegsId] )
            print(Q_mucil[airSegsId])
            print(organTypes[airSegsId])
            raise Exception
    else:
        Q_Exud = np.full(len(organTypes),0.)
        Q_mucil = np.full(len(organTypes),0.)

    if(len(outer_R_bc_sol[0]) == 0):
        outer_R_bc_sol = np.full((r.numComp,len(cell_volumes)), 0.)  
    if(len(outer_R_bc_wat) == 0):
        outer_R_bc_wat = np.full(cell_volumes.shape, 0.)        
            
    assert len(rs.rs.segments) == (len(rs.rs.nodes) -1)
    seg2cell = rs.rs.seg2cell
    cell2seg = rs.rs.cell2seg
    cellIds = np.fromiter(cell2seg.keys(), dtype=int)
    cellIds =  np.array([x for x in cellIds if x >= 0])#take out air segments
    

    
    N = 1#int(np.ceil(sim_time / dt))  # number of iterations
    
    """ simualtion loop """
    for i in range(N):

        t = rs_age + i * dt  # current simulation time
        
        weatherX = weather(t)

        
        
        comm.barrier()
            
        rs.Qlight = weatherX["Qlight"]
        rs.Csoil_seg = r.get_inner_solutes() * 1e3 # mol/cm3 to mmol/cm3 
        rx_old = 0
        new_soil_water_old = 0
        rhizoWAfter_old = np.full(len(organTypes),0.)
        n_iter = 0
        err = 1.e6 
        max_err = 1e-13 # ???
        max_iter = 100 #??
        while err > max_err and n_iter < max_iter:
            
            """ 1. xylem model """
            print('1. xylem model')
            
            ##
            # 1.1 plant water flow and photosynthesis
            ##
            
            # send soil concentration to plant:
            rsx = r.get_inner_heads(weather=weatherX)  # matric potential at the root soil interface, i.e. inner values of the cylindric models (not extrapolation to the interface!) [cm]
            
            start_time_plant = timeit.default_timer()
            comm.barrier()
            
            if (rank == 0):
                assert min(rs.Csoil_seg ) >= 0.
                try:
                    rs.solve_photosynthesis(sim_time_ = rs_age, 
                                sxx_=rsx, 
                                cells_ = False,#(i == 0),#for 1st computation, use cell data
                                ea_ = weatherX["ea"],#not used
                                es_=weatherX["es"],#not used
                                verbose_ = False, doLog_ = False,
                                TairC_= weatherX["TairC"],#not used
                                outputDir_= "./results/rhizoplantExud")
                except:
                    rs.solve_photosynthesis(sim_time_ = rs_age, 
                                sxx_=rsx, 
                                cells_ = False,#(i == 0),#for 1st computation, use cell data
                                ea_ = weatherX["ea"],#not used
                                es_=weatherX["es"],#not used
                                verbose_ = True, doLog_ = True,
                                TairC_= weatherX["TairC"],#not used
                                outputDir_= "./results/rhizoplantExud")
                    
                
            comm.barrier()
            
            ##
            # 1.2 get data (unlimited fluxes)
            ##
            
            seg_fluxes = np.array(rs.outputFlux)# [cm3/day] 
            seg_sol_fluxes = Q_Exud /dt# mol/day for segments
            seg_mucil_fluxes = Q_mucil/dt
            
            rs.time_plant_cumulW += (timeit.default_timer() - start_time_plant)
            
            if (rank == 0):
                soil_fluxes_ = rs.sumSegFluxes(seg_fluxes)  # [cm3/day]  per soil cell
                soil_fluxes = np.zeros(len(cell_volumes))
                print('soil_fluxes_',soil_fluxes_.keys(),  soil_fluxes_.items(),  soil_fluxes_.values(),  soil_fluxes_ )
                soil_fluxes[np.array(list(soil_fluxes_.keys()))] = np.array(list(soil_fluxes_.values())) #easier to handle array than dict. maybe change sumSegFluxes() function to choose type of output
            else:
                soil_fluxes = None
            soil_fluxes = comm.bcast(soil_fluxes, root=0)
            rs.outputFlux = comm.bcast(rs.outputFlux, root = 0) 
            
            
            rs.psiXyl = comm.bcast(rs.psiXyl, root = 0) 
            
            assert min(np.concatenate((seg_mucil_fluxes,seg_sol_fluxes))) >= 0. #currently, no net plant solute uptake
                
            """ 2. local 1D soil models (1DS)"""
            print('2. local 1D soil models (1DS)')
            
            ##
            # 2.1 distribute 3D flows between the 1DS
            #     use value per 1DS !!AT THE END OF THE TIME STEP!! => weight for @splitSoilVals()
            ##
            waterContent = r.getWaterVolumesCyl(doSum = False, reOrder = True)
            comp1content = r.getContentCyl(idComp=1, doSum = False, reOrder = True)
            comp2content = r.getContentCyl(idComp=2, doSum = False, reOrder = True)
            
            try:
                assert (waterContent[airSegsId] == 0.).all()
                assert (comp1content[airSegsId] == 0.).all()
                assert (comp2content[airSegsId] == 0.).all()
                assert waterContent.shape == (len(organTypes), )
            except:
                print(waterContent,waterContent.shape)
                raise Exception
            
            if rank == 0:
                if max(abs(outer_R_bc_wat )) > 0:
                    assert outer_R_bc_wat.shape == ( len(cell_volumes), )
                    proposed_outer_fluxes = r.splitSoilVals(outer_R_bc_wat / dt, waterContent) #cm3/day
                else:
                    proposed_outer_fluxes = np.full(len(organTypes), 0.)                
                if max(abs(outer_R_bc_sol[0] )) > 0:
                    proposed_outer_sol_fluxes = r.splitSoilVals(outer_R_bc_sol[0] / dt, comp1content)#mol/day
                else:
                    proposed_outer_sol_fluxes = np.full(len(organTypes), 0.)
                if max(abs(outer_R_bc_sol[1] )) > 0:
                    proposed_outer_mucil_fluxes = r.splitSoilVals(outer_R_bc_sol[1] / dt, comp2content)
                else:
                    proposed_outer_mucil_fluxes = np.full(len(organTypes), 0.)
            else:
                proposed_outer_fluxes = None
                proposed_outer_sol_fluxes = None
                proposed_outer_mucil_fluxes = None
            proposed_outer_fluxes = comm.bcast(proposed_outer_fluxes, root = 0)
            proposed_outer_sol_fluxes = comm.bcast(proposed_outer_sol_fluxes, root = 0)
            proposed_outer_mucil_fluxes = comm.bcast(proposed_outer_mucil_fluxes, root = 0)
            
            
            assert (np.array([len(seg_fluxes), len(proposed_outer_fluxes),len(seg_sol_fluxes),
                               len(proposed_outer_sol_fluxes), len(seg_mucil_fluxes),
                                len(proposed_outer_mucil_fluxes)]) == len(organTypes)).all()
            
            
            ##
            # 2.2 data before solve, for post proccessing
            # maybe move this part to within the solve function to go less often through the list of 1DS
            ##
            print('do reset')
            if n_iter > 0:
                r.reset() # go back to water and solute value at the BEGINING of the time step
            print('did reset')
            rhizoWBefore_ = r.getWaterVolumesCyl(doSum = False, reOrder = True)
            rhizoWBefore = sum(rhizoWBefore_) 
            
            rhizoTotCBefore_ = r.getTotCContentAll(doSum = False, reOrder = True)#np.array([ sum(cylC) for cylC in rhizoTotCBefore__]) 
            rhizoTotCBefore = sum(rhizoTotCBefore_) 
            
            soil_solute = np.array( [np.array(r.getC_rhizo(len(cell_volumes), idComp = idc + 1, konz = False)) for idc in range(r.numComp)])
            
            ##
            # 2.3A limit the net negative BCs
            # update the INNER BC for water as will be < 0  normally
            # update the OUTER BCs for solutes as will < 0 normally (necessary/TODO?)
            # TODO: necessary to limit the outer solute BC by the BC content? because we can also gain by reactions
            # for now no limit on the solute fluxes
            ##
            innerSurf = 2 * np.pi * r.seg_length * r.radii
            outerSurf = 2 * np.pi * r.seg_length * r.outer_radii
            Q_outer_totW = proposed_outer_fluxes * outerSurf * dt
            seg_fluxes_limited = np.maximum(seg_fluxes, -(rhizoWBefore_ - r.vg_soil.theta_R + Q_outer_totW )/(dt*innerSurf))
            seg_fluxes_limited[airSegsId] = seg_fluxes[airSegsId] # limitation only relevent for root segments belowground
            
            # proposed_outer_sol_fluxes_limited = np.maximum(proposed_outer_sol_fluxes, -(comp1content + Q_outer_totW )/(dt*innerSurf))
            # proposed_outer_mucil_fluxes_limited = np.maximum(proposed_outer_mucil_fluxes, -(comp2content + Q_outer_totW )/(dt*innerSurf))
            
            r.SinkLim1DS = abs(seg_fluxes_limited - seg_fluxes) # remember the error caused by the limitation
            
            ##
            # 2.3B simulation
            ##
            start_time_rhizo = timeit.default_timer()
            print("solve 1d soil", rank)
            
            r.solve(dt, n_iter,seg_fluxes_limited, proposed_outer_fluxes, seg_sol_fluxes,proposed_outer_sol_fluxes, 
                        seg_mucil_fluxes, proposed_outer_mucil_fluxes) # cm3/day or mol/day
            
            rs.time_rhizo_i += (timeit.default_timer() - start_time_rhizo)
        
            ##
            # 2.4 data after, for post proccessing
            # maybe move this part to within the solve function to go less often through the list of 1DS
            ##
            
            rhizoWAfter_ = r.getWaterVolumesCyl(doSum = False, reOrder = True)# per 1ds
            rhizoWAfter = sum(rhizoWAfter_) 
            
            rhizoTotCAfter_= r.getTotCContentAll(doSum = False, reOrder = True)# per 1ds 
            rhizoTotCAfter = sum(rhizoTotCAfter_)  
            
            
            ##
            # 2.5 error rates
            ##
            r.rhizoMassCError_abs = abs(rhizoTotCAfter - ( rhizoTotCBefore + sum(Q_Exud) + sum(Q_mucil) + sum(proposed_outer_sol_fluxes) *dt+ sum(proposed_outer_mucil_fluxes)*dt))
            r.rhizoMassCError_rel = abs(r.rhizoMassCError_abs/rhizoTotCAfter*100)
            
            # instantaneous mass water error for 1DS (with seg_fluxes_limited)
            errorsEachW = rhizoWAfter_ - ( rhizoWBefore_ + (seg_fluxes_limited + proposed_outer_fluxes)*dt)
            errorsEachW[airSegsId] = 0 # error evaluation only adapted for root segments belowground            
            r.rhizoMassWError_abs = sum(abs(errorsEachW[rhizoSegsId]))
            r.rhizoMassWError_rel = abs(r.rhizoMassWError_abs/sum(rhizoWAfter_[rhizoSegsId])*100)
            print("for limited flux: rhizoMassWError_abs, rel", r.rhizoMassWError_abs,r.rhizoMassWError_rel)
            if (r.rhizoMassWError_abs > 1e-9) or (r.rhizoMassCError_abs > 1e-9):# check directly for seg_fluxes_limited
                maxdiffrhizoC = np.where(errorsEachC == max(abs(errorsEachC)))
                print("rhizoMassWError_abs, rel", r.rhizoMassWError_abs,r.rhizoMassWError_rel  )
                print("rhizoMassCError_abs, rel", r.rhizoMassCError_abs,r.rhizoMassCError_rel  )
                print(errorsEachC, sum(abs(errorsEachC[np.where(organTypes == 2)])))
                print(maxdiffrhizoC, dt)
                raise Exception
                
            # mass water error to check when leaving fixed point iteration (with seg_fluxes)
            errorsEachW = rhizoWAfter_ - ( rhizoWBefore_ + (seg_fluxes + proposed_outer_fluxes)*dt)
            errorsEachW[airSegsId] = 0 # error evaluation only adapted for root segments belowground            
            r.rhizoMassWError_abs = sum(abs(errorsEachW[rhizoSegsId]))
            r.rhizoMassWError_rel = abs(r.rhizoMassWError_abs/sum(rhizoWAfter_[rhizoSegsId])*100)
            print("for proposed flux: rhizoMassWError_abs, rel", r.rhizoMassWError_abs,r.rhizoMassWError_rel, 
                    max(abs(seg_fluxes-seg_fluxes_limited)))
                
            errorsEachC = rhizoTotCAfter_ - ( rhizoTotCBefore_ + (seg_sol_fluxes+ proposed_outer_sol_fluxes+ seg_mucil_fluxes+ proposed_outer_mucil_fluxes)*dt)
            try:
                assert (errorsEachC[airSegsId] == 0.).all()
            except:
                print("errors in airSegsId",airSegsId)
                print(errorsEachC[airSegsId] )
                print("fluxes water", seg_fluxes[airSegsId], proposed_outer_fluxes[airSegsId]) 
                print(seg_sol_fluxes[airSegsId], proposed_outer_sol_fluxes[airSegsId],seg_mucil_fluxes[airSegsId], proposed_outer_mucil_fluxes[airSegsId])
                raise Exception
                
                         
            
            ##
            # 2.6 calculate initial 3DS net sources
            ##
            
            # mol per voxel, TODO: find a way to use one function for that and rhizoTotCAfter_
            new_soil_solute = np.array( [np.array(r.getC_rhizo(len(cell_volumes), idComp = idc + 1, konz = False)) for idc in range(r.numComp)])
            try:
                assert min(new_soil_solute.flatten()) >=0
            except:
                print("min(new_soil_solute)",min(new_soil_solute.flatten()),[min(nss) for nss in new_soil_solute])
                raise Exception
            
            
            #mol/day
            soil_source_sol = np.full(new_soil_solute.shape,0. )
            for nc in range(r.numComp):
                soil_source_sol[nc][cellIds] = np.array(new_soil_solute[nc][cellIds] - soil_solute[nc][cellIds] - outer_R_bc_sol[nc][cellIds])/dt

            soil_source_sol = comm.bcast(soil_source_sol, root = 0)    
            assert soil_source_sol.shape == (r.numComp, len(cell_volumes))

            
            
            """ 3. global soil models (3DS)"""
            
            ##
            # 3.1 data before, for post proccessing AND source adaptation
            ##
            if n_iter > 0:
                s.reset() #reset at the last moment: over functions use the solution/content at the end of the time step
            water_content = comm.bcast( np.array(s.getWaterContent()),root= 0)  # theta per cell [1]
            soil_water = np.multiply(water_content, cell_volumes)  # water per cell [cm3]
            buTotCBefore = comm.bcast(sum(s.getTotCContent()), root = 0) 
            
            soil_solute_content = comm.bcast(np.array([np.array(s.getContent(i+1, isDissolved = (i < r.numFluidComp))) for i in range(r.numComp)]),
                                                root = 0) # mol

            
            ##
            # 3.2 adapt and set sources
            # TODO: take into account BC != 0
            # TODO: currently, only limit water flow. not lear how to limit solute flow
            ##          
            lower_lim = np.ones(soil_fluxes.shape) * s.vg_soil.theta_R
            soil_fluxes_limited = np.maximum(soil_fluxes, -(water_content - lower_lim)/dt) # update to limit the net sinks 
            
            soil_sources_limited = np.concatenate((np.array([soil_fluxes_limited]),soil_source_sol ))
            soil_contents = np.concatenate((np.array([water_content]),soil_solute_content )) 
            assert soil_sources_limited.shape == (s.numComp+1, len(cell_volumes) )
            assert soil_contents.shape == (s.numComp+1, len(cell_volumes))
            
            for idComp in range(s.numComp+1):#mol/day
                if (max(abs(soil_sources_limited[idComp])) != 0.):
                    test_values = list(soil_sources_limited[idComp].copy())
                    test_keys = np.array([i for i in range(len(test_values))])
                    res = {}
                    for key in test_keys:
                        for value in test_values:
                            res[key] = value
                            test_values.remove(value)
                            break                        
                    write_file_float("setsource_"+str(idComp)+"_"+str(max_rank), res, directory_ =results_dir) 
                    s.setSource(res.copy(), eq_idx = idComp)  # [mol/day], in modules/richards.py
            r.SinkLim3DS = abs(soil_fluxes_limited - soil_fluxes) # at the end of the fixed point iteration, should be ~ 0        
            
            ##
            # 3.3 solve 3DS
            ##    
            
            start_time_3ds = timeit.default_timer()
            print("solve 3d soil")
            s.solve(dt)  # in modules/solverbase.py
            print("done")
            rs.time_3ds_i += (timeit.default_timer() - start_time_3ds)
            
            ##
            # 3.4 data after, for post proccessing 
            ##
            buTotCAfter =comm.bcast( sum(s.getTotCContent()) , root = 0)    
            water_content =comm.bcast( np.array(s.getWaterContent()), root = 0) 
            new_soil_water = np.multiply(water_content, cell_volumes)  # calculate net flux
            
            
            ##
            # 3.5 error rates
            ##
            s.bulkMassErrorWater_abs = abs(sum(new_soil_water) - (sum(soil_water) + sum(soil_fluxes)*dt))
            s.bulkMassErrorWater_rel = abs(s.bulkMassErrorWater_abs /sum(new_soil_water) )*100
            
            s.bulkMassErrorPlant_abs = abs(buTotCAfter - ( buTotCBefore + sum(Q_Exud) + sum(Q_mucil)))
            s.bulkMassErrorPlant_rel = abs(s.bulkMassErrorPlant_abs/buTotCAfter*100)
            s.bulkMassError1ds_abs = abs(buTotCAfter - ( buTotCBefore + sum(soil_source_sol.flatten())*dt))
            s.bulkMassError1ds_rel = abs(s.bulkMassError1ds_abs/buTotCAfter*100)
            
            print("errorCbalance soil 3d?",rank, buTotCAfter ,',', buTotCBefore ,',',  sum(Q_Exud) ,',',  sum(Q_mucil), 
                    ', soil_source_sol', sum(soil_source_sol.flatten())*dt,', s.bulkMassErrorPlant_abs', s.bulkMassErrorPlant_abs,
                    ', bulkMassError1ds_abs ',  s.bulkMassError1ds_abs , t)
            s.buTotCAfter = buTotCAfter
            s.buTotCBefore = buTotCBefore
            
            print("errorWbalance soil 3d?",rank, sum(new_soil_water) ,',', sum(soil_water) ,',',   sum(soil_fluxes)*dt,
                        ', bulkMassErrorWater_abs', s.bulkMassErrorWater_abs,', bulkMassErrorWater_rel',  s.bulkMassErrorWater_rel , t)
            
            ##
            # 3.6 get 1DS outer BC
            ##
            outer_R_bc_wat = new_soil_water - soil_water - soil_fluxes *dt # change in water per cell [cm3]
            soil_solute_content_new = comm.bcast(np.array([np.array(s.getContent(i+1, isDissolved = (i < r.numFluidComp))) for i in range(r.numComp)]), root = 0) # mol
            outer_R_bc_sol = np.array([soil_solute_content_new[i] - soil_solute_content[i] - soil_source_sol[i]*dt for i in range(r.numComp)])# mol
            assert outer_R_bc_sol.shape == (r.numComp, len(cell_volumes))

                
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
            
            
            """ 4. prints and evaluation of the iteration """
            # TODO: see which value are more sensible. very low rx in the leaves may make the level of change always low in the plant
            
            r.checkMassOMoleBalance2(soil_fluxes*0, soil_source_sol*0, dt,
                                    seg_fluxes =seg_fluxes*0, diff1d3dCW_rel_lim = np.Inf)
            rx = np.array(rs.psiXyl)
            errRxPlant = np.linalg.norm(rx - rx_old)
            errW1ds = np.linalg.norm(rhizoWAfter_[rhizoSegsId] - rhizoWAfter_old[rhizoSegsId])
            
            errW3ds = np.linalg.norm(new_soil_water - new_soil_water_old)
            errs = np.concatenate((np.array([errRxPlant, errW1ds, errW3ds]), r.SinkLim3DS,r.SinkLim1DS ))
            err = max(r.maxDiff1d3dCW_abs)
            
            rx_old = rx.copy()
            new_soil_water_old = new_soil_water.copy()
            rhizoWAfter_old = rhizoWAfter_.copy()
            
            if rank == 0:
                write_file_array("errorIteration"+str(max_rank), errs, directory_ =results_dir) 
                write_file_array("errorIteration1d3d"+str(max_rank), r.maxDiff1d3dCW_abs, directory_ =results_dir) 
                write_file_float("errorIteration_time"+str(max_rank), rs_age, directory_ =results_dir) 
                
            if False:#rank == 0:    
                write_file_float("soil_fluxes_"+str(max_rank), soil_fluxes, directory_ =results_dir) 
                write_file_array("new_soil_water_"+str(max_rank), new_soil_water, directory_ =results_dir) 
                write_file_array("soil_water_"+str(max_rank), soil_water, directory_ =results_dir) 
                write_file_array("outer_R_bc_wat_B"+str(max_rank), outer_R_bc_wat, directory_ =results_dir) 
            
            n_iter += 1
            #end iteration
        
        ####
        #   error rates    
        ####
        
        # error mass 1DS
        if (r.rhizoMassWError_abs > 1e-9) or (r.rhizoMassCError_abs > 1e-9):
            maxdiffrhizoC = np.where(errorsEachC == max(abs(errorsEachC)))
            print("rhizoMassWError_abs, rel", r.rhizoMassWError_abs,r.rhizoMassWError_rel  )
            print("rhizoMassCError_abs, rel", r.rhizoMassCError_abs,r.rhizoMassCError_rel  )
            print(errorsEachC, sum(abs(errorsEachC[np.where(organTypes == 2)])))
            print(maxdiffrhizoC, dt)
            raise Exception

        # error 3DS-1DS
        print('soil_fluxes',soil_fluxes, soil_source_sol, dt)
        r.checkMassOMoleBalance2(soil_fluxes*0, soil_source_sol*0, dt,
                                seg_fluxes =seg_fluxes*0)
                                
        if rank == 0:            
            for nc in range(r.numFluidComp, r.numComp):
                try:
                    assert (abs(outer_R_bc_sol[nc]) < 1e-16).all()
                except:
                    print("outer_R_bc_sol[nc][cellIds] != 0.",rank, nc+1, cellIds)
                    print(outer_R_bc_sol[nc])
                    print(soil_solute_content_new[nc] , soil_solute_content[nc] , soil_source_sol[nc]*dt)
                    raise Exception
                    
        # error mass plant    
        try:
            assert s.bulkMassErrorPlant_abs < 1.e-5
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
            
        # error mass 3DS
    
        try:
            # use soil_fluxes and not seg_fluxes. seg_fluxes includes air segments. sum(seg_fluxes) ~ 0.
            # maybe 0.1% of error is too large
            assert s.bulkMassErrorWater_abs < 1e-9
        except:
            print(new_soil_water , soil_water , seg_fluxes*dt)
            print(sum(new_soil_water) , sum(soil_water) , sum(soil_fluxes.values())*dt)
            print(s.bulkMassErrorWater_abs)
            print( s.bulkMassErrorWater_rel )
            print(new_soil_water)
            raise Exception
            
        try:
            assert (abs(((new_soil_water - soil_water ) - (soil_fluxes*dt + outer_R_bc_wat ))/((new_soil_water - soil_water ) )) < 0.1).all()
        except:
            print(new_soil_water - (soil_water +soil_fluxes* dt + outer_R_bc_wat ))
            raise Exception
            
        # end time step inner loop

    return outer_R_bc_sol, outer_R_bc_wat # first guess for next fixed point iteration
    #end of inner loop

