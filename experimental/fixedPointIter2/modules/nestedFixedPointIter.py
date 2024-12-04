"""
functions for the cylindrical coupling approach
"""
import numpy as np
import timeit
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank(); max_rank = comm.Get_size()
import sys
from functional.xylem_flux import *
from air_modelsPlant import AirSegment
import visualisation.vtk_plot as vtk
import plantbox as pb
#import evapotranspiration as evap
import timeit
import visualisation.vtk_plot as vp
import functional.van_genuchten as vg
from decimal import *
from functional.xylem_flux import sinusoidal2, sinusoidal
import helpfull
from helpfull import write_file_array, write_file_float
from FPItHelper import fixedPointIterationHelper
import FPItHelper
import PhloemPhotosynthesis
import printData
from functional.xylem_flux import sinusoidal2

from weatherFunctions import weather, weatherChange


def innerLoop(plantModel,rs_age, fpit_Helper, perirhizalModel , sim_time, dt, s):

    results_dir = perirhizalModel.results_dir
    N = int(np.ceil(sim_time / dt))  # number of iterations
    n_iter_inner_max = 0
    for Ni in range(N):
        fpit_Helper.n_iter3 = 0
        rs_age_i_dt = rs_age + Ni * dt  # current simulation time
        keepGoing = True
        # FPItHelper.storeOldMassData1d(perirhizalModel)
        # print("perirhizalModel.rhizoWBefore_eachCyl",perirhizalModel.rhizoWBefore_eachCyl)
        while keepGoing:
                                
            
            perirhizalModel.solve_gave_up = False # we reset solve_gave_up (== did dumux fail?) each time

            """ 1. xylem model """
            ##
            # 1.1 plant water flow (and photosynthesis)
            ##
            
            
            plantModel.time_start_plant = timeit.default_timer() # to compute time spent in water flow simulation
            
            
            PhloemPhotosynthesis.computeWaterFlow(fpit_Helper, perirhizalModel, plantModel, rs_age_i_dt, dt)
            
            # store time spent in computing plant water flow
            plantModel.time_plant_cumulW += (timeit.default_timer() - plantModel.time_start_plant)
            

            """ 2. local 1D soil models (1DS)"""
            ##
            # 2.1 distribute 3D flows between the 1DS
            #     use value per 1DS !!AT THE END OF THE TIME STEP!! => weight for @splitSoilVals()
            ##
            
            fpit_Helper.distribute3dto1dFlows(rs_age_i_dt, dt)
            
            
            ##
            # 2.2 reset 1d models + store data before solve, for post processing
            # 
            ##
            
            if ((fpit_Helper.n_iter3 > 0)|(fpit_Helper.n_iter2 > 0)|(
                    fpit_Helper.n_iter > 0)):
                #if rank == 0:
                #    print("\t\tdoreset")
                perirhizalModel.reset() # go back to water and solute value at the BEGINING of the time step
            else:
                FPItHelper.storeOldMassData1d(perirhizalModel)
                    
                
            

            ##
            # 2.3 simulation 1d models
            ##
            plantModel.time_start_rhizo = timeit.default_timer()

            
            if (perirhizalModel.mpiVerbose or (max_rank == 1)) and rank == 0:
                    print("\t\tsolve all 1d soils ")
                    
            perirhizalModel.solve(dt, 
                                  fpit_Helper.seg_fluxes , # inner BC water
                                  fpit_Helper.proposed_outer_fluxes, # outer BC water
                                  fpit_Helper.seg_sol_fluxes, # inner BC solute 1
                                  fpit_Helper.proposed_outer_sol_fluxes, # outer BC solute 1
                                  fpit_Helper.seg_mucil_fluxes, # inner BC solute 2
                                  fpit_Helper.proposed_outer_mucil_fluxes, # outer BC
                                  # solute 2
                                  # fpit_Helper.n_iter
                                 ) # cm3/day or mol/day
            
            if (perirhizalModel.mpiVerbose or (max_rank == 1)) and rank == 0:
                    print("\t\tsolve all 1d soils finished")
                    
            plantModel.time_rhizo_i += (timeit.default_timer() - plantModel.time_start_rhizo)
            
            fpit_Helper.storeLimitedFluxes() # store actual BC of 1d models (to compare with prescribed BC)
            
            
            
            ##
            # 2.4 data after, for post processing
            ##
            
            FPItHelper.storeNewMassData1d(perirhizalModel)
             
            
            
            ##
            # 2.5 error rates
            ##
            fpit_Helper.massBalanceError1d(dt)
            fpit_Helper.computeConvergence2()


            # print extra data for troubleshooting
            # TODO: finish adapting the name of the objects to print
            printData.printFPitDataInnerLoop(perirhizalModel, plantModel, fpit_Helper, rs_age_i_dt)
            ###########
            #       leave??
            ##########
            keepGoing, gaveUp = helpfull.continueLoop(perirhizalModel, n_iter=
                                fpit_Helper.n_iter3,
                                                 dt_inner= dt,
                         failedLoop = False, real_dt =float(Ni) * dt,doPrint = False,
                                                      name='inner_loopdata',
                         isInner=True,
                         plant=plantModel, FPIT_id = 3)
            if not perirhizalModel.doNestedFixedPointIter:
                keepGoing, gaveUp = False, False
            fpit_Helper.n_iter3 += 1

        n_iter_inner_max = max(n_iter_inner_max,fpit_Helper.n_iter3)
            
        # store transpiration and assimilation data    
        if (rank == 0):
            # root -soil exchange per root segment for water and solute 1 and 2
            plantModel.seg_fluxes0Cumul_inner += perirhizalModel.seg_fluxes_limited * dt
            plantModel.seg_fluxes1Cumul_inner += perirhizalModel.seg_fluxes_limited_sol_In * dt
            plantModel.seg_fluxes2Cumul_inner += perirhizalModel.seg_fluxes_limited_mucil_In * dt

            if perirhizalModel.doPhotosynthesis:
                plantModel.TranspirationCumul_inner += sum(np.array(plantModel.Ev) * dt) #transpiration [cm3/day] * day
                plantModel.AnCumul_inner += np.array(plantModel.An ) * (dt*24*3600) # //[mol CO2 m-2 s-1] to [mol CO2 m-2]
                    
                write_file_float("N_Transpiration_inner", sum(np.array(plantModel.Ev) * dt), directory_ =results_dir)
                
                write_file_array("N_Q_Ag_dot_inner", plantModel.AnCumul_inner/ (dt*24*3600), directory_ =results_dir)
                write_file_array("N_Q_Ag_dot_innerbis", plantModel.An, directory_ =results_dir)
            else:
                rootSegs = np.array(perirhizalModel.organTypes) == pb.root
                plantModel.TranspirationCumul_inner += sum(fpit_Helper.seg_fluxes[rootSegs]) * dt
                
            write_file_float("N_TranspirationCumul_inner", plantModel.TranspirationCumul_inner, directory_ =results_dir)
        
        dt_inner = float(Ni + 1) * float(dt)
        return gaveUp, dt_inner, n_iter_inner_max


def simulate_const(s, plantModel, sim_time, dt, rs_age, 
                   Q_plant,
                    perirhizalModel = [],
                    outer_R_bc_sol=[], #mol
                    outer_R_bc_wat = [], # cm3 water
                   seg_fluxes=[],
                   doMinimumPrint=True):
    """     
    simulates the coupled scenario       
        root architecture is not growing  
        conductivities are not changing over time
        
    s: soil model
    plantModel: MAppedPlant or (TODO adapt for ) MappedRoot object
    sim_time: total duration of the fixed point iteration [d]
    dt: time step of the fixed-point iteration [d]
    rs_age: plant age [d]
    Q_plant: plant solute flow (if exudation/mucilage release simulation)
    perirhizalModel: object containing all the 1d models
    computed at last time fixed-point iteration (first guess for new loop):
        outer_R_bc_sol: inter-cell flow of solutes in 3d soil
        outer_R_bc_wat: inter-cell flow of water in 3d soil
        seg_fluxes: plant-soil water exchanges (used to compute soil kr)
    doMinimumPrint: print extra information to file? (for debugging)
    """
    ######
    ###              store data which are constant during iteration loop
    ######
    results_dir = s.results_dir
    cell_volumes = s.getCellVolumes() #cm3    
    N = int(np.ceil(sim_time / dt))  # number of iterations
     
    # max num of iteration needed for the fixed-point iteration loop. used by @see helpfull::suggestNumStepsChange
    n_iter_inner_max = 0  

    ### plant shape
    organTypes = np.asarray(perirhizalModel.organTypes, int)
    cylVol = perirhizalModel.getVolumesCyl(doSum = False, reOrder = True) # volume of perirhizal zones [cm3]
    if (rank == 0) and perirhizalModel.debugMode:
        print("cylVol",cylVol)
        print("perirhizalModel.rhizoVol",perirhizalModel.rhizoVol)
        print("perirhizalModel.eidx_all_",perirhizalModel.eidx_all_)
        
    airSegsId = perirhizalModel.airSegs # id of plant segments WITHOUT perirhizal zones
    rhizoSegsId = perirhizalModel.rhizoSegsId # id of plant segments WITH perirhizal zones   
    cellIds = perirhizalModel.cellWithRoots # only id of cells with roots
    emptyCells = np.array(list(set([xsoil for xsoil in range(s.numberOfCellsTot)]) - set(cellIds))) # -1 (air) not included
    if (not doMinimumPrint) or perirhizalModel.debugMode:
        # for debugging
        write_file_array('emptyCells',emptyCells,directory_ =results_dir, fileType = '.csv')
        write_file_array('cellWithRoots',cellIds,directory_ =results_dir, fileType = '.csv')
    
    ######
    ###              flow, 1st guess for iteration loop
    ######
    # get root exudation and mucilage release
    Q_Exud_i = Q_plant[0].copy(); Q_mucil_i = Q_plant[1].copy() #mol  
    if len(Q_Exud_i) > 0: # resize vector if needed
        Q_Exud_i.resize(len(organTypes), refcheck=False) #, refcheck=False for cProfile
        Q_mucil_i.resize(len(organTypes), refcheck=False)
        try:
            # make sure that no exudation in shoot segments or roots aboveground
            assert (Q_Exud_i[airSegsId] == 0).all()
            assert (Q_mucil_i[airSegsId] == 0).all()
        except:
            print(Q_Exud_i[airSegsId] )
            print(Q_mucil_i[airSegsId])
            print(organTypes[airSegsId])
            raise Exception
    else:# first loop: create array of 0
        Q_Exud_i = np.full(len(organTypes),0.)
        Q_mucil_i = np.full(len(organTypes),0.)

    # first loop: create array of 0
    if(len(outer_R_bc_sol[0]) == 0):
        outer_R_bc_sol = np.full((perirhizalModel.numSoluteComp,s.numberOfCellsTot), 0.)  
    if(len(outer_R_bc_wat) == 0):
        outer_R_bc_wat = np.full(cell_volumes.shape, 0.)     
        
        
        
    """ simualtion loop """
    for Ni in range(N):

        rs_age_i_dt = rs_age + Ni * dt  # current simulation time
        
        perirhizalModel.weatherX = weather(simDuration = rs_age_i_dt, dt = dt,
                                           hp =  max([tempnode[2] for tempnode in plantModel.get_nodes()]) /100., #canopy height [m]
                                           spellData= perirhizalModel.spellData)
        if rank==0:
            # transpiration = plantModel.maxTranspiration * min(rs_age_i_dt/plantModel.maxTranspirationAge,1.) *sinusoidal2(rs_age_i_dt, dt) # just for printing: it is recomputed during @see computeWaterFlow()
            # , transpiration: {transpiration:.2e} cm3/d
            print(f"\n\ninner loop step: {Ni}/{N}. current simulation time: {rs_age_i_dt:.2f} day, Qlight: {perirhizalModel.weatherX['Qlight']:.2e} cm3/d")
        
        
        
        weatherChange(rs_age_i_dt, perirhizalModel, s) # implement sudden change in temperature, soil wat. content ext...
        
        if perirhizalModel.doPhotosynthesis: # data needed for photosynthesis
            PhloemPhotosynthesis.computeAtmosphereData(plantModel, perirhizalModel)
            
        ####
        #       reset inner loop data (create XXX_old vals to evaluate convergence)
        ####
        # simple class to store fixed-point iteration data and some usefull functions
        fpit_Helper = fixedPointIterationHelper(s, perirhizalModel, plantModel, 
                                                seg_fluxes, outer_R_bc_wat, 
                                                outer_R_bc_sol, cylVol, Q_Exud_i, Q_mucil_i,
                                                dt, # sub-time step
                                                sim_time, # sim_time = N * dt
                                                emptyCells)

        # looping 1
        keepGoing = True
        while keepGoing:

            # looping 2
            #dt_inner = min(fpit_Helper.dt_inner2, fpit_Helper.dt_inner)
            keepGoingInner = True
            fpit_Helper.n_iter2 = 0
                # nested inner loop
            while keepGoingInner :
                assert perirhizalModel.dt_inner2 <= dt

                #if rank == 0:
                #    print("enter inner loop")
                failedInnerLoop, real_dt, n_iter_inner_max2= innerLoop(
                    plantModel = plantModel, rs_age = rs_age_i_dt,
                    fpit_Helper  =fpit_Helper,perirhizalModel = perirhizalModel ,
                    sim_time = dt, dt = perirhizalModel.dt_inner2, s=s)
                #if rank == 0:
                #    print("leave inner loop")
                keepGoingInner = (failedInnerLoop & (
                            fpit_Helper.n_iter2 < perirhizalModel.k_iter))
            
                # we leave the loop either because dumux threw an error or because
                # we reached dt_inner
                try:
                    assert ((abs(real_dt - dt) < perirhizalModel.dt_inner2 ) or
                            failedInnerLoop)
                except:
                    print('real_dtinner',real_dt ,dt, perirhizalModel.dt_inner2 ,
                          failedInnerLoop)
                    write_file_array("real_dtinner_error", np.array([dt,real_dt ,
                                                                     perirhizalModel.dt_inner2 ,
                                                                     failedInnerLoop,rs_age_i_dt,Ni]),
                                     directory_ =results_dir, fileType = '.csv') 
                    raise Exception
                    
                if perirhizalModel.doNestedFixedPointIter: 
                    fpit_Helper.n_iter2 += 1
                        
                    perirhizalModel.dt_inner2 = helpfull.suggestNumStepsChange(dt,
                                                                perirhizalModel.dt_inner2 ,
                                                                  failedInnerLoop, perirhizalModel,
                                                                  n_iter_inner_max2)
                
                if keepGoingInner:
                    
                    if rank==0:
                        print("error too high, restart loop with dt_inner2",
                              perirhizalModel.dt_inner2)
                    # reset data to the beginning of the iteration loops
                    helpfull.resetPlantWFlux(plantModel, perirhizalModel)
                    
            assert not failedInnerLoop

                
            """ 3. global soil models (3DS)"""
            
            ##
            # 3.1 data before, for post proccessing AND source adaptation
            ##
            
            
            
            if (fpit_Helper.n_iter > 0) :
                s.reset() #reset at the last moment: over functions use the solution/content at the end of the time step
            else:
                FPItHelper.storeOldMassData3d(s,perirhizalModel)
            
            
            
            #######
            #
            # 1d1dFlow (do after the reset of the soil data)
            #
            #######
            
            if perirhizalModel.do1d1dFlow:
                fpit_Helper.do1d1dFlow()
            
            ##
            # 3.2 adapt and set sources for 3d soil using data from 1d models
            # TODO: take into account BC != 0 (rain, evaporation)
            ##          
            
            
            fpit_Helper.compute3DSource(dt)
            
            
            ##
            # 3.3 solve 3DS
            ##    
            
            plantModel.time_start_3ds = timeit.default_timer()
            if (perirhizalModel.mpiVerbose ) and (rank == 0):
                print("solve 3d soil" )
                
            fpit_Helper.solve3DS(dt)
            
            
            if (perirhizalModel.mpiVerbose ) and (rank == 0):
                print("solve 3d soil finished" )
                
            plantModel.time_3ds_i += (timeit.default_timer() - plantModel.time_start_3ds)
            
            ##
            # 3.4 data after, for post proccessing 
            ##
            
            
            FPItHelper.storeNewMassData3d(s,perirhizalModel)
            
            perirhizalModel.check1d3dDiff() # just to get error value, will not throw an error
            
            
            

            ##
            # 3.5.A get 1DS outer BC (from limited-flux: used at next iteration)
            ##
            
            # get flux and source data directly from dumux. 
            # TODO: do the same to get directly change rate of 1d model
            outer_R_bc = s.getFlux_10c() # < 0 means net sink, > 0 means net source
            bulkSoil_sources = s.getSource_10c() # < 0 means net sink, > 0 means net source
            
            if rank == 0:
                fpit_Helper.outer_R_bc_wat = outer_R_bc[0]# [cm3] 
                fpit_Helper.sources_wat_from3d =  bulkSoil_sources[0]# cm3

                fpit_Helper.outer_R_bc_sol = outer_R_bc[1:] # mol
                fpit_Helper.sources_sol_from3d =  bulkSoil_sources[1:] # mol

                assert fpit_Helper.outer_R_bc_sol.shape == (perirhizalModel.numSoluteComp, s.numberOfCellsTot)
            ##
            # 3.6 mass or content balance error 3d 
            ##
            fpit_Helper.massBalanceError3d(dt)
                    
            
            
            
            """ 4. prints and evaluation of the iteration """
            
            #perirhizalModel.check1d3dDiff() # just to get error value, will not throw an error
            
            fpit_Helper.computeConvergence() # evaluate convergence and get other error metrics
            
            # print extra data for troubleshooting
            # TODO: finish adapting the name of the objects to print
            printData.printFPitData(perirhizalModel, s, plantModel, fpit_Helper, rs_age_i_dt)
            
            
            ##
            # 3.5.B set 1DS outer BC to 0 where needed
            ##
            
            for nc in range(perirhizalModel.numDissolvedSoluteComp, perirhizalModel.numSoluteComp):
                fpit_Helper.outer_R_bc_sol[nc][:] = 0.# to not have balance error as 3d flux
                # all changes in cellIds for 3D soil is cause by biochemical reactions computed in 1D models.
                # thus, for elements which do not flow (range(perirhizalModel.numDissolvedSoluteComp, perirhizalModel.numSoluteComp)), 
                # there are no changes
                # except those prescribed by the 1d model.
            
             

            fpit_Helper.n_iter += 1
            # store max number of iteration which occured in the inner iteration loop
            # for adaptation of the next inner time step (@see helpfull::suggestNumStepsChange())
            n_iter_inner_max = max(n_iter_inner_max,fpit_Helper.n_iter)


            # for debug
            if False:#not perirhizalModel.doMinimumPrint:
                datas = []
                datasName = [ ]
                if rank == 0:   
                    # to check where a specific (problematic) segment is
                    #is569 = np.array([i for i in range(len(fpit_Helper.rsx_old))]) == 569 
                    datas = [#is569,
                             plantModel.psiXyl, 
                              fpit_Helper.rsx_old, 
                              np.array(fpit_Helper.seg_fluxes), 
                             fpit_Helper.proposed_outer_fluxes
                            ]
                    datasName = [ #'is569',
                                "psiXyl",
                                 "rsi", "flux_in", 
                                 "flux_out"
                                ]
                printData.doVTPplots(str(int(rs_age*10))+"_"+str(fpit_Helper.n_iter), #indx of vtp plot
                                    perirhizalModel, plantModel,s, perirhizalModel.getSoilTextureAndShape(), 
                                    datas, datasName, initPrint=False, doSolutes = perirhizalModel.doSoluteFlow)
            ###########
            #       leave??
            ##########
            keepGoing, failedLoop = helpfull.continueLoop(perirhizalModel, n_iter=
                                fpit_Helper.n_iter,
                                                 dt_inner= dt,
                         failedLoop = False, real_dt =float(Ni) * dt,
                                                          name='inner_loopdata',
                         isInner=True,
                         plant=plantModel)
        ####
        #   error rates    
        ####
    
        # => NO need to recompute it here I think
        # error 3DS-1DS
        # perirhizalModel.check1d3dDiff()
        
        
        
        # did we leave the inner loop because n_iter == k_iter (failedLoop = True)?
        if (failedLoop):# no need to go on, leave inner loop now and reset lower time step
            if rank == 0:
                print('Failed, no need to go on, leave inner loop now and reset lower time step')
            break
        # end outer loop
            
    
    dt_inner =float(Ni+1)*float( dt) # get the real simetime if sim_time / dt != int
    
    # fluxes == first guess for next fixed point iteration
    return outer_R_bc_sol, fpit_Helper.outer_R_bc_wat, fpit_Helper.seg_fluxes, dt_inner, failedLoop, n_iter_inner_max


