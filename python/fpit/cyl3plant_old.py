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
import functional.van_genuchten as vg
from scenario_setup import weather
from decimal import *
from functional.xylem_flux import sinusoidal2, sinusoidal

from scenario_setup import write_file_array, write_file_float


def continueLoop_(rs,n_iter, dt_inner,failedLoop,real_dtinner, name):
    sumDiff1d3dCW_rel = rs.sumDiff1d3dCW_rel[:(rs.numFluidComp+1)]
    sumDiff1d3dCW_rel = np.where(np.isnan(sumDiff1d3dCW_rel),0.,sumDiff1d3dCW_rel)
    cL = ((np.floor(rs.err) > max_err) or (abs(rs.rhizoMassWError_abs) > 1e-13) or (abs(rs.rhizoMassCError_abs) > 1e-9) or (max(abs(rs.errDiffBCs)) > 1.) or  rs.solve_gave_up or (max(abs(sumDiff1d3dCW_rel))>1)) and (n_iter < max_iter)

    comm.barrier()
    print('continue loop?',rank,cL)
    comm.barrier()
    return cL
        
def simulate_const(s, rs, sim_time, dt, rs_age, Q_plant,
                    r = [],
                    outer_R_bc_sol=[], #mol
                    outer_R_bc_wat = [], seg_fluxes=[],
                    results_dir = './results/',
                  adaptRSI  = True, plantType = "plant",
                  k_iter_ = 100,lightType_ = "", outer_n_iter = 0,continueLoop =continueLoop_):
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

    if r.mpiVerbose or (max_rank == 1):
        comm.barrier()
        print("cyl3plant:cell_volumes", rank)
        comm.barrier()
    cell_volumes = comm.bcast(s.getCellVolumes() , root = 0) #cm3

    if r.mpiVerbose or (max_rank == 1):
        comm.barrier()
        print("cyl3plant:GOTcell_volumes", rank)
        comm.barrier()
    airSegsId = r.airSegs
    #np.array(list(set(np.concatenate((np.array([]),#to avoid error if there are no air segs
    #                                              r.cell2seg.get(-1),
    #                                              np.where(np.array(r.organTypes) != 2)[0])) )))#aboveground
    rhizoSegsId = r.rhizoSegsId
    # np.array([i for i in range(len(organTypes)) if i not in airSegsId])
    
    # check that rhizoSegsId and airSegsId are as expected
    local_isRootSeg = np.array([not isinstance(cyl,AirSegment) for cyl in np.array(r.cyls)])
    global_isRootSeg = r.getXcyl(local_isRootSeg, doSum = False, reOrder = True)
    assert (global_isRootSeg[rhizoSegsId]).all()
    
    Q_Exud = Q_plant[0].copy(); Q_mucil = Q_plant[1].copy() #mol/day
    if len(Q_Exud) > 0:

        Q_Exud.resize(len(organTypes))
        Q_mucil.resize(len(organTypes))
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
    

    
    N = int(np.ceil(sim_time / dt))  # number of iterations
    skip = 1  # 3 * 6  # for output and results, skip iteration
    
    """ simualtion loop """
    for Ni in range(N):

        rs_age_i_dt = rs_age + Ni * dt  # current simulation time
        
        
        if r.mpiVerbose or (max_rank == 1):
            comm.barrier()
            print("simualtion loop, Ni, N:",Ni, N,rs_age_i_dt)
            comm.barrier()
        r.weatherX = weather(rs_age_i_dt)

        
        
        comm.barrier()
        if lightType_ == "":
            rs.Qlight = r.weatherX["Qlight"]
        elif lightType_ == "nolight":
            rs.Qlight = 0.
        else:
            raise Exception
        rs.Csoil_seg = r.get_inner_solutes() * 1e3 # mol/cm3 to mmol/cm3 
        solution0_3ds_old = np.array(s.getSolution(0))
        solution0_1ds_old = np.array([cyl.getSolution(0) for cyl in r.cyls],dtype=object)
        
        rx_old = 0
        new_soil_water_old = 0
        soil_solute_content_new_old = 0
        rhizoWAfter_old = np.full(len(organTypes),0.)
        rhizoTotCAfter_old = np.full(len(organTypes),0.)
        seg_fluxes_old = 0 
        proposed_outer_fluxes_old = 0
        proposed_outer_sol_fluxes_old = 0
        proposed_outer_mucil_fluxes_old = 0
        outer_R_bc_wat_old =  outer_R_bc_wat.copy()
        outer_R_bc_sol_old = outer_R_bc_sol.copy()
        n_iter = 0
        r.err = 1.e6 
        max_err = 1.# ???
        max_iter = k_iter_#100 #??
        rsx_set = r.get_inner_heads(weather=r.weatherX)# matric potential at the segment-exterior interface, i.e. inner values of the (air or soil) cylindric models 
        rsx_old = rsx_set.copy()
        r.sumDiff1d3dCW_absBU = r.sumDiff1d3dCW_abs # to go from cumulative to instantenuous 1d3d error
        r.sumDiff1d3dCW_relBU = r.sumDiff1d3dCW_rel # to go from cumulative to instantenuous 1d3d error
        waterContentOld = r.getWaterVolumesCyl(doSum = False, reOrder = True)
        # weightBefore = True
        rsx_old_  = r.get_inner_heads(weather=r.weatherX) 
        r.rhizoMassWError_abs =1.# 
        r.rhizoMassCError_abs =1.# 
        r.errDiffBCs = np.array([1.])
        r.solve_gave_up = False
        r.diff1d3dCurrant_rel =1e6
        
        while  continueLoop(r,n_iter, dt,False,Ni * dt,'fpit_loopdata'):
            #( (np.floor(r.err) > max_err) or (abs(r.rhizoMassWError_abs) > 1e-13) 
            #or (abs(r.rhizoMassCError_abs) > 1e-9) or (max(abs(r.errDiffBCs)) > 1.) 
            #or  r.solve_gave_up or (max(abs(sumDiff1d3dCW_rel))>1)) and (n_iter < max_iter) :
            r.solve_gave_up = False
            """ 1. xylem model """
            
            if r.mpiVerbose or (max_rank == 1):
                comm.barrier()
                print('1. xylem model', rank)
                comm.barrier()
            
            ##
            # 1.1 plant water flow and photosynthesis
            ##
            
            # send soil concentration to plant:
            
            start_time_plant = timeit.default_timer()
            comm.barrier()

            if r.SRIBefore or (r.beforeAtNight and (r.weatherX["Qlight"] == 0.)) :
                rsx_set = rsx_old_

            # s.vg_soil.Ksat [cm/d] conductivity AT SATURATION. decrease with wat. content, via krw value 
            krw = r.getKrw()# - 
            deltaR_cm =  r.getDeltaR()
            if r.l_ks == "dx_2":
                dist_factor =  deltaR_cm
            elif r.l_ks == "dx":
                dist_factor = deltaR_cm * 2
            elif r.l_ks == "root":
                dist_factor =  np.array(r.radii)# * 100.
            else:
                raise Exception
                
            if rank == 0:
                #soilKIn = Ksat * krw  / deltaR_cm # [cm/d] * [-] / [cm] = day-1
                soilKIn =np.divide(vg.hydraulic_conductivity(rsx_set, r.vg_soil),dist_factor)
                #np.array(r.radii) * 100.) # y use the root radius?
                soilKOut = s.vg_soil.Ksat  /dist_factor# (np.array(r.radii) * 100.)#deltaR_cm # [cm/d]  / [cm]  = day-1
                soilK = soilKIn
                
                if( len(seg_fluxes) > 0.) and not (r.SRIBefore  or (r.beforeAtNight and (r.weatherX["Qlight"] == 0.))):
                    soilK[np.where(seg_fluxes> 0. ) ] = soilKOut[np.where(seg_fluxes > 0. ) ]
                    
                if len(r.airSegs) > 0:   
                    soilK[r.airSegs] = np.Inf # works when sending it to C++?
                # write_file_array("fpit_deltaR_cm", deltaR_cm, directory_ =results_dir, fileType = '.csv') 
                write_file_array("fpit_krw", krw, directory_ =results_dir, fileType = '.csv') 
                write_file_array("fpit_dist_factor", dist_factor, directory_ =results_dir, fileType = '.csv') 
                write_file_array("fpit_soilK", soilK, directory_ =results_dir, fileType = '.csv') 
                #print('soilK',soilK,'soilKIn',soilKIn,'soilKOut', soilKOut,'perimeter_cm',perimeter_cm,
                #      'krw',krw,'Ksat',Ksat,s.vg_soil.Ksat)
                # soil_k = np.divide(vg.hydraulic_conductivity(rsx_set, r.vg_soil), r.radii) # y use the root radius? 
                write_file_array("fpit_soilKrwKs", krw*s.vg_soil.Ksat, directory_ =results_dir, fileType = '.csv') 
                write_file_array("fpit_soilKrwKsBis", vg.hydraulic_conductivity(rsx_set, r.vg_soil), directory_ =results_dir, fileType = '.csv') 
                # write_file_array("fpit_soilKBis", soil_k, directory_ =results_dir, fileType = '.csv') 
            if( plantType == "plant") and (rank == 0):
                
                try:                    
                    assert min(rs.Csoil_seg ) >= 0.

                    write_file_array("fpit_rsxUsed",np.array(rsx_set),directory_ =results_dir, fileType = '.csv')
                    write_file_float("fpit_weatherX",r.weatherX,directory_ =results_dir)
                    rs.solve_photosynthesis(sim_time_ = rs_age_i_dt, 
                                sxx_=rsx_set, 
                                cells_ = False,#(i == 0),#for 1st computation, use cell data
                                ea_ = r.weatherX["ea"],#not used
                                es_=r.weatherX["es"],#not used
                                verbose_ = False, doLog_ = False,
                                TairC_= r.weatherX["TairC"],#not used
                                            soil_k_ = soilK, # [day-1]
                                outputDir_= "./results/rhizoplantExud")
                    seg_fluxes = np.array(rs.outputFlux)# [cm3/day] 
                    TransRate = sum(np.array(rs.Ev)) #transpiration [cm3/day] 
                    write_file_array("fpit_Ev",np.array(rs.Ev),directory_ =results_dir, fileType = '.csv')
                    write_file_array("fpit_Jw",np.array(rs.Jw),directory_ =results_dir, fileType = '.csv')
                    write_file_array("fpit_fw",np.array(rs.fw),directory_ =results_dir, fileType = '.csv')#pg
                    write_file_array("fpit_pg",np.array(rs.pg),directory_ =results_dir, fileType = '.csv')
                    write_file_array("fpit_n_iter",np.array([ n_iter,rs.loop , r.solve_gave_up]), directory_ =results_dir, fileType = '.csv') 

                    write_file_array('fpit_transRate',np.array([TransRate,TransRate*dt]), directory_ =results_dir, fileType = '.csv' )
                    write_file_array("fpit_errPhoto", np.array(rs.maxErr) , directory_ =results_dir, fileType = '.csv') 
                    write_file_array("fpit_errPhotoAbs", np.array(rs.maxErrAbs) , directory_ =results_dir, fileType = '.csv') 
                    write_file_array("fpit_organTypes", organTypes, directory_ =results_dir, fileType = '.csv') 
            
                    if (plantType == "plant") and (rank == 0):
                        leavesSegs = np.where(organTypes ==4)
                        fluxes_leaves = seg_fluxes[leavesSegs]
                        if (min(rs.Ev) < 0) or (min(rs.Jw) < 0) or (min(fluxes_leaves)<-1e-15):
                            print("leaf looses water", min(rs.Ev),min(rs.Jw), min(fluxes_leaves))
                            print("seg_fluxes",seg_fluxes,"leavesSegs", leavesSegs)                
                            raise Exception
                except:
                    rs.minLoop = 2
                    rs.maxLoop = 5
                    rs.solve_photosynthesis(sim_time_ = rs_age_i_dt, 
                                sxx_=rsx_set, 
                                cells_ = False,#(i == 0),#for 1st computation, use cell data
                                ea_ = r.weatherX["ea"],#not used
                                es_=r.weatherX["es"],#not used
                                verbose_ = True, doLog_ = True,
                                TairC_= r.weatherX["TairC"],#not used
                                            soil_k_ = soilK, # [day-1]
                                outputDir_= ".")
                    raise Exception

            elif (rank == 0):
                transpiration = 6. *  sinusoidal2(rs_age, dt)
                rx = rs.solve(rs_age, 
                             -transpiration, 0., 
                             rsx_set, cells = False, 
                              soil_k = [])
                rs.psiXyl = rx
                seg_fluxes = np.array(rs.segFluxes(rs_age_i_dt, rx, rsx_set, 
                                                   False, False, #approx, cells
                                                   []))
                rs.outputFlux = seg_fluxes
            elif( plantType == "plant") and (rank > 0):
                seg_fluxes = None
            else :
                rs.psiXyl = None
                rs.outputFlux = None
                seg_fluxes = None
            
            comm.barrier()
            
            if (plantType == "plant") and (rank == 0):
                leavesSegs = np.where(organTypes ==4)
                fluxes_leaves = seg_fluxes[leavesSegs]
                if (min(rs.Ev) < 0) or (min(rs.Jw) < 0) or (min(fluxes_leaves)<-1e-15):
                    print("leaf looses water", min(rs.Ev),min(rs.Jw), min(fluxes_leaves))
                    print("seg_fluxes",seg_fluxes,"leavesSegs", leavesSegs)                
                    raise Exception
            
            ##
            # 1.2 get data (unlimited fluxes)
            ##
            
            
            seg_sol_fluxes = Q_Exud /dt# mol/day for segments
            seg_mucil_fluxes = Q_mucil/dt
            
            rs.time_plant_cumulW += (timeit.default_timer() - start_time_plant)

            
            if r.mpiVerbose or (max_rank == 1):
                comm.barrier()
                print("cyl3plant:seg_fluxes", rank)
                comm.barrier()
            seg_fluxes = comm.bcast(seg_fluxes, root=0)
            rs.outputFlux = comm.bcast(rs.outputFlux, root = 0) 
            
            
            rs.psiXyl = comm.bcast(rs.psiXyl, root = 0) 
            
            if r.mpiVerbose or (max_rank == 1):
                comm.barrier()
                print("cyl3plant:GOTseg_fluxes", rank)
                comm.barrier()
            assert min(np.concatenate((seg_mucil_fluxes,seg_sol_fluxes))) >= 0. #currently, no net plant solute uptake
                
            """ 2. local 1D soil models (1DS)"""
            comm.barrier()
            if r.mpiVerbose or (max_rank == 1):
                print(rank, '2. local 1D soil models (1DS)')
            comm.barrier()
            if r.mpiVerbose or (max_rank == 1):
                print(rank,'2.1 distribute 3D flows between the 1DS')
            ##
            # 2.1 distribute 3D flows between the 1DS
            #     use value per 1DS !!AT THE END OF THE TIME STEP!! => weight for @splitSoilVals()
            ##
            if r.mpiVerbose or (max_rank == 1):
                print(rank,'2.1 distribute 3D flows between the 1DS')
            if r.weightBefore  or (r.beforeAtNight and (r.weatherX["Qlight"] == 0.)):
                if r.mpiVerbose or (max_rank == 1):
                    print(rank,'waterContent = waterContentOld')
                waterContent = waterContentOld
            else:
                if r.mpiVerbose or (max_rank == 1):
                    print(rank,'waterContent = r.getWaterVolumesCyl(doSum = False, reOrder = True)')
                waterContent = r.getWaterVolumesCyl(doSum = False, reOrder = True)
            
            if r.mpiVerbose or (max_rank == 1):
                print(rank,'get comp1content')
            comp1content = r.getContentCyl(idComp=1, doSum = False, reOrder = True)
            if r.mpiVerbose or (max_rank == 1):
                print(rank,'get comp2content')
            comp2content = r.getContentCyl(idComp=2, doSum = False, reOrder = True)
            if r.mpiVerbose or (max_rank == 1):
                print(rank,'check comp1content',len(airSegsId))
            
            if len(airSegsId)>0:
                try:
                    assert (waterContent[airSegsId] == 0.).all()
                    assert (comp1content[airSegsId] == 0.).all()
                    assert (comp2content[airSegsId] == 0.).all()
                    assert waterContent.shape == (len(organTypes), )
                    if r.mpiVerbose or (max_rank == 1):
                        print(rank,'2.1 asserts',(waterContent[airSegsId] == 0.).all(),(comp1content[airSegsId] == 0.).all(),
                          (comp2content[airSegsId] == 0.).all(),waterContent.shape == (len(organTypes), ))
                except:
                    print(rank,'len(airSegsId)>0', '(waterContent[airSegsId] != 0.).all()', waterContent,
                          comp1content,comp2content,airSegsId,waterContent.shape)
                    raise Exception
                    
            comm.barrier()
            if r.mpiVerbose or (max_rank == 1):
                print(rank, 'todo: splitSoilValsA')
            comm.barrier()
            if r.mpiVerbose or (max_rank == 1):
                print(rank, 'splitSoilValsB')
            if rank == 0:
                print(rank, 'get proposed_fluxes')
                if max(abs(outer_R_bc_wat )) > 0:
                    assert outer_R_bc_wat.shape == ( len(cell_volumes), )
                    proposed_outer_fluxes = r.splitSoilVals(outer_R_bc_wat / dt, waterContent) #cm3/day
                else:
                    proposed_outer_fluxes = np.full(len(organTypes), 0.)   
                print(rank, 'got proposed_outer_fluxes')
                if max(abs(outer_R_bc_sol[0] )) > 0:
                    proposed_outer_sol_fluxes = r.splitSoilVals(outer_R_bc_sol[0] / dt, comp1content)#mol/day
                else:
                    proposed_outer_sol_fluxes = np.full(len(organTypes), 0.)
                print(rank, 'got proposed_outer_sol_fluxes')
                if max(abs(outer_R_bc_sol[1] )) > 0:
                    proposed_outer_mucil_fluxes = r.splitSoilVals(outer_R_bc_sol[1] / dt, comp2content)
                else:
                    proposed_outer_mucil_fluxes = np.full(len(organTypes), 0.)
                print(rank, 'got proposed_outer_mucil_fluxes')
            else:
                proposed_outer_fluxes = None
                proposed_outer_sol_fluxes = None
                proposed_outer_mucil_fluxes = None
                if r.mpiVerbose or (max_rank == 1):
                    print(rank, 'set proposed_fluxes to none')
            if r.mpiVerbose or (max_rank == 1):
                print(rank,'to comm.barrier')
                comm.barrier()
            if r.mpiVerbose or (max_rank == 1):
                print(rank, 'did_comm.bcast(fluxes)')
                comm.barrier()
            if r.mpiVerbose or (max_rank == 1):
                print(rank,'left comm.barrier')
                
            proposed_outer_fluxes = comm.bcast(proposed_outer_fluxes, root = 0)
            proposed_outer_sol_fluxes = comm.bcast(proposed_outer_sol_fluxes, root = 0)
            proposed_outer_mucil_fluxes = comm.bcast(proposed_outer_mucil_fluxes, root = 0)
            
            
            if r.mpiVerbose or (max_rank == 1):
                comm.barrier()
                print("cyl3plant:GOTproposed_outer_fluxes", rank)
                comm.barrier()
            try:
                assert (np.array([len(seg_fluxes), len(proposed_outer_fluxes),len(seg_sol_fluxes),
                               len(proposed_outer_sol_fluxes), len(seg_mucil_fluxes),
                                len(proposed_outer_mucil_fluxes)]) == len(organTypes)).all()
            except:
                print('issue length arrays')
                print(np.array([len(seg_fluxes), len(proposed_outer_fluxes),len(seg_sol_fluxes),
                               len(proposed_outer_sol_fluxes), len(seg_mucil_fluxes),
                                len(proposed_outer_mucil_fluxes)]) , len(organTypes))
                raise Exception
            
            
            ##
            # 2.2 data before solve, for post proccessing
            # maybe move this part to within the solve function to go less often through the list of 1DS
            ##
            
            if r.mpiVerbose or (max_rank == 1):
                comm.barrier()
                print(rank, '2.2 data before solve, for post proccessing')
                comm.barrier()
                
            if (n_iter > 0) :
                r.reset() # go back to water and solute value at the BEGINING of the time step
                
            solution0_1ds_new = np.array([cyl.getSolution(0) for cyl in r.cyls],dtype=object)
            
            if r.mpiVerbose or (max_rank == 1):
                comm.barrier()
                print('did reset')#,rank,solution0_1ds_new , solution0_1ds_old)
                comm.barrier()
            try:
                diff_solution0_1ds = np.concatenate(np.array([solution0_1ds_new[i] == solution0_1ds_old[i] for i in range(len(solution0_1ds_new))],dtype=object))
                assert diff_solution0_1ds.all()
            except:
                print('issue with assert 1dsnewold',rank)
                print(rank,solution0_1ds_new , solution0_1ds_old)
                print(rank, diff_solution0_1ds)
                print(rank, r.eidx)
                raise Exception
            
            if r.mpiVerbose or (max_rank == 1):
                comm.barrier()
                print('getWaterVolumesCyl')
                comm.barrier()
                
            rhizoWBefore_ = r.getWaterVolumesCyl(doSum = False, reOrder = True) #cm3
            rhizoWBefore = sum(rhizoWBefore_) 
            
            
            if r.mpiVerbose or (max_rank == 1):
                comm.barrier()
                print('getTotCContentAll')
                comm.barrier()
                
            rhizoTotCBefore_ = r.getTotCContentAll(doSum = False, reOrder = True)
            rhizoTotCBefore = sum(rhizoTotCBefore_) 
            
            
            if r.mpiVerbose or (max_rank == 1):
                comm.barrier()
                print('getC_rhizo')
                comm.barrier()
            if r.mpiVerbose or (max_rank == 1):
                print('get soil_solute')
            soil_solute = np.array( [np.array(r.getC_rhizo(len(cell_volumes), idComp = idc + 1, konz = False)) for idc in range(r.numComp)])
            if r.mpiVerbose or (max_rank == 1):
                print('got soil_solute')
                comm.barrier()
            ##
            # 2.3A 1st limit to the net negative BCs
            # update the INNER BC for water as will be < 0  normally
            # update the OUTER BCs for solutes as will < 0 normally (necessary/TODO?)
            # TODO: necessary to limit the outer solute BC by the BC content? because we can also gain by reactions
            # for now no limit on the solute fluxes
            ##
            if r.mpiVerbose or (max_rank == 1):
                print('2.3A 1st limit to the net negative BCs')
            cylVolume =  np.pi *(np.array( r.outer_radii)*np.array( r.outer_radii )- np.array( r.radii) * np.array( r.radii))* np.array( r.seg_length)
            assert ((rhizoWBefore_ - r.vg_soil.theta_R * cylVolume)[rhizoSegsId] >=0).all()
            
            if r.mpiVerbose or (max_rank == 1):
                print('did assert ((rhizoWBefore_ - r.vg_soil.theta_R * cylVolume)[rhizoSegsId] >=0).all()')
                comm.barrier()

            Q_outer_totW = proposed_outer_fluxes * dt
            # seg_fluxes_limited = np.maximum(seg_fluxes, -(rhizoWBefore_ - r.vg_soil.theta_R * cylVolume+ Q_outer_totW )/dt)
            # if len(airSegsId)>0:
            #     seg_fluxes_limited[airSegsId] = seg_fluxes[airSegsId] # limitation only relevent for root segments belowground
            
            #omax = rho_ * krw * kc * ((h - criticalPressure_) / dz )* pos0
            # proposed_outer_sol_fluxes_limited = np.maximum(proposed_outer_sol_fluxes, -(comp1content + Q_outer_totW )/dt)
            # proposed_outer_mucil_fluxes_limited = np.maximum(proposed_outer_mucil_fluxes, -(comp2content + Q_outer_totW )/dt)
            #print(rank, '1rst seg_fluxes_limited',seg_fluxes_limited[rhizoSegsId],'diff', seg_fluxes_limited[rhizoSegsId] - seg_fluxes[rhizoSegsId])
    
            ##
            # 2.3B simulation
            ##
            if r.mpiVerbose or (max_rank == 1):
                print('2.3B simulation')
            start_time_rhizo = timeit.default_timer()
            if r.mpiVerbose or (max_rank == 1):
                comm.barrier()
                print("solve 1d soil", rank)
                comm.barrier()
            #seg_fluxes_limited
            
            r.solve(dt, n_iter,seg_fluxes , proposed_outer_fluxes, seg_sol_fluxes,proposed_outer_sol_fluxes, 
                        seg_mucil_fluxes, proposed_outer_mucil_fluxes) # cm3/day or mol/day
            
            if r.mpiVerbose or (max_rank == 1):
                comm.barrier()
                print("solve 1d soil_finished", rank)
                comm.barrier()
            rs.time_rhizo_i += (timeit.default_timer() - start_time_rhizo)
            for lId, cyl in enumerate(r.cyls):
                if not isinstance(cyl, AirSegment):
                    gId = r.eidx[lId]
                    
                    pHead = np.array(cyl.getSolutionHead()).flatten()
                    write_file_array("watercontentcyl"+str(gId),cyl.getWaterContent(), 
                                     directory_ =results_dir, allranks = True)
                    write_file_array("pressureHeadcyl"+str(gId),pHead, 
                                     directory_ =results_dir, allranks = True)
                    write_file_array("coordcyl"+str(gId), cyl.getDofCoordinates().flatten(), 
                                     directory_ =results_dir, allranks = True)
                    for ccc in range(r.numComp):
                        sol0 = np.array(cyl.getSolution(ccc)).flatten()
                        write_file_array("solution"+str(ccc)+"_"+str(gId)+"", 
                                     sol0, 
                                     directory_ =results_dir, allranks = True)
                    if max(pHead) > 0:
                        print('issue phead',gId,rank, pHead, sol0 )
                        raise Exception

            if r.mpiVerbose or (max_rank == 1):
                comm.barrier()
                print("share seg_fluxes_limited", rank)# get 2nd limitation, gather and cast among thread 
                comm.barrier()
            seg_fluxes_limited = r.getXcyl(r.seg_fluxes_limited, doSum = False, reOrder = True) 
            seg_fluxes_limited_Out = r.getXcyl(r.seg_fluxes_limited_Out, doSum = False, reOrder = True) 
            seg_fluxes_limited_sol_Out = r.getXcyl(r.seg_fluxes_limited_sol_Out, doSum = False, reOrder = True) 
            seg_fluxes_limited_mucil_Out = r.getXcyl(r.seg_fluxes_limited_mucil_Out, doSum = False, reOrder = True) 
            seg_fluxes_limited_sol_In = r.getXcyl(r.seg_fluxes_limited_sol_In, doSum = False, reOrder = True)
            seg_fluxes_limited_mucil_In = r.getXcyl(r.seg_fluxes_limited_mucil_In, doSum = False, reOrder = True) 
            
                
            if len(airSegsId)>0:                
                try:
                    assert (seg_fluxes_limited[airSegsId] == seg_fluxes[airSegsId]).all()
                except:
                    print('seg_fluxes_limited vs seg_flux', seg_fluxes_limited[airSegsId] - seg_fluxes[airSegsId])
                    raise Exception
            
            r.SinkLim1DS = seg_fluxes_limited - seg_fluxes # remember the error caused by the limitation
            r.OutLim1DS = seg_fluxes_limited_Out - proposed_outer_fluxes # remember the error caused by the limitation
            r.InOutBC_Cdiff = np.array([seg_fluxes_limited_sol_Out-proposed_outer_sol_fluxes,
                                            seg_fluxes_limited_mucil_Out-proposed_outer_mucil_fluxes,
                                           seg_fluxes_limited_sol_In-seg_sol_fluxes,
                                            seg_fluxes_limited_mucil_In-seg_mucil_fluxes])
            
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
            # instantaneous mass water error for 1DS (with seg_fluxes_limited)
            errorsEachC = rhizoTotCAfter_ - ( rhizoTotCBefore_ + (seg_fluxes_limited_sol_In+ seg_fluxes_limited_mucil_In+\
                                                                  seg_fluxes_limited_sol_Out+ seg_fluxes_limited_mucil_Out)*dt)
            if len(airSegsId)>0:                
                try:
                    assert (errorsEachC[airSegsId] == 0.).all()
                except:
                    print("errors in airSegsId",airSegsId)
                    print('errorsEachC',errorsEachC[airSegsId] )
                    print("fluxes solutes",seg_fluxes_limited_sol_In[airSegsId], 
                          seg_fluxes_limited_mucil_In[airSegsId],seg_fluxes_limited_sol_Out[airSegsId], seg_fluxes_limited_mucil_Out[airSegsId])
                    print('rhizoTotCAfter_, rhizoTotCBefore_', rhizoTotCAfter_, rhizoTotCBefore_)
                    raise Exception
                    
            r.rhizoMassCError_absLim = sum(abs(errorsEachC[rhizoSegsId]))
            if rhizoTotCAfter != 0:
                r.rhizoMassCError_relLim = abs(r.rhizoMassCError_absLim/rhizoTotCAfter*100)
            else:
                r.rhizoMassCError_relLim = np.nan
            
            # mass water error to check when leaving fixed point iteration (with seg_fluxes)
            errorsEachC = rhizoTotCAfter_ - ( rhizoTotCBefore_ + (seg_sol_fluxes+ proposed_outer_sol_fluxes+ seg_mucil_fluxes+ proposed_outer_mucil_fluxes)*dt)
            if len(airSegsId)>0:                
                try:
                    assert (errorsEachC[airSegsId] == 0.).all()
                except:
                    print("errors in airSegsId",airSegsId)
                    print('errorsEachC',errorsEachC[airSegsId] )
                    print('seg_sol_fluxes',seg_sol_fluxes[airSegsId], proposed_outer_sol_fluxes[airSegsId],
                          seg_mucil_fluxes[airSegsId], proposed_outer_mucil_fluxes[airSegsId])
                    raise Exception
            r.rhizoMassCError_abs  = sum(abs(errorsEachC[rhizoSegsId]))
            if rhizoTotCAfter != 0:
                r.rhizoMassCError_rel = abs(r.rhizoMassCError_abs/rhizoTotCAfter*100)
            else:
                r.rhizoMassCError_rel = np.nan
            
            # instantaneous mass water error for 1DS (with seg_fluxes_limited)
            errorsEachW = rhizoWAfter_ - ( rhizoWBefore_ + (seg_fluxes_limited + seg_fluxes_limited_Out)*dt)
            if len(airSegsId)>0:
                errorsEachW[airSegsId] = 0 # error evaluation only adapted for root segments belowground            
            r.rhizoMassWError_absLim = sum(abs(errorsEachW[rhizoSegsId]))
            r.rhizoMassWError_relLim = abs(r.rhizoMassWError_absLim/sum(rhizoWAfter_[rhizoSegsId])*100)

            
                
            # mass water error to check when leaving fixed point iteration (with seg_fluxes)
            errorsEachW = rhizoWAfter_ - ( rhizoWBefore_ + (seg_fluxes + proposed_outer_fluxes)*dt)
            if len(airSegsId)>0:
                errorsEachW[airSegsId] = 0 # error evaluation only adapted for root segments belowground            
            r.rhizoMassWError_abs = sum(abs(errorsEachW[rhizoSegsId]))
            r.rhizoMassWError_rel = abs(r.rhizoMassWError_abs/sum(rhizoWAfter_[rhizoSegsId])*100)
            if r.mpiVerbose or (max_rank == 1):
                print(rank, "for proposed flux: rhizoMassWError_abs, rel", r.rhizoMassWError_abs,r.rhizoMassWError_rel, 
                    max(abs(seg_fluxes-seg_fluxes_limited)))
                
            
            if False:#r.rhizoMassCError_absLim > 1e-13:
                print('r.rhizoMassCError_absLim > 1e-13',r.rhizoMassCError_absLim,r.rhizoMassWError_absLim )
                print('rhizoTotCAfter_',sum(rhizoTotCAfter_), sum( rhizoTotCBefore_),sum(rhizoTotCAfter_)- sum( rhizoTotCBefore_),
                     'seg_fluxes_limited_sol_In',sum( seg_fluxes_limited_sol_In), 'seg_fluxes_limited_mucil_In',sum(seg_fluxes_limited_mucil_In),
                      'seg_fluxes_limited_sol_Out', sum(seg_fluxes_limited_sol_Out),'seg_fluxes_limited_mucil_Out',sum(seg_fluxes_limited_mucil_Out),
                     'proposed','seg_sol_fluxes', sum(seg_sol_fluxes),'proposed_outer_sol_fluxes',sum(proposed_outer_sol_fluxes), 
                      'seg_mucil_fluxes', sum( seg_mucil_fluxes),'proposed_outer_mucil_fluxes',sum( proposed_outer_mucil_fluxes),' Q_Exud,Q_mucil',sum(Q_Exud),sum(Q_mucil),
                     'diff',sum(seg_sol_fluxes-seg_fluxes_limited_sol_In),sum(proposed_outer_sol_fluxes-seg_fluxes_limited_sol_Out),
                      sum(seg_mucil_fluxes-seg_fluxes_limited_mucil_In),sum(proposed_outer_mucil_fluxes-seg_fluxes_limited_mucil_Out))
                print('seg_fluxes , proposed_outer_fluxes',sum(seg_fluxes ), sum(proposed_outer_fluxes),
                      'seg_fluxes_limited, seg_fluxes_limited_Out',sum(seg_fluxes_limited),sum( seg_fluxes_limited_Out),
                     'diff',sum(seg_fluxes-seg_fluxes_limited),sum(seg_fluxes_limited_Out-proposed_outer_fluxes))
                
                         
            
            ##
            # 2.6 calculate initial 3DS net sources
            ##
            
            # mol per voxel, TODO: find a way to use one function for that and rhizoTotCAfter_
            r.new_soil_solute = np.array( [np.array(r.getC_rhizo(len(cell_volumes), idComp = idc + 1, konz = False)) for idc in range(r.numComp)])
            try:
                assert min(r.new_soil_solute.flatten()) >=0
            except:
                print("min(r.new_soil_solute)",min(r.new_soil_solute.flatten()),[min(nss) for nss in r.new_soil_solute])
                # raise Exception
            
            
            #mol/day
            soil_source_sol = np.full(r.new_soil_solute.shape,0. )
            for nc in range(r.numComp):
                soil_source_sol[nc][cellIds] = np.array(r.new_soil_solute[nc][cellIds] - soil_solute[nc][cellIds] - outer_R_bc_sol[nc][cellIds])/dt

            if r.mpiVerbose or (max_rank == 1):
                comm.barrier()
                print("cyl3plant:soil_source_sol", rank)
                comm.barrier()
            soil_source_sol = comm.bcast(soil_source_sol, root = 0)    
            
            if r.mpiVerbose or (max_rank == 1):
                comm.barrier()
                print("cyl3plant:GOTsoil_source_sol", rank)
                comm.barrier()
                
            assert soil_source_sol.shape == (r.numComp, len(cell_volumes))

            # mabe check here that sum(soil_source_sol) == Qmucil + Qexud
            s.errSoil_source_sol_abs = sum(soil_source_sol.flatten()) - (sum(Q_Exud) + sum(Q_mucil))/dt
            if (sum(Q_Exud) + sum(Q_mucil))/dt != 0.:
                s.errSoil_source_sol_rel = abs(s.errSoil_source_sol_abs/((sum(Q_Exud) + sum(Q_mucil))/dt)*100)
            else:
                s.errSoil_source_sol_rel = np.nan
                
            """ 3. global soil models (3DS)"""
            
            ##
            # 3.1 data before, for post proccessing AND source adaptation
            ##
            if (n_iter > 0) :
                s.reset() #reset at the last moment: over functions use the solution/content at the end of the time step
            solution0_3ds_new = np.array(s.getSolution(0))
            assert (solution0_3ds_new == solution0_3ds_old).all()
            
            if r.mpiVerbose or (max_rank == 1):
                comm.barrier()
                print("cyl3plant:water_contentsoil_water", rank)
                comm.barrier()
            water_content = comm.bcast( np.array(s.getWaterContent()),root= 0)  # theta per cell [1]
            soil_water = np.multiply(water_content, cell_volumes)  # water per cell [cm3]
            buTotCBefore = comm.bcast(sum(s.getTotCContent()), root = 0) 
            
            soil_solute_content = comm.bcast(np.array([np.array(
                s.getContent(i+1, isDissolved = (i < r.numFluidComp))) for i in range(r.numComp)]),
                                                root = 0) # mol

            
            if r.mpiVerbose or (max_rank == 1):
                comm.barrier()
                print("cyl3plant:GOTwater_contentsoil_water", rank)
                comm.barrier()
            ##
            # 3.2 adapt and set sources
            # TODO: take into account BC != 0 (rain, evaporation)
            # TODO: currently, only limit water flow. not clear how to limit solute flow
            ##          
            
            if (rank == 0):
                soil_fluxes_ = rs.sumSegFluxes(seg_fluxes)  # [cm3/day]  per soil cell
                soil_fluxes = np.zeros(len(cell_volumes))
                #easier to handle array than dict. maybe change sumSegFluxes() function to choose type of output
                soil_fluxes[np.array(list(soil_fluxes_.keys()))] = np.array(list(soil_fluxes_.values())) 
            else:
                soil_fluxes = None
            if r.mpiVerbose or (max_rank == 1):
                comm.barrier()
                print("cyl3plant:soil_fluxes", rank)
                comm.barrier()
            soil_fluxes = comm.bcast(soil_fluxes, root=0)
            if r.mpiVerbose or (max_rank == 1):
                comm.barrier()
                print("cyl3plant:GOTsoil_fluxes", rank)
                comm.barrier()
            if (rank == 0):
                soil_fluxes_limited_ = rs.sumSegFluxes(seg_fluxes_limited)  # [cm3/day]  per soil cell
                soil_fluxes_limited = np.zeros(len(cell_volumes))
                #easier to handle array than dict. maybe change sumSegFluxes() function to choose type of output
                soil_fluxes_limited[np.array(list(soil_fluxes_limited_.keys()))] = np.array(list(soil_fluxes_limited_.values()))
            else:
                soil_fluxes_limited = None
            if r.mpiVerbose or (max_rank == 1):
                comm.barrier()
                print("cyl3plant:soil_fluxes_limited", rank)
                comm.barrier()
            soil_fluxes_limited = comm.bcast(soil_fluxes_limited, root=0)
            if r.mpiVerbose or (max_rank == 1):
                comm.barrier()
                print("cyl3plant:GOTsoil_fluxes_limited", rank)
                comm.barrier()
            
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
                    write_file_float("setsource_"+str(idComp), res, directory_ =results_dir) 
                    s.setSource(res.copy(), eq_idx = idComp)  # [mol/day], in modules/richards.py
            r.SinkLim3DS = abs(soil_fluxes_limited - soil_fluxes) # at the end of the fixed point iteration, should be ~ 0        
            
            ##
            # 3.3 solve 3DS
            ##    
            
            start_time_3ds = timeit.default_timer()
            if r.mpiVerbose or (max_rank == 1):
                print("solve 3d soil", rank)
            k_soil_solve = 0
            redoSolve = True
            maxRelShift = s.MaxRelativeShift
            while redoSolve:
                s.ddt = 1.e-5
                try:
                    
                    if r.mpiVerbose or (max_rank == 1):
                        comm.barrier()
                        print("entering the s.solve", rank)
                    s.solve(dt, maxDt = 250/(3600*24))  # in modules/solverbase.py
                    if r.mpiVerbose or (max_rank == 1):
                        print("leaving the s.solve", rank)
                        comm.barrier()
                    redoSolve = False
                    # newton parameters are re-read at each 'solve()' calls
                    s.setParameter("Newton.MaxRelativeShift", str(s.MaxRelativeShift))# reset value
                except Exception as err:
                    if r.mpiVerbose or (max_rank == 1):
                        print(rank, f"Unexpected {err=}, {type(err)=}")
                    if k_soil_solve < 5:
                        if r.mpiVerbose or (max_rank == 1):
                            print(rank,
                              'soil.solve() failed. Increase NewtonMaxRelativeShift from',
                              maxRelShift,'to',maxRelShift*10.)
                        maxRelShift *= 10.
                        # newton parameters are re-read at each 'solve()' calls
                        s.setParameter("Newton.MaxRelativeShift", str(maxRelShift))
                        s.reset()
                        k_soil_solve += 1
                    else:
                        raise Exception
                
            if r.mpiVerbose or (max_rank == 1):
                print("done")
            rs.time_3ds_i += (timeit.default_timer() - start_time_3ds)
            
            ##
            # 3.4 data after, for post proccessing 
            ##
            if r.mpiVerbose or (max_rank == 1):
                comm.barrier()
                print("cyl3plant:buTotCAfter", rank)
                comm.barrier()
            buTotCAfter =comm.bcast( sum(s.getTotCContent()) , root = 0)    
            water_content =comm.bcast( np.array(s.getWaterContent()), root = 0) 
            new_soil_water = np.multiply(water_content, cell_volumes)  # calculate net flux
            
            if r.mpiVerbose or (max_rank == 1):
                comm.barrier()
                print("cyl3plant:GOTbuTotCAfter", rank)
                comm.barrier()
            
            ##
            # 3.5 error rates 
            ##
            s.bulkMassErrorWater_absLim = abs(sum(new_soil_water) - (sum(soil_water) + sum(soil_fluxes_limited)*dt))
            s.bulkMassErrorWater_abs = abs(sum(new_soil_water) - (sum(soil_water) + sum(soil_fluxes)*dt))
            s.bulkMassErrorWater_rel = abs(s.bulkMassErrorWater_abs /sum(new_soil_water) )*100
            
            s.bulkMassCErrorPlant_absReal = buTotCAfter - ( buTotCBefore + sum(Q_Exud) + sum(Q_mucil))
            s.bulkMassCErrorPlant_abs = abs(buTotCAfter - ( buTotCBefore + sum(Q_Exud) + sum(Q_mucil)))
            if buTotCAfter > 0:
                s.bulkMassCErrorPlant_rel = abs(s.bulkMassCErrorPlant_abs/buTotCAfter*100)
            else:
                s.bulkMassCErrorPlant_rel = np.nan
            s.bulkMassCError1ds_abs = abs(buTotCAfter - ( buTotCBefore + sum(soil_source_sol.flatten())*dt))
            if buTotCAfter > 0:
                s.bulkMassCError1ds_rel = abs(s.bulkMassCError1ds_abs/buTotCAfter*100)
            else:
                s.bulkMassCError1ds_rel = np.nan
            
            if rank == 0:
                print("errorCbalance soil 3d?",rank, buTotCAfter ,',', buTotCBefore ,',',  sum(Q_Exud) ,',',  sum(Q_mucil), 
                        ', soil_source_sol', sum(soil_source_sol.flatten())*dt,', s.bulkMassCErrorPlant_abs', s.bulkMassCErrorPlant_abs,
                        ', bulkMassCError1ds_abs ',  s.bulkMassCError1ds_abs , rs_age_i_dt)
            s.buTotCAfter = buTotCAfter
            s.buTotCBefore = buTotCBefore
            
            if rank == 0:
                print("errorWbalance soil 3d?",rank, sum(new_soil_water) ,',', sum(soil_water) ,',',   sum(soil_fluxes)*dt,
                            ', bulkMassErrorWater_abs', s.bulkMassErrorWater_abs,', bulkMassErrorWater_absLim',  s.bulkMassErrorWater_absLim , rs_age_i_dt)

            ##
            # 3.6 get 1DS outer BC (from limited-flux: used at next iteration)
            ##
            outer_R_bc_wat = new_soil_water - soil_water - soil_fluxes_limited *dt # change in water per cell [cm3]
            if r.mpiVerbose or (max_rank == 1):
                comm.barrier()
                print("cyl3plant:soil_solute_content_new", rank)
                comm.barrier()
            soil_solute_content_new = comm.bcast(np.array([np.array(s.getContent(i+1, isDissolved = (i < r.numFluidComp))) for i in range(r.numComp)]), root = 0) # mol
            if r.mpiVerbose or (max_rank == 1):
                comm.barrier()
                print("cyl3plant:GOTsoil_solute_content_new", rank)
                comm.barrier()
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
                                    seg_fluxes =seg_fluxes*0, diff1d3dCW_abs_lim = np.Inf) # just to get error value, will not throw an error
            rx = np.array(rs.psiXyl)
            
            rsx = r.get_inner_heads(weather=r.weatherX)  # matric potential at the segment-exterior interface, i.e. inner values of the (air or soil) cylindric models (not extrapolation to the interface!) [cm]
            if rank == 0:
                rsx_divide = np.where(rsx!=0.,rsx,1.)
                errWrsiAll = abs((rsx - rsx_old)/rsx_divide)*100.
            else:
                errWrsiAll = None
            errWrsiAll = comm.bcast(errWrsiAll, root =0)
            try:
                errWrsi = max(errWrsiAll)#np.linalg.norm(errWrsiAll)
            except:
                print('errWrsi issue',rsx , rsx_old,rsx_divide)
                raise Exception
            
                
            if adaptRSI:
                rsx_set = (rsx+rsx_old)/2
            else:
                rsx_set = rsx
            rsx_old = rsx_set.copy()
            
            
            diffBCS1dsFluxIn =   seg_fluxes  - seg_fluxes_old   #only for water as plant exud is outside of loop
            
            
            
            seg_fluxes_old = seg_fluxes.copy()
            diffBCS1dsFluxOut =   proposed_outer_fluxes  - proposed_outer_fluxes_old 
            diffBCS1dsFluxOut_sol =   proposed_outer_sol_fluxes  - proposed_outer_sol_fluxes_old 
            diffBCS1dsFluxOut_mucil =   proposed_outer_mucil_fluxes  - proposed_outer_mucil_fluxes_old 
            
            proposed_outer_fluxes_old = proposed_outer_fluxes.copy()
            proposed_outer_sol_fluxes_old =proposed_outer_sol_fluxes.copy()
            proposed_outer_mucil_fluxes_old =proposed_outer_mucil_fluxes.copy()
            
            diffouter_R_bc_wat =   outer_R_bc_wat  - outer_R_bc_wat_old 
            outer_R_bc_wat_old = outer_R_bc_wat.copy()
            
            diffouter_R_bc_sol =   outer_R_bc_sol  - outer_R_bc_sol_old 
            outer_R_bc_sol_old = outer_R_bc_sol.copy()
            
            rx_divide = np.where(rx !=0, rx, 1.)
            errRxPlant = max(abs((rx - rx_old)/rx_divide)*100.)#np.linalg.norm(abs((rx - rx_old)/rx_divide)*100)
            rx_old = rx.copy()
            errW1ds = np.linalg.norm(rhizoWAfter_[rhizoSegsId] - rhizoWAfter_old[rhizoSegsId])
            errC1ds = np.linalg.norm(rhizoTotCAfter_[rhizoSegsId] - rhizoTotCAfter_old[rhizoSegsId])
            
            rhizoWAfter_old = rhizoWAfter_.copy()
            rhizoTotCAfter_old = rhizoTotCAfter_.copy()
            
            errBCS1dsFluxIn =max(abs((diffBCS1dsFluxIn/ np.where(seg_fluxes,seg_fluxes,1.))*100))# np.linalg.norm(diffBCS1dsFluxIn)
            errBCS1dsFluxOut = max(abs((diffBCS1dsFluxOut/ np.where(proposed_outer_fluxes,proposed_outer_fluxes,1.))*100))#np.linalg.norm(diffBCS1dsFluxOut)
            errBCS1dsFluxOut_sol = max(abs((diffBCS1dsFluxOut_sol/ np.where(proposed_outer_sol_fluxes,proposed_outer_sol_fluxes,1.))*100))
            errBCS1dsFluxOut_mucil = max(abs((diffBCS1dsFluxOut_mucil/ np.where(proposed_outer_mucil_fluxes,proposed_outer_mucil_fluxes,1.))*100))
            errOuter_R_bc_wat = max(abs((diffouter_R_bc_wat/ np.where(outer_R_bc_wat,outer_R_bc_wat,1.))*100))
            #errOuter_R_bc_sol = max(abs((diffouter_R_bc_sol[:r.numFluidComp].reshape(-1)/ np.where(outer_R_bc_sol[:r.numFluidComp].reshape(-1),
            #                                                                                       outer_R_bc_sol[:r.numFluidComp].reshape(-1),1.))*100))
            #mucilage flow is basically nill
            errOuter_R_bc_sol = max(abs((diffouter_R_bc_sol[:1].reshape(-1)/ np.where(outer_R_bc_sol[:1].reshape(-1),
                                                                                                   outer_R_bc_sol[:1].reshape(-1),1.))*100))
            r.errDiffBCs = np.array([errBCS1dsFluxIn,
                                  errBCS1dsFluxOut, 
                                  errOuter_R_bc_wat,
                                  errBCS1dsFluxOut_sol, 
                                  #errBCS1dsFluxOut_mucil,#mucilage flow is basically nill
                                  errOuter_R_bc_sol
                                    ])
            #sum(abs(diffBCS1dsFluxIn)), sum(abs(diffBCS1dsFluxOut)),sum(abs(diffouter_R_bc_wat)),
            #sum(abs(diffBCS1dsFluxOut_sol.reshape(-1))), sum(abs(diffBCS1dsFluxOut_mucil)),sum(abs(diffouter_R_bc_sol.reshape(-1)))])
            if r.mpiVerbose or (max_rank == 1):
                comm.barrier()
                print("cyl3plant:errDiffBCs", rank)
                comm.barrier()
            r.errDiffBCs = comm.bcast(r.errDiffBCs,root= 0)
            if r.mpiVerbose or (max_rank == 1):
                comm.barrier()
                print("cyl3plant:GOTerrDiffBCs", rank)
                comm.barrier()
            errW3ds = np.linalg.norm(new_soil_water - new_soil_water_old)
            errC3ds = np.linalg.norm(soil_solute_content_new - soil_solute_content_new_old)
            
            new_soil_water_old = new_soil_water.copy()
            soil_solute_content_new_old = soil_solute_content_new.copy()
            
            if r.mpiVerbose or (max_rank == 1):
                comm.barrier()
                print("cyl3plant:errsandCo", rank)
                comm.barrier()
            r.rhizoMassWError_abs = comm.bcast(r.rhizoMassWError_abs,root= 0)
            r.rhizoMassCError_abs = comm.bcast(r.rhizoMassCError_abs,root= 0)
            
            r.err = comm.bcast(max(errRxPlant,r.rhizoMassWError_rel, errWrsi, errW3ds,errC1ds, errC3ds, 
                                   s.bulkMassCErrorPlant_rel, s.bulkMassCError1ds_rel),root= 0)
            r.maxDiff1d3dCW_abs = comm.bcast(r.maxDiff1d3dCW_abs,root= 0)
            r.diff1d3dCurrant = max(r.sumDiff1d3dCW_abs - r.sumDiff1d3dCW_absBU) # to not depend on cumulative error
            r.diff1d3dCurrant_rel =max(r.sumDiff1d3dCW_rel - r.sumDiff1d3dCW_relBU) # to not depend on cumulative error
            if r.mpiVerbose or (max_rank == 1):
                comm.barrier()
                print("cyl3plant:solve_gave_up", rank)
                comm.barrier()
            r.solve_gave_up = (np.array(comm.bcast(comm.gather(r.solve_gave_up ,root = 0),root = 0))).any()
            if r.mpiVerbose or (max_rank == 1):
                comm.barrier()
                print("cyl3plant:solve_gave_up", rank)
                comm.barrier()
            
            if r.mpiVerbose or (max_rank == 1):
                comm.barrier()
                print("cyl3plant:GOTerrsandCo", rank)
                comm.barrier()
            r.errs =np.array([errRxPlant, errW1ds, errW3ds,errC1ds, errC3ds, 
                            max(r.SinkLim3DS),max(abs(r.SinkLim1DS)),max(abs(r.OutLim1DS)),
                            max(abs(r.InOutBC_Cdiff.reshape(-1))),
                            max(r.maxDiff1d3dCW_abs), 
                            errWrsi,#maxDiff1d3dCW_absBU, 
                            s.bulkMassErrorWater_abs,s.bulkMassErrorWater_absLim,
                            r.rhizoMassWError_absLim,r.rhizoMassWError_abs,
                            s.bulkMassCError1ds_abs,s.bulkMassCErrorPlant_abs,
                            r.rhizoMassCError_absLim,r.rhizoMassCError_abs,
                            sum(abs(diffBCS1dsFluxIn)), sum(abs(diffBCS1dsFluxOut)),sum(abs(diffouter_R_bc_wat)),
                            sum(abs(diffBCS1dsFluxOut_sol.reshape(-1))),sum(abs(diffBCS1dsFluxOut_mucil)),
                              sum(abs(diffouter_R_bc_sol.reshape(-1))), 
                           r.diff1d3dCurrant,r.diff1d3dCurrant_rel, r.rhizoMassWError_rel,r.err ])
            
            rhizoWaterPerVoxel = r.getWaterVolumesCyl(doSum = False, reOrder = True)
            
            theta3ds = s.getWaterContent()# proposed_outer_mucil_fluxes
            if  (n_iter % skip == 0) and (rank == 0):
                write_file_array("fpit_sumDiff1d3dCW_relBU", r.sumDiff1d3dCW_relBU, directory_ =results_dir, fileType = '.csv') 
                write_file_array("fpit_sumDiff1d3dCW_absBU", r.sumDiff1d3dCW_absBU, directory_ =results_dir, fileType = '.csv') 
                write_file_array("fpit_errbulkMass",np.array([s.bulkMassCErrorPlant_abs, s.bulkMassCErrorPlant_rel, #not cumulative 
                                                s.bulkMassCError1ds_abs, s.bulkMassCError1ds_rel, 
                                                s.bulkMassErrorWater_abs,s.bulkMassErrorWater_rel,
                                                s.bulkMassCErrorPlant_absReal,
                                                s.errSoil_source_sol_abs, s.errSoil_source_sol_rel]), directory_ =results_dir, fileType = '.csv') 
                write_file_array("fpit_errorMassRhizo", np.array([r.rhizoMassCError_abs, r.rhizoMassCError_rel,
                                                        r.rhizoMassWError_abs, r.rhizoMassWError_rel]), 
                                 directory_ =results_dir, fileType = '.csv')# not cumulativecumulative (?)
                write_file_array("fpit_errDiffBCs", r.errDiffBCs, directory_ =results_dir, fileType = '.csv') 
                write_file_array("fpit_diffBCS1dsFluxOut_sol", diffBCS1dsFluxOut, directory_ =results_dir, fileType = '.csv') 
                write_file_array("fpit_diffBCS1dsFluxOut_mucil", diffBCS1dsFluxOut, directory_ =results_dir, fileType = '.csv') 
                
                write_file_array("fpit_errBCS1dsFluxOut_mucilBis", abs((diffBCS1dsFluxOut_mucil/ np.where(proposed_outer_mucil_fluxes,proposed_outer_mucil_fluxes,1.))*100), directory_ =results_dir, fileType = '.csv') 
                
                write_file_array("fpit_diffBCS1dsFluxInBis", abs((diffBCS1dsFluxIn/ np.where(seg_fluxes,seg_fluxes,1.))*100), directory_ =results_dir, fileType = '.csv') 
                write_file_array("fpit_diffBCS1dsFluxOutBis", abs((diffBCS1dsFluxOut/ np.where(proposed_outer_fluxes,proposed_outer_fluxes,1.))*100), directory_ =results_dir, fileType = '.csv') 
                write_file_array("fpit_diffouter_R_bc_watBis", abs((diffBCS1dsFluxOut_sol/ np.where(proposed_outer_sol_fluxes,proposed_outer_sol_fluxes,1.))*100), directory_ =results_dir, fileType = '.csv') 
                
                write_file_array("fpit_diffBCS1dsFluxIn", diffBCS1dsFluxIn, directory_ =results_dir, fileType = '.csv') 
                write_file_array("fpit_diffBCS1dsFluxOut", diffBCS1dsFluxOut, directory_ =results_dir, fileType = '.csv') 
                write_file_array("fpit_diffouter_R_bc_wat", diffouter_R_bc_wat, directory_ =results_dir, fileType = '.csv') 
                write_file_array("fpit_diffBCS1dsFluxIn_rhizoSegsId", diffBCS1dsFluxIn[rhizoSegsId], directory_ =results_dir, fileType = '.csv') 
                write_file_array("fpit_diffBCS1dsFluxOut_rhizoSegsId", diffBCS1dsFluxOut[rhizoSegsId], directory_ =results_dir, fileType = '.csv') 
                
                
                if (plantType != "plant") :
                    write_file_array('fpit_transRate',np.array([transpiration]), directory_ =results_dir, fileType = '.csv' )
                    write_file_array("fpit_n_iter",np.array([ n_iter, r.solve_gave_up ]), directory_ =results_dir, fileType = '.csv') 

                write_file_array('fpit_watVolTheta',np.array([sum(new_soil_water),np.mean(theta3ds)]), directory_ =results_dir, 
                                 fileType = '.csv' )
                write_file_array("fpit_errorMassW1d", np.array(errorsEachW), directory_ =results_dir, fileType = '.csv') 
                write_file_array("fpit_errorMassC1d", np.array(errorsEachC), directory_ =results_dir, fileType = '.csv') 
                write_file_array("fpit_error", r.errs, directory_ =results_dir, fileType = '.csv') 
                write_file_array("fpit_error1d3d", r.maxDiff1d3dCW_abs, directory_ =results_dir, fileType = '.csv') 
                write_file_array("fpit_time", np.array([rs_age,rs.Qlight,sim_time,dt,Ni,n_iter]), directory_ =results_dir ) 
                write_file_array("fpit_SinkLim3DS", r.SinkLim3DS, directory_ =results_dir, fileType = '.csv') 
                write_file_array("fpit_SinkLim1DS", r.SinkLim1DS, directory_ =results_dir, fileType = '.csv') 
                
                
                write_file_array("fpit_sol_Out_diff", seg_fluxes_limited_sol_Out-proposed_outer_sol_fluxes, directory_ =results_dir, 
                                 fileType = '.csv') 
                write_file_array("fpit_mucil_Out_diff",  seg_fluxes_limited_mucil_Out-proposed_outer_mucil_fluxes, directory_ =results_dir, 
                                 fileType = '.csv') 
                write_file_array("fpit_sol_In_diff", seg_fluxes_limited_sol_In-seg_sol_fluxes, directory_ =results_dir, 
                                 fileType = '.csv') 
                write_file_array("fpit_mucil_In_diff",  seg_fluxes_limited_mucil_In-seg_mucil_fluxes, directory_ =results_dir, 
                                 fileType = '.csv') 
                write_file_array("fpit_InOutBC_Cdiff", r.InOutBC_Cdiff.reshape(-1), directory_ =results_dir, 
                                 fileType = '.csv') 
                
                write_file_array("fpit_seg_fluxes_limited", seg_fluxes_limited, directory_ =results_dir, fileType = '.csv') 
                write_file_array("fpit_seg_fluxes", seg_fluxes, directory_ =results_dir, fileType = '.csv') 
                write_file_array("fpit_seg_fluxes_limited_Out",seg_fluxes_limited_Out, directory_ =results_dir, fileType = '.csv')
                write_file_array("fpit_proposed_outer_fluxes", proposed_outer_fluxes, directory_ =results_dir, fileType = '.csv')
                
                # solutes. only limited vs unlimited for the rhizo: not relevent for soil source
                # not sure i need to have it for the inner BC as plant suc flow outside of simloop and we should only have 
                # exud, never uptake.
                #
                write_file_array("fpit_seg_fluxes_limited_sol_Out", seg_fluxes_limited_sol_Out, directory_ =results_dir, fileType = '.csv') 
                write_file_array("fpit_seg_fluxes_sol_Out", proposed_outer_sol_fluxes, directory_ =results_dir, fileType = '.csv') 
                
                write_file_array("fpit_seg_fluxes_limited_mucil_Out", seg_fluxes_limited_mucil_Out, directory_ =results_dir, fileType = '.csv') 
                write_file_array("fpit_seg_fluxes_mucil_Out", proposed_outer_mucil_fluxes, directory_ =results_dir, fileType = '.csv') 
                
                write_file_array("fpit_seg_fluxes_limited_sol_In", seg_fluxes_limited_sol_In, directory_ =results_dir, fileType = '.csv') 
                write_file_array("fpit_seg_fluxes_sol_In", seg_sol_fluxes, directory_ =results_dir, fileType = '.csv') 
                
                write_file_array("fpit_seg_fluxes_limited_mucil_In", seg_fluxes_limited_mucil_In, directory_ =results_dir, fileType = '.csv') 
                write_file_array("fpit_seg_fluxes_mucil_In", seg_mucil_fluxes, directory_ =results_dir, fileType = '.csv') 
                
                
                write_file_array("fpit_soil_fluxes_limited", soil_fluxes_limited, directory_ =results_dir, fileType = '.csv') 
                write_file_array("fpit_soil_fluxes", soil_fluxes, directory_ =results_dir, fileType = '.csv') 
                
                for nc in range(r.numComp):# normally all 0 for nc >= numFluidComp
                    write_file_array("fpit_outer_R_bc_sol"+str(nc+1), outer_R_bc_sol[nc], directory_ =results_dir, fileType = '.csv')  
                    write_file_array("fpit_diffouter_R_bc_sol"+str(nc+1), diffouter_R_bc_sol[nc], directory_ =results_dir, fileType = '.csv') 
                    
                    write_file_array("fpit_diffouter_R_bc_solBis"+str(nc+1), abs((diffouter_R_bc_sol[nc]/ np.where(outer_R_bc_sol[nc],outer_R_bc_sol[nc],1.))*100), directory_ =results_dir, fileType = '.csv') 
                    write_file_float("fpit_diffouter_R_bc_solBismax"+str(nc+1),max( abs((diffouter_R_bc_sol[nc]/ np.where(outer_R_bc_sol[nc],outer_R_bc_sol[nc],1.))*100)), directory_ =results_dir) 
                
                write_file_array("fpit_outer_R_bc_wat", outer_R_bc_wat, directory_ =results_dir, fileType = '.csv')  
                
                write_file_array("fpit_all1d3dDiff",r.all1d3dDiff, directory_ =results_dir, fileType = '.csv') 
                write_file_array("fpit_all1d3dDiffBis",r.all1d3dDiff[cellIds], directory_ =results_dir, fileType = '.csv') 
                
                write_file_array("fpit_outer_R_bc_watBis", outer_R_bc_wat[cellIds], directory_ =results_dir) 
                write_file_array("fpit_psi_sri_error", errWrsiAll, directory_ =results_dir, fileType = '.csv') 
                write_file_array("fpit_psi_sri_errorRoot", errWrsiAll[rhizoSegsId], directory_ =results_dir, fileType = '.csv') 
                write_file_array("fpit_psi_sri_set", rsx_set, directory_ =results_dir, fileType = '.csv') 
                write_file_array("fpit_psi_sri_real", rsx, directory_ =results_dir, fileType = '.csv') 
                write_file_array("fpit_new_soil_water", new_soil_water, directory_ =results_dir, fileType = '.csv') 
                write_file_array("fpit_rhizoWaterPerVoxel", rhizoWaterPerVoxel, directory_ =results_dir, fileType = '.csv') 
                write_file_array("fpit_rx", rx, directory_ =results_dir, fileType = '.csv') 
                write_file_array("fpit_sumErrors1ds3ds", np.concatenate((r.sumDiff1d3dCW_abs, r.sumDiff1d3dCW_rel)), directory_ =results_dir, fileType = '.csv')
                write_file_array("fpit_maxErrors1ds3ds", np.concatenate((r.maxDiff1d3dCW_abs, r.maxDiff1d3dCW_rel)), directory_ =results_dir, fileType = '.csv')# 
                
                write_file_array("fpit_rhizoWAfter", rhizoWAfter_[rhizoSegsId] , directory_ =results_dir, fileType = '.csv') 
                write_file_array("fpit_rhizoTotCAfter", rhizoTotCAfter_[rhizoSegsId] , directory_ =results_dir, fileType = '.csv') 

            if r.mpiVerbose or (max_rank == 1):
                comm.barrier()
                print('end iteration', rank, n_iter, r.err,r.maxDiff1d3dCW_abs)
                comm.barrier()
            n_iter += 1
            #end iteration
        if rank == 0:
            write_file_array("N_seg_fluxes",np.array(rs.outputFlux), directory_ =results_dir)
            write_file_array("N_soil_fluxes",soil_fluxes, directory_ =results_dir)
            write_file_array("N_error", r.errs, directory_ =results_dir) 
            write_file_array("N_error", r.errs, directory_ =results_dir, fileType = '.csv') 
            write_file_array("N_n_iter",np.array([ n_iter,N ]), directory_ =results_dir)
        
            
        ####
        #   error rates    
        ####
        

        # error 3DS-1DS
        if r.mpiVerbose or (max_rank == 1):
            comm.barrier()
            print('error 3DS-1DS', rank)
        r.checkMassOMoleBalance2(soil_fluxes*0, soil_source_sol*0, dt,
                                seg_fluxes =seg_fluxes*0, diff1d3dCW_abs_lim = np.Inf)
        if r.mpiVerbose or (max_rank == 1):
            print('finished  error 3DS-1DS', rank)
            comm.barrier()
                                
        if rank == 0:            
            for nc in range(r.numFluidComp, r.numComp):
                try:
                    assert (abs(outer_R_bc_sol[nc][cellIds]) < 1e-16).all()
                except:
                    print("outer_R_bc_sol[nc][cellIds] != 0.",rank, nc+1, cellIds)
                    print(outer_R_bc_sol[nc][cellIds])
                    print(soil_solute_content_new[nc] , soil_solute_content[nc] , soil_source_sol[nc]*dt)
                    raise Exception
                    
        if r.mpiVerbose or (max_rank == 1):
            comm.barrier()
            print('end time step inner loop')
            comm.barrier()
        
        
        failedLoop = continueLoop(r,0, dt, False,Ni * dt,'fpit_testdata')
        if r.mpiVerbose or (max_rank == 1):
            print('left iteration', rank, n_iter,Ni,'/',N, r.err, max(r.maxDiff1d3dCW_abs), r.rhizoMassWError_abs,'failed?', failedLoop)
            comm.barrier()
        if (failedLoop):# no need to go on, leave inner loop now and reset lower time step
            print('Failed, no need to go on, leave inner loop now and reset lower time step')
            break
        # end time step inner loop
    dt_inner =float(Ni+1)*float( dt) # get the real simetime if sim_time / dt != int
    if r.mpiVerbose or (max_rank == 1):
        print('end of inner loop, failed?',failedLoop, n_iter,Ni,'/',N, dt_inner, dt)
    
    return outer_R_bc_sol, outer_R_bc_wat, seg_fluxes, dt_inner, failedLoop# fluxes == first guess for next fixed point iteration
    #end of inner loop
