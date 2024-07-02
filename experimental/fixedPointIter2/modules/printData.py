
import numpy as np
import pandas as pd
import timeit
import pandas as pd
import matplotlib; matplotlib.use('agg')
import sys;
import os
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank(); max_rank = comm.Get_size()
import timeit
import plantbox as pb

from air_modelsPlant import AirSegment
from helpfull import write_file_array, write_file_float, div0, div0f
import vtk_plot_adapted as vp

def initialPrint(plant):
    """ put headings on files which will be filled later """
    columnNames = np.array([# non-convergence metrics
                "errRxPlant","errW1ds","errW3ds","errC1ds","errC3ds","errWrsi","errWrsiRealInput",
                "errBCS1dsFluxIn","errBCS1dsFluxOut","errBCS1dsFluxOut_sol","errBCS1dsFluxOut_mucil",
                "errOuter_R_bc_wat","errOuter_R_bc_sol",
                # realised vs prescribed fluxes and sinks
                "SinkLim3DS",
                "SinkLim1DS","OutLim1DS",
                "OutBCsol_diff","OutBCmucil_diff","InBCsol_diff","InBC_mucildiff",
                # 1d-3d differences/errors
                "sumDiff1d3dCW_abs","sumDiff1d3dCW_rel","diff1d3dCurrant","diff1d3dCurrant_rel",
                "maxDiff1d3dCW_abs","maxDiff1d3dCW_rel","maxdiff1d3dCurrant","maxdiff1d3dCurrant_rel",
                # mass balance error 3d model
                "bulkMassErrorWater_abs","bulkMassErrorWater_rel",
                "bulkMassErrorWater_absLim","bulkMassErrorWater_relLim",
                "bulkMassCErrorPlant_abs","bulkMassCError1ds_abs",
                "bulkMassCErrorPlant_rel","bulkMassCError1ds_rel",
                # mass balance error 1d models
                "rhizoMassWError_relLim","rhizoMassCError_relLim",
                "rhizoMassWError_rel","rhizoMassCError_rel",
                # summary metric
                "err"])
    write_file_array("inner_error", columnNames, directory_ =plant.results_dir, fileType = '.csv')
    
        
def printPlantShape(rs,r):
    """ store plant shape data after doing growth simulation"""
    results_dir = rs.results_dir
    write_file_array('seg2cell_keys',rs.seg2cell,directory_ =results_dir, 
                             fileType = '.csv')
    write_file_array('seg2cell_vals',np.array(list(rs.seg2cell.values())),
                     directory_ =results_dir, fileType = '.csv')
    write_file_array("organTypes", np.array(rs.organTypes), directory_ =results_dir)
    write_file_array("subTypes", np.array(r.rs.subTypes), directory_ =results_dir)
    length_Segs = np.array(r.rs.segLength())
    write_file_array("length_Segs", length_Segs, directory_ =results_dir)
    write_file_array("root_radii", np.array(r.rs.radii), directory_ =results_dir)
    orgs = r.plant.getOrgans(-1, False)
    parentorgid = np.array([org.getParent().getId() for org in orgs])
    parentNidx = np.array([org.getParent().getNodeId(org.parentNI) for org in orgs])
    id_orgs = np.array([org.getId() for org in orgs])
    ot_orgs = np.array([org.organType() for org in orgs])
    st_orgs = np.array([org.getParameter("subType") for org in orgs])
    lenOrg = np.array([org.getLength(False) for org in orgs])   
    write_file_array("id_orgs", id_orgs, directory_ =results_dir)
    write_file_array("ot_orgs", ot_orgs, directory_ =results_dir)
    write_file_array("st_orgs", st_orgs, directory_ =results_dir)  
    write_file_array("lenOrg", lenOrg, directory_ =results_dir) 
    write_file_array("parentNidx", parentNidx, directory_ =results_dir) 

    write_file_array("nodes_X",
                     np.array([tempnode[0] for tempnode in r.get_nodes()]), 
                     directory_ =results_dir)
    write_file_array("nodes_Y", np.array([tempnode[1] for tempnode in r.get_nodes()]), directory_ =results_dir)
    write_file_array("nodes_Z", np.array([tempnode[2] for tempnode in r.get_nodes()]), directory_ =results_dir)
    idPerNode = np.concatenate([
        np.full(org.getNumberOfNodes()-1,org.getId()) for org in orgs]).reshape(-1)
    globalNodeId = np.concatenate([org.getNodeIds()[1:] for org in orgs]).reshape(-1)
    write_file_array("orgidPerNode", idPerNode, directory_ =results_dir)
    write_file_array("globalNodeId", globalNodeId, directory_ =results_dir)
    volOrg = np.array([org.orgVolume(-1,False) for org in orgs]) 
    write_file_array("volOrg", volOrg, directory_ =results_dir)
    
def printDiff1d3d(perirhizalModel, s):
    """print differences between 1d and 3d soil models"""
    results_dir = perirhizalModel.results_dir
    write_file_array("content_diff1d3dabs"+str(0), 
                         perirhizalModel.allDiff1d3dCW_abs[0], directory_ =results_dir, fileType = '.csv')
    write_file_array("content_diff1d3drel"+str(0), 
                     perirhizalModel.allDiff1d3dCW_rel[0], directory_ =results_dir, fileType = '.csv')
    for nc in range(perirhizalModel.numSoluteComp):# normally all 0 for nc >= numDissolvedSoluteComp
            write_file_array("content_diff1d3dabs"+str(nc+1), 
                             perirhizalModel.allDiff1d3dCW_abs[nc+1], directory_ =results_dir, fileType = '.csv')
            write_file_array("content_diff1d3drel"+str(nc+1), 
                             perirhizalModel.allDiff1d3dCW_rel[nc+1], directory_ =results_dir, fileType = '.csv')

            write_file_array("sol_content3d"+str(nc+1), 
                             s.getContent(nc+1), 
                             directory_ =results_dir, fileType = '.csv')  
            write_file_array("sol_content1d"+str(nc+1), 
                             perirhizalModel.getContentCyl(idComp = nc+1, doSum = False, reOrder = True), 
                             directory_ =results_dir, fileType = '.csv') 
            

    write_file_array("pHead", np.array(s.getSolutionHead()), directory_ =results_dir, fileType = '.csv') 
    write_file_array("theta", np.array(s.getWaterContent()), directory_ =results_dir, fileType = '.csv') 
    write_file_array("watVol", np.array(s.getWaterVolumes()), directory_ =results_dir, fileType = '.csv') 
    for i in range(perirhizalModel.numDissolvedSoluteComp):

        write_file_array("Soil_solute_conc"+str(i+1), 
                         np.array(s.getSolution(i+1)).flatten()* perirhizalModel.molarDensityWat_m3/1e6, # mol C/mol w * molW/m3 W * m3 W/cm3 W
                         directory_ =results_dir, fileType = '.csv') 
    for i in range(perirhizalModel.numDissolvedSoluteComp, perirhizalModel.numSoluteComp):

        write_file_array("Soil_solute_conc"+str(i+1), 
                         np.array(s.getSolution(i+1)).flatten()* perirhizalModel.bulkDensity_m3 /1e6 , 
                         directory_ =results_dir, fileType = '.csv') 

            
def printTimeAndError(rs, rs_age):
    """
        get sum of the (resp. max) error/difference between 1d and 3d model for each element (water + solutes)
    """
    results_dir = rs.results_dir
    write_file_array("time", np.array([rs_age,0.]), directory_ =results_dir)
    write_file_array("sumErrors1ds3ds", np.concatenate((rs.sumDiff1d3dCW_abs, rs.sumDiff1d3dCW_rel)), directory_ =results_dir, fileType = '.csv')
    write_file_array("maxErrors1ds3ds", np.concatenate((rs.maxDiff1d3dCW_abs, rs.maxDiff1d3dCW_rel)), directory_ =results_dir, fileType = '.csv') 
    
    
def printOutput(rs_age, perirhizalModel, phloemDataStorage, plantModel):
    """ print and save outputs of phloem flow and photosynthesis"""
    results_dir = perirhizalModel.results_dir
    
    print("\n\n\n\t\t", int(rs_age//1),"d", int(rs_age%1*24),"h", int(rs_age%1*24%1*60),"mn")
    write_file_float("trans", plantModel.TranspirationCumul, directory_ =results_dir)
    if perirhizalModel.doPhloemFlow:
        print(round(plantModel.Qlight *1e6),"mumol m-2 s-1")
        print("Error in Suc_balance:\n\tabs (mmol) {:5.2e}\trel (-) {:5.2e}".format(phloemDataStorage.error_st_abs,
                                                                    phloemDataStorage.error_st_rel))
        print("Error in photos:\n\tabs (cm3/day) {:5.2e}".format(abs(sum(plantModel.outputFlux))))
        print("C_ST (mol ml-1):\n\tmean {:.2e}\tmin  {:5.2e} at {:d} segs \tmax  {:5.2e}".format(np.mean(phloemDataStorage.C_ST), min(phloemDataStorage.C_ST),
                                                                                                  len(np.where(phloemDataStorage.C_ST == min(phloemDataStorage.C_ST) )[0]), max(phloemDataStorage.C_ST)))        
        print("C_me (mol ml-1):\n\tmean {:.2e}\tmin  {:5.2e}\tmax  {:5.2e}".format(np.mean(phloemDataStorage.C_meso), min(phloemDataStorage.C_meso), max(phloemDataStorage.C_meso)))        
        print('Q_X (mol Suc): \n\tST   {:.2e}\tmeso {:5.2e}\tin   {:5.2e}'.format(sum(phloemDataStorage.Q_ST), sum(phloemDataStorage.Q_meso), phloemDataStorage.Q_in))
        print('\tRm   {:.2e}\tGr   {:.2e}\tExud {:5.2e}'.format(sum(phloemDataStorage.Q_Rm), sum(phloemDataStorage.Q_Gr), sum(phloemDataStorage.Q_Exud)))

        if min(phloemDataStorage.C_ST) < 0.0:
            print("min(C_ST) < 0.0", min(phloemDataStorage.C_ST),np.mean(phloemDataStorage.C_ST),max(phloemDataStorage.C_ST))
            raise Exception

        tiproots, tipstems, tipleaves = plantModel.get_organ_segments_tips()
        write_file_array("root_segments_tips",tiproots, 
                         directory_ =results_dir)
        write_file_array("rhizoSegsId", np.array(perirhizalModel.rhizoSegsId), 
                         directory_ =results_dir)
        write_file_array("Q_ST", phloemDataStorage.Q_ST, directory_ =results_dir)#mmol
        write_file_array("C_ST", phloemDataStorage.C_ST, directory_ =results_dir)#mmol/cm3
        write_file_array("C_meso", phloemDataStorage.C_meso, directory_ =results_dir)
        write_file_array("Q_meso", phloemDataStorage.Q_meso, directory_ =results_dir)


        write_file_array("Q_S_ST", phloemDataStorage.Q_S_ST, directory_ =results_dir)#mmol
        write_file_array("C_S_ST", phloemDataStorage.C_S_ST, directory_ =results_dir)#mmol/cm3
        write_file_array("C_S_meso", phloemDataStorage.C_S_meso, directory_ =results_dir)
        write_file_array("Q_S_meso", phloemDataStorage.Q_S_meso, directory_ =results_dir)

        write_file_array("Q_Rm", phloemDataStorage.Q_Rm, directory_ =results_dir)
        write_file_array("Q_Mucil", phloemDataStorage.Q_Mucil, directory_ =results_dir)
        write_file_array("Q_Exud", phloemDataStorage.Q_Exud, directory_ =results_dir)
        write_file_array("Q_Gr", phloemDataStorage.Q_Gr, directory_ =results_dir)
        
        
        write_file_float("Q_Exud_i", sum(phloemDataStorage.Q_Exud_i_seg), directory_ =results_dir)
        write_file_float("Q_Mucil_i", sum(phloemDataStorage.Q_Mucil_i_seg), directory_ =results_dir)
        write_file_float("Q_Exud_tot", phloemDataStorage.Q_Exud_cumul, directory_ =results_dir)
        write_file_float("Q_Mucil_tot", phloemDataStorage.Q_Mucil_cumul, directory_ =results_dir)
        
        write_file_float("Q_Ag", phloemDataStorage.Q_in, directory_ =results_dir)
        write_file_array("psiXyl", plantModel.psiXyl, directory_ =results_dir)
        write_file_array("Fpsi", plantModel.Fpsi, directory_ =results_dir)
        write_file_array("fw", plantModel.fw, directory_ =results_dir)
        write_file_array("gco2", plantModel.gco2, directory_ =results_dir)
        write_file_array("Q_Ag_dot", plantModel.AgPhl, directory_ =results_dir)
        write_file_array("C_rsi", np.array(plantModel.Csoil_seg ), 
                         directory_ =results_dir)#mmol/cm3
        write_file_array("transrate",plantModel.Jw, directory_ =results_dir, fileType = '.csv')
        
        write_file_array("errorsPlant", np.array([phloemDataStorage.error_st_abs,
                                                  phloemDataStorage.error_st_rel,#cumulative
                                                 abs(sum(plantModel.outputFlux))]), 
                         directory_ =results_dir, fileType = '.csv') #not cumulative

def printCylData(perirhizalModel, rs_age):
    """ print data of each perirhizal model"""
    results_dir = perirhizalModel.results_dir
    
    for lId, cyl in enumerate(perirhizalModel.cyls):
        if not isinstance(cyl, AirSegment):
            gId = perirhizalModel.eidx[lId]

            pHead = np.array(cyl.getSolutionHead()).flatten()
            write_file_float("Cyl_time_"+str(gId),rs_age, 
                                 directory_ =results_dir+'cyl_val/', allranks = True)
            write_file_array("Cyl_watercontent_"+str(gId),cyl.getWaterContent(), 
                             directory_ =results_dir+'cyl_val/', allranks = True)
            write_file_array("Cyl_pressureHead_"+str(gId),pHead, 
                             directory_ =results_dir+'cyl_val/', allranks = True)
            write_file_array("Cyl_coord_"+str(gId),
                         cyl.getDofCoordinates().flatten(), 
                         directory_ =results_dir+'cyl_val/', allranks = True)
            write_file_array("Cyl_cellVol_"+str(gId),
                             cyl.getCellVolumes().flatten(), 
                             directory_ =results_dir+'cyl_val/', allranks = True)
            if perirhizalModel.doSoluteFlow:
                for ccc in range(perirhizalModel.numSoluteComp):
                    sol0 = np.array(cyl.getContent(ccc +1)).flatten()
                    write_file_array("Cyl_content"+str(ccc+1)+"_"+str(gId)+"", 
                                 sol0, 
                                 directory_ =results_dir+'cyl_val/', allranks = True)
                    if min(sol0) < 0.:
                        print('issue sol0',sol0, 'solid',ccc+1)
                        raise Exception
            if max(pHead) > 0:
                print('issue phead',gId,rank, pHead )
                raise Exception
                
def getAndPrintErrorRates(perirhizalModel, plantModel, s, phloemData):
    """
        define and save other error rates
    """
    results_dir = perirhizalModel.results_dir
    write_file_array("endFpitLoop_error", perirhizalModel.errs, 
                     directory_ =results_dir, fileType = '.csv') 
    write_file_array("endFpitLoop_data", 
                     np.array([perirhizalModel.n_iter, perirhizalModel.err, 
                               perirhizalModel.rhizoMassWError_abs,
                               perirhizalModel.dt_inner]), 
                     directory_ =results_dir, fileType = '.csv')

    totC3dAfter = sum(s.getTotCContent())  

    soil_water3dAfter = sum(np.multiply(np.array(s.getWaterContent()), s.getCellVolumes()))    

    if rank == 0:
        if perirhizalModel.doSoluteFlow:
            s.bulkMassErrorCumul_abs = abs((totC3dAfter - ( s.totC3dInit + 
                                        sum(phloemData.Q_Exud) + 
                                        sum(phloemData.Q_Mucil))))
            if totC3dAfter != 0:
                s.bulkMassErrorCumul_rel = abs(s.bulkMassErrorCumul_abs/totC3dAfter*100)
            else:
                s.bulkMassErrorCumul_rel =np.nan
        else:
            s.bulkMassErrorCumul_abs = np.nan
            s.bulkMassErrorCumul_rel =np.nan
        # if we have a free flow BC at the bottom, that could increase the error
        # ideally, I should add the flow at the bellow BC here
        s.bulkMassErrorWaterCumul_abs = abs(soil_water3dAfter - ( s.buWSoilInit - plantModel.TranspirationCumul_eval))
        s.bulkMassErrorWaterCumul_rel = abs(s.bulkMassErrorWaterCumul_abs/soil_water3dAfter*100)
    else:
        s.bulkMassErrorCumul_abs = None
        s.bulkMassErrorCumul_rel = None
        s.bulkMassErrorWaterCumul_abs = None
        s.bulkMassErrorWaterCumul_rel = None

    if perirhizalModel.doSoluteFlow:
        #absolute and relative (%) error
        write_file_array("errorsBulkSoilCumulative", np.array([
                                            s.bulkMassErrorCumul_abs,s.bulkMassErrorCumul_rel,#cumulative
                                            s.bulkMassErrorWaterCumul_abs,s.bulkMassErrorWaterCumul_rel]),
                         directory_ =results_dir, fileType = '.csv')#cumulative

        

def doVTPplots(vtpindx, perirhizalModel, plantModel, s, 
                soilTextureAndShape,datas = [], datasName=[], 
                initPrint=False, doSolutes= False):
    """ vtp plots for water and C components in 3d soil and at root-soil
        interface 
    """
    results_dir = perirhizalModel.results_dir
    if rank == 0:
        vp.plot_plant(plantModel.plant,p_name =datasName,
                            vals =datas, 
                            filename =results_dir+"vtpvti/plantAt"+ str(vtpindx), 
                      range_ = [0,5000])
    if not initPrint:
        cell_volumes = s.getCellVolumes()
        # otherwise the shoot looks weird
        periodic = False
        
        
        extraArray_ = perirhizalModel.soilModel.getWaterContent()
        extraArrayName_ = "theta (cm3/cm3)"
        
        getConcentration = True # True: get concentration, False: get content
        rsi_Conce = perirhizalModel.get_concentrationOrContent(0, getConcentration)
        
        # TODO: adapt to have plot_plants_and_soil (i.e., with leaf surface)
        vp.plot_roots_and_soil(perirhizalModel.mappedSegments(),extraArrayName_,rsi_Conce, s, periodic, 
                               soilTextureAndShape['min_b'],
                               soilTextureAndShape['max_b'],
                               soilTextureAndShape['cell_number'], 
                filename="vtpvti/soil_rx"+ str(vtpindx),sol_ind =-1,
                               extraArray = extraArray_, extraArrayName = extraArrayName_,
                interactiveImage=False)  # VTK vizualisation
        
        if doSolutes:        
            for i in range(1, perirhizalModel.numComp):
                extraArray_ = perirhizalModel.soilModel.getSolution(i) * perirhizalModel.phaseDensity(i)/1e6

                if not getConcentration:
                    extraArrayName_ = "C"+str(i)+" mol"
                else:
                    extraArrayName_ = "[C"+str(i)+"] (mol/cm3)"
                    if i <= perirhizalModel.numDissolvedSoluteComp:
                        extraArray_ /= (perirhizalModel.soilModel.getWaterContent() *cell_volumes) #cm3 water
                    else: 
                        extraArray_ /= cell_volumes #cm3


                vp.plot_roots_and_soil(perirhizalModel.mappedSegments(),
                                       extraArrayName_ ,
                                       perirhizalModel.get_concentrationOrContent(i , getConcentration), s, 
                                       periodic,
                                   soilTextureAndShape['min_b'],
                                   soilTextureAndShape['max_b'],
                                       soilTextureAndShape['cell_number'], 
                        filename="vtpvti/C"+str(i)+'_'+ str(vtpindx), 
                                       sol_ind =-1,
                                       extraArray = extraArray_, 
                        extraArrayName = extraArrayName_,
                    interactiveImage=False)  # VTK vizualisation
    if rank == 0:
        print('did VTP print in file',results_dir+"vtpvti/plantAt"+ str(vtpindx) )
        
        
def errorWeatherChange(results_dir, cyl, pheadOld,nc_content, nc_content_new, nc_molFr):
    print('the solute content changed',rank, cyl.gId,'nc_content',nc_content,nc_content_new,'vols',
                          cyl_cell_volumes,'newWatMol',newWatMol,'nc_molFr',nc_molFr)
    print('evaluation error',
         np.maximum(abs(
        nc_content.reshape(-1) - nc_content_new.reshape(-1)
         ))
         )
    write_file_array('errorChange_phead'+str(cyl.gId),
                     pheadOld, 
                     directory_ =results_dir,  allranks = True)
    write_file_array('errorChange_phead'+str(cyl.gId),
                     cyl.getSolutionHead(), 
                     directory_ =results_dir,  allranks = True)
    write_file_array('errorChange_dofCoord'+str(cyl.gId),cyl.getDofCoordinates(), 
                     directory_ =results_dir, fileType = '.csv', allranks = True)
    write_file_array('errorChange_points'+str(cyl.gId),cyl.getPoints(), 
                     directory_ =results_dir, fileType = '.csv', allranks = True)
    write_file_array('errorChange_oldC'+str(cyl.gId),nc_content, 
                     directory_ =results_dir, fileType = '.csv', allranks = True)
    write_file_array('errorChange_newC'+str(cyl.gId), nc_content_new,
                     directory_ =results_dir, fileType = '.csv', allranks = True)
    write_file_array('errorChange_vols'+str(cyl.gId), cyl_cell_volumes,
                     directory_ =results_dir, fileType = '.csv', allranks = True )
    write_file_array('errorChange_newWatMol'+str(cyl.gId), newWatMol,
                     directory_ =results_dir, fileType = '.csv', allranks = True)
    write_file_array('errorChange_nc_molFr'+str(cyl.gId), nc_molFr,
                     directory_ =results_dir, fileType = '.csv', allranks = True)
    
    
def printFPitData(perirhizalModel, s, plantModel, fpit_Helper, rs_age_i_dt):
    """
        TODO: fix name objects
    """
    results_dir = perirhizalModel.results_dir
    write_file_array("inner_error", perirhizalModel.errs, directory_ =results_dir, fileType = '.csv') 
    results_dir = results_dir+'fpit/'
    if not perirhizalModel.doMinimumPrint:
        for lId, cyl in enumerate(perirhizalModel.cyls):
            if not isinstance(cyl, AirSegment):
                gId = perirhizalModel.eidx[lId]

                pHead = np.array(cyl.getSolutionHead()).flatten()
                write_file_array("fpit_Cyl_watercontent_"+str(gId),cyl.getWaterContent(), 
                                 directory_ =results_dir+'cyl_val/', allranks = True)
                write_file_array("fpit_Cyl_pressureHead_"+str(gId),pHead, 
                                 directory_ =results_dir+'cyl_val/', allranks = True)
                if max(pHead) > 0:
                    print('issue phead',gId,rank, pHead)
                    raise Exception
                    

            
        if rank == 0.:
        
            write_file_array("fpit_kr_soil", fpit_Helper.soilK,
                         directory_ =results_dir, fileType = '.csv') 
                         
            krplant =np.array([ plantModel.kr_f(rs_age_i_dt - plantModel.rs.nodeCTs[int(plantModel.rs.segments[seg_ind].y)], 
                        int(plantModel.rs.subTypes[seg_ind]), 
                        2, 
                        seg_ind, 
                        # cells = False to add
                        ) for seg_ind in range(len(plantModel.rs.segments))])  # c++ conductivity call back functions
            write_file_array("fpit_kr_plant",krplant,
                         directory_ =results_dir, fileType = '.csv') 
                         
                         
            write_file_array("fpit_seg_fluxes_limited", perirhizalModel.seg_fluxes_limited,
                         directory_ =results_dir, fileType = '.csv') 
            write_file_array("fpit_seg_fluxes", np.array(fpit_Helper.seg_fluxes),
                         directory_ =results_dir, fileType = '.csv') 
            write_file_array("fpit_seg_fluxes_limited_Out",perirhizalModel.seg_fluxes_limited_Out,
                         directory_ =results_dir, fileType = '.csv')
            write_file_array("fpit_proposed_outer_fluxes", fpit_Helper.proposed_outer_fluxes,
                         directory_ =results_dir, fileType = '.csv')

            write_file_array("fpit_errorsEach1DSWLim",
                             perirhizalModel.errorsEachWLim,
                             directory_ =results_dir, fileType = '.csv')
            write_file_array("fpit_errorsEach1DSW",
                             perirhizalModel.errorsEachW,
                             directory_ =results_dir, fileType = '.csv')

            write_file_array("fpit_psi_sri_oldvsnew",
                             perirhizalModel.errWrsis,
                             directory_ =results_dir, fileType = '.csv')

            write_file_array("fpit_psi_sri_realvsinput",
                             perirhizalModel.errWrsiRealInputs,
                             directory_ =results_dir, fileType = '.csv')
            write_file_array("fpit_psi_sri_input", fpit_Helper.rsx_input, directory_ =results_dir, fileType = '.csv') 
            write_file_array("fpit_psi_sri_real", fpit_Helper.rsx_old, directory_ =results_dir, fileType = '.csv')   
        
            write_file_array("fpit_errbulkMass",
                         np.array([s.bulkMassCErrorPlant_abs,
                                   s.bulkMassCErrorPlant_rel, #not cumulative
                                   s.bulkMassCError1ds_abs,
                                   s.bulkMassCError1ds_rel, 
                                   s.bulkMassErrorWater_abs,
                                   s.bulkMassErrorWater_rel,
                                   s.bulkMassErrorWater_absLim,
                                   s.bulkMassErrorWater_relLim,
                                   s.errSoil_source_sol_abs, 
                                   s.errSoil_source_sol_rel]), 
                         directory_ =results_dir, fileType = '.csv')   
            write_file_array("fpit_errorMassRhizo", np.array([perirhizalModel.rhizoMassCError_abs, perirhizalModel.rhizoMassCError_rel,
                                                perirhizalModel.rhizoMassWError_abs, perirhizalModel.rhizoMassWError_rel]), 
                         directory_ =results_dir, fileType = '.csv')# not cumulativecumulative (?)  



            write_file_array("fpit_sol_content_diff1d3dabs"+str(0), perirhizalModel.allDiff1d3dCW_abs[0],
                         directory_ =results_dir, fileType = '.csv')
            write_file_array("fpit_sol_content_diff1d3drel"+str(0), perirhizalModel.allDiff1d3dCW_rel[0], 
                         directory_ =results_dir, fileType = '.csv')

            write_file_array("fpit_outer_R_bc_wat", fpit_Helper.outer_R_bc_wat, directory_ =results_dir, fileType = '.csv') 
            write_file_array("fpit_rx", fpit_Helper.rx_old, directory_ =results_dir, fileType = '.csv') 

                         
            write_file_array("thetaCyl_4splitSoilVals", fpit_Helper.thetaCyl_4splitSoilVals,
                             directory_ =perirhizalModel.results_dir, fileType = '.csv')
                             
                             
            write_file_array("fpit_diffBCS1dsFluxOut", fpit_Helper.diffBCS1dsFluxOut,
                             directory_ =perirhizalModel.results_dir, fileType = '.csv')
                             
            write_file_array("fpit_flow1d1dw", 
                             perirhizalModel.flow1d1d_w, 
                             directory_ =results_dir, fileType = '.csv') 
            write_file_array("fpit_flow1d1dsol", 
                             perirhizalModel.flow1d1d_sol, 
                             directory_ =results_dir, fileType = '.csv') 
            write_file_array("fpit_flow1d1dmucil", 
                             perirhizalModel.flow1d1d_mucil, 
                             directory_ =results_dir, fileType = '.csv') 
        if False:# (not doMinimumPrint): this part still needs to be updated
                
            write_file_array("fpit_soil_fluxes_limited", soil_fluxes_limited, 
                             directory_ =results_dir, fileType = '.csv') 
            write_file_array("fpit_soil_fluxes", soil_fluxes,
                             directory_ =results_dir, fileType = '.csv') 
                         


            if (plantType != "plant") :
                write_file_array('fpit_transRate',np.array([transpiration]), directory_ =results_dir,
                                 fileType = '.csv' )
                write_file_array("fpit_n_iter",np.array([ n_iter, perirhizalModel.solve_gave_up ]), 
                                 directory_ =results_dir, fileType = '.csv') 



            # solutes. only limited vs unlimited for the rhizo: not relevent for soil source
            # not sure i need to have it for the inner BC as plant suc flow outside of simloop and we should only have 
            # exud, never uptake.
            #
            write_file_array("fpit_seg_fluxes_limited_sol_Out", seg_fluxes_limited_sol_Out,
                             directory_ =results_dir, fileType = '.csv') 
            write_file_array("fpit_seg_fluxes_sol_Out", proposed_outer_sol_fluxes,
                             directory_ =results_dir, fileType = '.csv') 

            write_file_array("fpit_seg_fluxes_limited_mucil_Out",
                             seg_fluxes_limited_mucil_Out,
                             directory_ =results_dir, fileType = '.csv') 
            write_file_array("fpit_seg_fluxes_mucil_Out", proposed_outer_mucil_fluxes,
                             directory_ =results_dir, fileType = '.csv') 

            write_file_array("fpit_seg_fluxes_limited_sol_In", seg_fluxes_limited_sol_In, 
                             directory_ =results_dir, fileType = '.csv') 
            write_file_array("fpit_seg_fluxes_sol_In", seg_sol_fluxes,
                             directory_ =results_dir, fileType = '.csv') 

            write_file_array("fpit_seg_fluxes_limited_mucil_In", seg_fluxes_limited_mucil_In,
                             directory_ =results_dir, fileType = '.csv') 
            write_file_array("fpit_seg_fluxes_mucil_In", seg_mucil_fluxes,
                             directory_ =results_dir, fileType = '.csv') 



            

            for nc in range(perirhizalModel.numSoluteComp):

                write_file_array("fpit_sol_content_diff1d3dabs"+str(nc+1), perirhizalModel.allDiff1d3dCW_abs[nc+1],
                                 directory_ =results_dir, fileType = '.csv')
                write_file_array("fpit_sol_content_diff1d3drel"+str(nc+1), perirhizalModel.allDiff1d3dCW_rel[nc+1], 
                                 directory_ =results_dir, fileType = '.csv')
                write_file_array("fpit_sol_content_3dabs"+str(nc+1), perirhizalModel.contentIn3d[nc+1],
                                 directory_ =results_dir, fileType = '.csv')
                write_file_array("fpit_sol_content_1dabs"+str(nc+1), perirhizalModel.contentIn1d[nc+1],
                                 directory_ =results_dir, fileType = '.csv')

                write_file_array("fpit_sol_content3d"+str(nc+1), 
                                 s.getContent(nc+1),
                                 directory_ =results_dir, fileType = '.csv')  
                write_file_float("fpit_sol_content1d"+str(nc+1), 
                                 sum(perirhizalModel.getContentCyl(idComp = nc+1, doSum = False, reOrder = True)), 
                                 directory_ =results_dir)#, fileType = '.csv')  
                write_file_array("fpit_outer_R_bc_sol"+str(nc+1), outer_R_bc_sol[nc],
                                 directory_ =results_dir, fileType = '.csv')  
                write_file_array("fpit_diffouter_R_bc_sol"+str(nc+1), diffouter_R_bc_sol[nc],
                                 directory_ =results_dir, fileType = '.csv') 

                write_file_array("fpit_bulkMassCError1dsAll_real"+str(nc+1), s.bulkMassCError1dsAll_real[nc], directory_ =results_dir, fileType = '.csv') 
                write_file_array("fpit_sources_sol_real"+str(nc+1), sources_sol[nc], directory_ =results_dir, fileType = '.csv') 


     

            write_file_array("fpit_all1d3dDiff",perirhizalModel.all1d3dDiff, directory_ =results_dir, fileType = '.csv') 
            write_file_array("fpit_all1d3dDiffBis",perirhizalModel.all1d3dDiff[cellIds], directory_ =results_dir, fileType = '.csv') 

            write_file_array("fpit_new_soil_water", new_soil_water, directory_ =results_dir, fileType = '.csv') 
            write_file_array("fpit_rhizoWaterPerVoxel", rhizoWaterPerVoxel, directory_ =results_dir, fileType = '.csv') 
            
            write_file_array("fpit_sumErrors1ds3dsAbs", perirhizalModel.sumDiff1d3dCW_abs, directory_ =results_dir, fileType = '.csv')
            write_file_array("fpit_sumErrors1ds3dsRel", perirhizalModel.sumDiff1d3dCW_rel, directory_ =results_dir, fileType = '.csv')
            write_file_array("fpit_maxErrors1ds3dsAbs", perirhizalModel.maxDiff1d3dCW_abs, directory_ =results_dir, fileType = '.csv')
            write_file_array("fpit_maxErrors1ds3dsRel", perirhizalModel.maxDiff1d3dCW_rel, directory_ =results_dir, fileType = '.csv')

            write_file_array("fpit_rhizoWAfter", rhizoWAfter_[rhizoSegsId] , directory_ =results_dir, fileType = '.csv') 
            write_file_array("fpit_rhizoTotCAfter", rhizoTotCAfter_[rhizoSegsId] , directory_ =results_dir, fileType = '.csv')
            
        

def WarnErrorIncrease( IdComponentWithError, perirhizalModel):
    """
        print data for debug and throw error
        @param: index of the ocmponent(s) with a high relative and absolute error
        TODO: make better error message
    """
    results_dir = perirhizalModel.results_dir

    print('perirhizalModel.sumDiff1d3dCW_rel , perirhizalModel.sumDiff1d3dCW_relOld',
          perirhizalModel.sumDiff1d3dCW_rel , perirhizalModel.sumDiff1d3dCW_relOld, perirhizalModel.diff1d3dCurrant_rel,
         np.floor(max(perirhizalModel.sumDiff1d3dCW_rel - perirhizalModel.sumDiff1d3dCW_relOld)),
          max(perirhizalModel.sumDiff1d3dCW_rel - perirhizalModel.sumDiff1d3dCW_relOld),

          perirhizalModel.sumDiff1d3dCW_rel - perirhizalModel.sumDiff1d3dCW_relOld, 
          ' perirhizalModel.sumDiff1d3dCW_abs', perirhizalModel.sumDiff1d3dCW_abs, 
          IdComponentWithError,  perirhizalModel.sumDiff1d3dCW_abs[IdComponentWithError] )
    write_file_array("error_sumDiff1d3dCW_relOld",
                     perirhizalModel.sumDiff1d3dCW_relOld, 
                     directory_ =results_dir, 
                     fileType = '.csv') 
    write_file_array("error_sumDiff1d3dCW_absOld",
                     perirhizalModel.sumDiff1d3dCW_absOld, 
                     directory_ =results_dir, 
                     fileType = '.csv') 
    write_file_array("error_sumDiff1d3dCW_rel", perirhizalModel.sumDiff1d3dCW_rel, 
                     directory_ =results_dir, fileType = '.csv') 
    write_file_array("error_sumDiff1d3dCW_abs", perirhizalModel.sumDiff1d3dCW_abs, 
                     directory_ =results_dir, fileType = '.csv') 
    write_file_float("error_diff1d3dCurrant_rel",
                     perirhizalModel.diff1d3dCurrant_rel, 
                     directory_ =results_dir) 

    write_file_array("error_maxDiff1d3dCW_relOld",
                     perirhizalModel.maxDiff1d3dCW_relOld, 
                     directory_ =results_dir, 
                     fileType = '.csv') 
    write_file_array("error_maxDiff1d3dCW_absOld",
                     perirhizalModel.maxDiff1d3dCW_absOld, 
                     directory_ =results_dir, 
                     fileType = '.csv') 
    write_file_array("error_maxDiff1d3dCW_rel", perirhizalModel.maxDiff1d3dCW_rel, 
                     directory_ =results_dir, fileType = '.csv') 
    write_file_array("error_maxDiff1d3dCW_abs", perirhizalModel.maxDiff1d3dCW_abs, 
                     directory_ =results_dir, fileType = '.csv') 
    write_file_float("error_maxdiff1d3dCurrant_rel",
                     perirhizalModel.maxdiff1d3dCurrant_rel, 
                     directory_ =results_dir) 

    write_file_array("fpit_sol_content_diff1d3dabs"+str(0), perirhizalModel.allDiff1d3dCW_abs[0], 
                     directory_ =results_dir, fileType = '.csv')
    write_file_array("fpit_sol_content_diff1d3drel"+str(0), perirhizalModel.allDiff1d3dCW_rel[0], 
                     directory_ =results_dir, fileType = '.csv')

    for nc in range(perirhizalModel.numSoluteComp):# normally all 0 for nc >= numDissolvedSoluteComp
            write_file_array("error_sol_content_diff1d3dabs"+str(nc+1), 
                             perirhizalModel.allDiff1d3dCW_abs[nc+1], directory_ =results_dir, fileType = '.csv')
            write_file_array("error_sol_content_diff1d3drel"+str(nc+1), 
                             perirhizalModel.allDiff1d3dCW_rel[nc+1], directory_ =results_dir, fileType = '.csv')
    raise Exception