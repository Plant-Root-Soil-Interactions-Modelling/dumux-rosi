
import numpy as np
import pandas as pd
import timeit
import pandas as pd
import matplotlib; matplotlib.use('agg')
import sys;
import os
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank(); max_rank = comm.Get_size()
import timeit

from air_modelsPlant import AirSegment
from helpfull import write_file_array, write_file_float, div0, div0f

def initialPrint(plant):
    """ put headings on files which will be filled later """
    errs = np.array(["errRxPlant", "errW1ds", "errW3ds","errC1ds", "errC3ds",
                     "max(r.SinkLim3DS)","max(r.SinkLim1DS)","max(abs(r.OutLim1DS))",
                     "max(abs(r.InOutBC_Cdiff))",
                     "max(r.maxDiff1d3dCW_abs)",
                     "errWrsi",
                     "bulkMassErrorWater_abs","bulkMassErrorWater_absLim",
                     "rhizoMassWError_absLim","rhizoMassWError_abs",
                     "bulkMassErrorC_abs","bulkMassCErrorPlant",
                     "rhizoMassCError_absLim","rhizoMassCError_abs",
                     "sum(abs(diffBCS1dsFluxIn))", "sum(abs(diffBCS1dsFluxOut))",
                     "sum(abs(diffouter_R_bc_wat))",
                     "sum(abs(diffBCS1dsFluxOut_sol))",
                     "sum(abs(diffBCS1dsFluxOut_mucil))","sum(abs(diffouter_R_bc_sol))",
                     "diff1d3dCurrant","diff1d3dCurrant_rel","rhizoMassWError_rel",'err'])
    write_file_array("OuterSuccess_error", errs, directory_ =plant.results_dir, fileType = '.csv')
    
    if not plant.doMinimumPrint:# print data during fixed-point iteration?
        write_file_array("N_error", errs, directory_ =plant.results_dir, fileType = '.csv')
        write_file_array("fpit_error", errs, directory_ =plant.results_dir, fileType = '.csv') 
        
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
                             s.getContent(nc+1, nc < s.numDissolvedSoluteComp), 
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
    results_dir = rs.results_dir
    write_file_array("time", np.array([rs_age,0.]), directory_ =results_dir)
    write_file_array("sumErrors1ds3ds", np.concatenate((rs.sumDiff1d3dCW_abs, rs.sumDiff1d3dCW_rel)), directory_ =results_dir, fileType = '.csv')
    write_file_array("maxErrors1ds3ds", np.concatenate((rs.maxDiff1d3dCW_abs, rs.maxDiff1d3dCW_rel)), directory_ =results_dir, fileType = '.csv') 
    
    
def printPhloemFlowOutput(rs_age, perirhizalModel, phloemDataStorage, plantModel):
    """ print and save outputs of phloem flow and photosynthesis"""
    results_dir = perirhizalModel.results_dir
    
    print("\n\n\n\t\tat ", int(rs_age//1),"d", int(rs_age%1*24),"h", int(rs_age%1*24%1*60),"mn",  
                  round(plantModel.Qlight *1e6),"mumol m-2 s-1")
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
    write_file_float("Q_Exud", phloemDataStorage.Q_Exud, directory_ =results_dir)
    write_file_float("Q_Mucil", phloemDataStorage.Q_Mucil, directory_ =results_dir)
    write_file_float("Q_Exud_tot", phloemDataStorage.Q_Exud_inflate, directory_ =results_dir)
    write_file_float("Q_Mucil_tot", phloemDataStorage.Q_Mucil_inflate, directory_ =results_dir)
    
    write_file_float("Q_Ag", phloemDataStorage.Q_in, directory_ =results_dir)
    write_file_array("psiXyl", plantModel.psiXyl, directory_ =results_dir)
    write_file_array("Fpsi", plantModel.Fpsi, directory_ =results_dir)
    write_file_array("fw", plantModel.fw, directory_ =results_dir)
    write_file_array("gco2", plantModel.gco2, directory_ =results_dir)
    write_file_array("Q_Ag_dot", plantModel.AgPhl, directory_ =results_dir)
    write_file_array("C_rsi", np.array(plantModel.Csoil_seg ), 
                     directory_ =results_dir)#mmol/cm3
    write_file_float("trans", plantModel.TranspirationCumul, directory_ =results_dir)
    write_file_array("transrate",plantModel.Jw, directory_ =results_dir, fileType = '.csv')
    
    write_file_array("errorsPlant", np.array([phloemDataStorage.error_st_abs,
                                              phloemDataStorage.error_st_rel,#cumulative
                                             abs(sum(plantModel.outputFlux))]), 
                     directory_ =results_dir, fileType = '.csv') #not cumulative

def printCylData(perirhizalModel):
    """ print data of each perirhizal model"""
    results_dir = perirhizalModel.results_dir
    
    for lId, cyl in enumerate(perirhizalModel.cyls):
        if not isinstance(cyl, AirSegment):
            gId = perirhizalModel.eidx[lId]

            pHead = np.array(cyl.getSolutionHead()).flatten()
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
            for ccc in range(perirhizalModel.numSoluteComp):
                sol0 = np.array(cyl.getContent(ccc +1, ccc < perirhizalModel.numDissolvedSoluteComp)).flatten()
                write_file_array("Cyl_content"+str(ccc+1)+"_"+str(gId)+"", 
                             sol0, 
                             directory_ =results_dir+'cyl_val/', allranks = True)
            if max(pHead) > 0:
                print('issue phead',gId,rank, pHead, sol0 )
                raise Exception
                
def getAndPrintErrorRates(perirhizalModel, plantModel, s, phloemData):
    """
        define and save other error rates
    """
    results_dir = perirhizalModel.results_dir
    write_file_array("OuterSuccess_error", perirhizalModel.errs, 
                     directory_ =results_dir, fileType = '.csv') 
    write_file_array("OuterSuccess_data", 
                     np.array([perirhizalModel.n_iter, perirhizalModel.err, 
                               perirhizalModel.rhizoMassWError_abs,
                               perirhizalModel.dt_inner]), 
                     directory_ =results_dir, fileType = '.csv')

    buTotCAfter = sum(s.getTotCContent())   #0 get stuck here

    buWAfter = sum(np.multiply(np.array(s.getWaterContent()), s.getCellVolumes()))    

    if rank == 0:
        if perirhizalModel.doSoluteFlow:
            s.bulkMassErrorCumul_abs = abs((buTotCAfter - ( s.buTotCSoilInit + 
                                        phloemData.Q_Exud_inflate + 
                                        phloemData.Q_Mucil_inflate)))
            if buTotCAfter != 0:
                s.bulkMassErrorCumul_rel = abs(s.bulkMassErrorCumul_abs/buTotCAfter*100)
            else:
                s.bulkMassErrorCumul_rel =np.nan
        else:
            s.bulkMassErrorCumul_abs = np.nan
            s.bulkMassErrorCumul_rel =np.nan
        # if we have a free flow BC at the bottom, that could increase the error
        # ideally, I should add the flow at the bellow BC here
        s.bulkMassErrorWaterCumul_abs = abs(buWAfter - ( s.buWSoilInit - plantModel.TranspirationCumul_eval))
        s.bulkMassErrorWaterCumul_rel = abs(s.bulkMassErrorWaterCumul_abs/buWAfter*100)
    else:
        s.bulkMassErrorCumul_abs = None
        s.bulkMassErrorCumul_rel = None
        s.bulkMassErrorWaterCumul_abs = None
        s.bulkMassErrorWaterCumul_rel = None

    if perirhizalModel.doSoluteFlow:
        #absolute and relative (%) error
        write_file_array("errorsBulkSoil", np.array([s.bulkMassCErrorPlant_abs, s.bulkMassCErrorPlant_rel, #not cumulative 
                                            s.bulkMassCError1ds_abs, s.bulkMassCError1ds_rel, 
                                            s.bulkMassErrorCumul_abs,s.bulkMassErrorCumul_rel,#cumulative
                                            s.bulkMassErrorWater_abs,s.bulkMassErrorWater_rel, #not cumulative
                                            s.bulkMassErrorWaterCumul_abs,s.bulkMassErrorWaterCumul_rel]),
                         directory_ =results_dir, fileType = '.csv')#cumulative
        write_file_array("errorMassRhizo", np.array([perirhizalModel.rhizoMassCError_abs, perirhizalModel.rhizoMassCError_rel,
                                                    perirhizalModel.rhizoMassWError_abs, perirhizalModel.rhizoMassWError_rel]), directory_ =results_dir)# not cumulativecumulative (?)


def doVTPplots(perirhizalModel, s, soilTextureAndShape):
    """ vtp plots for water and C components in 3d soil and at root-soil
        interface
    """
    results_dir = perirhizalModel.results_dir
    # otherwise the shoot looks weird
    periodic = False
    
    
    extraArray_ = perirhizalModel.soilModel.getWaterContent()
    extraArrayName_ = "theta (cm3/cm3)"
    
    rp = perirhizalModel.get_concentration(0, konz)
    
    # TODO: adapt to have plot_plants_and soil
    vp.plot_roots_and_soil(perirhizalModel.mappedSegments(),extraArrayName_,rp, s, periodic, 
                           soilTextureAndShape['min_b'],
                           soilTextureAndShape['max_b'],
                           soilTextureAndShape['cell_number'], 
            filename="soil_rx",sol_ind =-1,
                           extraArray = extraArray_, extraArrayName = extraArrayName_,
            interactiveImage=False)  # VTK vizualisation
    print("idcomp_done", 0, rank)
    for i in range(1, perirhizalModel.numComp):
        print("idcomp", i, rank)
        extraArray_ = perirhizalModel.soilModel.getSolution(i) * perirhizalModel.phaseDensity(i)/1e6

        if not konz:
            extraArrayName_ = "C"+str(i)+" mol"
        else:
            extraArrayName_ = "[C"+str(i)+"] (mol/cm3)"
            if i <= perirhizalModel.numDissolvedSoluteComp:
                extraArray_ /= (perirhizalModel.soilModel.getWaterContent() *cell_volumes) #cm3 water
            else: 
                extraArray_ /= cell_volumes #cm3


        vp.plot_roots_and_soil(perirhizalModel.mappedSegments(),
                               extraArrayName_ ,
                               perirhizalModel.get_concentration(i , konz), s, 
                               periodic,
                           soilTextureAndShape['min_b'],
                           soilTextureAndShape['max_b'],
                               soilTextureAndShape['cell_number'], 
                filename="C"+str(i), 
                               sol_ind =-1,
                               extraArray = extraArray_, 
                extraArrayName = extraArrayName_,
            interactiveImage=False)  # VTK vizualisation
        print("idcomp_done", i, rank)
        
        
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
    
    
def printFPitData(perirhizalModel, s, plantModel, fpit_Helper):
    write_file_array("inner_error", perirhizalModel.errs, directory_ =results_dir, fileType = '.csv') 
    if  (not doMinimumPrint):
        write_file_array("fpit_errbulkMass",
                         np.array([s.bulkMassCErrorPlant_abs,
                                   s.bulkMassCErrorPlant_rel, #not cumulative
                                   s.bulkMassCError1ds_abs,
                                   s.bulkMassCError1ds_rel, 
                                   s.bulkMassErrorWater_abs,
                                   s.bulkMassErrorWater_rel,
                                   s.bulkMassCErrorPlant_absReal,
                                   s.errSoil_source_sol_abs, 
                                   s.errSoil_source_sol_rel]), 
                         directory_ =results_dir, fileType = '.csv') 
        write_file_array("fpit_errorMassRhizo", np.array([perirhizalModel.rhizoMassCError_abs, perirhizalModel.rhizoMassCError_rel,
                                                perirhizalModel.rhizoMassWError_abs, perirhizalModel.rhizoMassWError_rel]), 
                         directory_ =results_dir, fileType = '.csv')# not cumulativecumulative (?)
        
        write_file_array("fpit_flow1d1dw", 
                         perirhizalModel.flow1d1d_w, 
                         directory_ =results_dir, fileType = '.csv') 
        write_file_array("fpit_flow1d1dsol", 
                         perirhizalModel.flow1d1d_sol, 
                         directory_ =results_dir, fileType = '.csv') 
        write_file_array("fpit_flow1d1dmucil", 
                         perirhizalModel.flow1d1d_mucil, 
                         directory_ =results_dir, fileType = '.csv') 


        if (plantType != "plant") :
            write_file_array('fpit_transRate',np.array([transpiration]), directory_ =results_dir,
                             fileType = '.csv' )
            write_file_array("fpit_n_iter",np.array([ n_iter, perirhizalModel.solve_gave_up ]), 
                             directory_ =results_dir, fileType = '.csv') 


        write_file_array("fpit_seg_fluxes_limited", seg_fluxes_limited,
                         directory_ =results_dir, fileType = '.csv') 
        write_file_array("fpit_seg_fluxes", np.array(seg_fluxes),
                         directory_ =results_dir, fileType = '.csv') 
        write_file_array("fpit_seg_fluxes_limited_Out",seg_fluxes_limited_Out,
                         directory_ =results_dir, fileType = '.csv')
        write_file_array("fpit_proposed_outer_fluxes", proposed_outer_fluxes,
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


        write_file_array("fpit_soil_fluxes_limited", soil_fluxes_limited, 
                         directory_ =results_dir, fileType = '.csv') 
        write_file_array("fpit_soil_fluxes", soil_fluxes,
                         directory_ =results_dir, fileType = '.csv') 

        write_file_array("fpit_sol_content_diff1d3dabs"+str(0), perirhizalModel.allDiff1d3dCW_abs[0],
                         directory_ =results_dir, fileType = '.csv')
        write_file_array("fpit_sol_content_diff1d3drel"+str(0), perirhizalModel.allDiff1d3dCW_rel[0], 
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
                             s.getContent(nc+1, nc < s.numDissolvedSoluteComp),
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




        write_file_array("fpit_outer_R_bc_wat", outer_R_bc_wat, directory_ =results_dir, fileType = '.csv')  

        write_file_array("fpit_all1d3dDiff",perirhizalModel.all1d3dDiff, directory_ =results_dir, fileType = '.csv') 
        write_file_array("fpit_all1d3dDiffBis",perirhizalModel.all1d3dDiff[cellIds], directory_ =results_dir, fileType = '.csv') 

        write_file_array("fpit_psi_sri_input", rsx_input, directory_ =results_dir, fileType = '.csv') 
        write_file_array("fpit_psi_sri_real", rsx, directory_ =results_dir, fileType = '.csv') 
        
        write_file_array("fpit_new_soil_water", new_soil_water, directory_ =results_dir, fileType = '.csv') 
        write_file_array("fpit_rhizoWaterPerVoxel", rhizoWaterPerVoxel, directory_ =results_dir, fileType = '.csv') 
        write_file_array("fpit_rx", rx, directory_ =results_dir, fileType = '.csv') 
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