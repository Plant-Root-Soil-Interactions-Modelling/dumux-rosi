
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
                                             abs(sum(plantModel.outputFlux))]), directory_ =results_dir, fileType = '.csv') #not cumulative

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