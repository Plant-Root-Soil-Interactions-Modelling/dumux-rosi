
    
import numpy as np
import pandas as pd
import timeit
import pandas as pd
import matplotlib; matplotlib.use('agg')
import sys;
import os
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank(); max_rank = comm.Get_size()
import timeit
    
from functional.xylem_flux import sinusoidal
from helpfull import sinusoidal3
from air_modelsPlant import AirSegment

    
def weather(simDuration, dt, spellData, hp:float=1):
        if simDuration == 0.:
            raise Exception
        Qnigh = 0; Qday = 960e-6 
        
        if  ((spellData['condition'] == "wet") or (simDuration <= spellData['spellStart']) or (simDuration > spellData['spellEnd'])):
            Tnigh = 15.8; Tday = 22
            RHday = 0.6; RHnigh = 0.88
            Pair = 1010.00 #hPa
            pmean = -100.
            cs = 350e-6
        elif spellData['condition'] == "dry":
            Tnigh = 20.7; Tday = 30.27
            RHday = 0.44; RHnigh = 0.78
            Pair = 1070.00 #hPa
            pmean = -450.
            cs = 350e-6
        else:
            print('spellData',spellData)
            raise Exception
            
        coefhours = sinusoidal3(simDuration + 0.5, dt)
        RH_ = RHnigh + (RHday - RHnigh) * coefhours
        TairC_ = Tnigh + (Tday - Tnigh) * coefhours
        Q_ = Qnigh + (Qday - Qnigh) * coefhours
        
        es =  6.112 * np.exp((17.67 * TairC_)/(TairC_ + 243.5))
        ea = es*RH_
        assert ea < es
        
        assert ((RH_ > 0) and(RH_ < 1))
        bl_thickness = 1/1000 #1mm * m_per_mm
        diffusivity= 2.5e-5#m2/sfor 25*C
        rbl =bl_thickness/diffusivity #s/m 13
        
        Kcanopymean = 1e-1 # m2/s
        meanCanopyL = (2/3) * hp /2
        rcanopy = meanCanopyL/Kcanopymean
        windSpeed = 2 #m/s
        zmzh = 2 #m
        karman = 0.41 #[-]
        
        rair = 1
        if hp > 0:
            rair = np.log((zmzh - (2/3)*hp)/(0.123*hp)) * np.log((zmzh - (2/3)*hp)/(0.1*hp)) / (karman*karman*windSpeed)

        weatherVar = {'TairC' : TairC_,'TairK' : TairC_ + 273.15,'Pair':Pair,"es":es,
                        'Qlight': Q_,'rbl':rbl,'rcanopy':rcanopy,'rair':rair,"ea":ea,
                        'cs':cs, 'RH':RH_, 'p_mean':pmean
                     }
        return weatherVar
    
def checkDoWeatherChange(perirhizalModel, rs_age_i_dt):
    """
        check if need to do the weather change
        @param: perirhizalModel object containing all the 1d models
    """
    # using scenario with a dry spell?
    if (perirhizalModel.spellData['scenario'] != 'none') and (perirhizalModel.spellData['scenario'] != 'baseline'):
        # did we just enter or leave the spell?            
        if ((rs_age_i_dt > perirhizalModel.spellData['spellEnd']) and (not perirhizalModel.leftSpell)):
            if rank == 0:
                print("leftSpell", rs_age_i_dt, perirhizalModel.spellData['spellStart'], perirhizalModel.spellData['spellEnd'])
            perirhizalModel.leftSpell = True
            return True
        if ((rs_age_i_dt > perirhizalModel.spellData['spellStart']) and (not perirhizalModel.enteredSpell)):
            if rank == 0:
                print("enteredSpell", rs_age_i_dt, perirhizalModel.spellData['spellStart'], perirhizalModel.spellData['spellEnd'])
            perirhizalModel.enteredSpell = True
            return True
    return False
    
    
def weatherChange(rs_age_i_dt, perirhizalModel, s):
    """
        function to change water content if we leave or enter spell period.
        deletethe water and solute gradient in the perirhizal zones
        @param: rs_age_i_dt. plant age [d]
        @param: perirhizalModel object containing all the 1d models
        @param: soil model
    """
    results_dir = perirhizalModel.results_dir
    if checkDoWeatherChange(perirhizalModel, rs_age_i_dt):
        cell_volumes = comm.bcast(s.getCellVolumes() , root = 0) #cm3    
        # get new wat. pot
        pheadinit_cm =  perirhizalModel.weatherX['p_mean'] # get new mean wat. pot [cm]
        cellsZ = comm.bcast(np.array( [ loc[2] for loc in s.getCellCenters()]))#cm
        meanZ = np.average(cellsZ)
        pheadinit_cm_all = pheadinit_cm - (cellsZ - meanZ) # get wat. pot. for each soil layer
        pheadinit_Pa = s.to_pa(pheadinit_cm_all)

        perirhizalModel.check1d3dDiff() # get error before changing data (for troubleshooting)
        if rank ==0:
            print('weather::weatherChange(): error before change',perirhizalModel.sumDiff1d3dCW_rel ,
                  perirhizalModel.sumDiff1d3dCW_abs,'pheadinit_cm',pheadinit_cm)
            
        #  save dissolved solute content (no need to change anything for solutes in soil phase)
        nc_content = np.array([comm.bcast(s.getContent(nc+1), root=0)  for nc in range(perirhizalModel.numDissolvedSoluteComp)])
        # send new wat. pot to dumux
        s.base.setSolution(pheadinit_Pa,0 )
        # new water content. [cm3 wat/cm3 scv] * [cm3 scv] * [m3/cm3] * [mol/m3 wat] = mol wat
        newWatMol = (comm.bcast(s.getWaterContent(),root = 0) * cell_volumes) * (1/1e6) * s.molarDensityWat_m3 
        # new solute mole fraction (mol C / mol water)
        nc_molFr =np.array( [nc_c/newWatMol for nc_c in nc_content])
        # set solution
        for nc in range(perirhizalModel.numDissolvedSoluteComp):
            s.base.setSolution(nc_molFr[nc],nc+1 )
        # get new content (check if same as old content)
        nc_content_new = np.array([comm.bcast(s.getContent(nc+1), root=0)  for nc in range(perirhizalModel.numDissolvedSoluteComp)])# mol
        
        if not perirhizalModel.doMinimumPrint:
            write_file_array('pheadinit_cm',pheadinit_cm_all, directory_ =results_dir, fileType = '.csv')
            write_file_array('newWatMol',newWatMol, directory_ =results_dir, fileType = '.csv')
            write_file_array('solContentBeforechange',nc_content, directory_ =results_dir, fileType = '.csv')
            write_file_array('solContentAfterchange',nc_content_new, directory_ =results_dir, fileType = '.csv')
            
        # same process for the cylinders
        for locIdCyl, cyl in enumerate(perirhizalModel.cyls):
            if not isinstance(cyl, AirSegment):
                pheadOld = cyl.getSolutionHead()
                cellId = perirhizalModel.seg2cell[cyl.gId]
                cyl_cell_volumes = cyl.getCellVolumesCyl()  #cm3 scv                            
                pheadinit_PaCyl = np.full(len(cyl_cell_volumes),pheadinit_Pa[cellId])
                
                nc_content = np.array([cyl.getContent(nc+1)  for nc in range(perirhizalModel.numDissolvedSoluteComp)])# mol
                cyl.base.setSolution(pheadinit_PaCyl,0 )

                newTheta = cyl.getWaterContent()
                newWatMol = (cyl.getWaterContent() * cyl_cell_volumes) * (1/1e6) * s.molarDensityWat_m3 # mol water
                nc_molFr =np.array( [nc_c/newWatMol for nc_c in nc_content]) # new solute mol fraction
                # set new mole fraction to keep same content
                for nc in range(perirhizalModel.numDissolvedSoluteComp):
                    cyl.base.setSolution(nc_molFr[nc],nc+1 )
                # new content (for check)
                nc_content_new = np.array([cyl.getContent(nc+1)  for nc in range(perirhizalModel.numDissolvedSoluteComp)])# mol

                try:
                    # allow for a small change
                    assert (abs(
                        nc_content.reshape(-1) - nc_content_new.reshape(-1)
                    ) <= np.minimum(abs(nc_content.reshape(-1)),
                                    abs( nc_content_new.reshape(-1)))*1e-6  
                           ).all()
                except:
                    printData.errorWeatherChange(results_dir, cyl, pheadOld,
                                                 nc_content, nc_content_new, nc_molFr)
                    raise Exception
                    
        s.save()# <= so that stay in the current weather state
        perirhizalModel.save() 
        s.buWSoilInit = sum(np.multiply(comm.bcast(np.array(s.getWaterContent()), root=0), cell_volumes)) # to evaluate cumulative water balance error
        # reset the inter-cell water fluxes to 0
        outer_R_bc_wat = np.full(cell_volumes.shape, 0.)
            
        # store old error rate and compare with the new one. 
        # allow for small increase
        beforeChange_sumDiff1d3dCW_abs = perirhizalModel.sumDiff1d3dCW_abs
        perirhizalModel.check1d3dDiff() 
        if rank ==0:
            print('weather::weatherChange(): error after change',
                  perirhizalModel.sumDiff1d3dCW_rel, perirhizalModel.sumDiff1d3dCW_abs )
        assert (perirhizalModel.sumDiff1d3dCW_abs <= beforeChange_sumDiff1d3dCW_abs*10).all()