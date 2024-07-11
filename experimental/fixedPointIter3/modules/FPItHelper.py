
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
import PhloemPhotosynthesis
import helpfull

class fixedPointIterationHelper():
    """
        little class to store data and do checks for the 
        inner iteration loop
    """
    def __init__(self, s, perirhizalModel, plantModel, 
                seg_fluxes, 
                 outer_R_bc_wat, outer_R_bc_sol, 
                 cylVol, 
                 Q_Exud_i, Q_mucil_i, dt,sim_time, emptyCells):
        self.s = s
        self.perirhizalModel = perirhizalModel
        self.plantModel = plantModel
        self.numSegs = len(np.asarray(plantModel.rs.organTypes, int))
        
        self.airSegsId = self.perirhizalModel.airSegs
        self.cylVol = cylVol
        self.emptyCells = emptyCells
        self.cell_volumes = comm.bcast(s.getCellVolumes() , root = 0)
        self.cellIds = perirhizalModel.cellWithRoots # only id of cells with roots
        
        
        # plant-soil solute flow, defined outside of iteration loop
        self.Q_Exud_i = Q_Exud_i
        self.Q_mucil_i = Q_mucil_i
        self.seg_sol_fluxes = Q_Exud_i /sim_time# mol/day for segments
        self.seg_mucil_fluxes = Q_mucil_i/sim_time
        self.sim_time = sim_time
        self.initializeOldData(seg_fluxes, outer_R_bc_wat, outer_R_bc_sol)
        
    def initializeOldData(self, seg_fluxes, outer_R_bc_wat, outer_R_bc_sol):
        """ create old data object to evaluate convergence """
        # data
        self.rx_old = 0 # plant water potential [cm]
        self.soil_water3dAfter_old = 0 # 
        self.totC3dAfter_eachVoxeleachComp_old = 0 # colute content [mol]
        self.rhizoWAfter_eachCyl_old = np.full(self.numSegs,0.) # water content in 1d models at the end of the time step
        self.rhizoWAfter_eachCyl = self.perirhizalModel.getWaterVolumesCyl(doSum = False, reOrder = True) #cm3 water per 1d model
        self.rhizoTotCAfter_eachCyl_old = np.full(self.numSegs,0.) # solute content in 1d models at the end of the time step
        self.seg_fluxes_old = 0 # plant-soil water exchanges
        self.proposed_outer_fluxes_old = 0 # 3d to 1d water flux
        self.proposed_outer_sol_fluxes_old = 0 # 3d to 1d small solute flux (1st dissolved component)
        self.proposed_outer_mucil_fluxes_old = 0 # 3d to 1d large solute (==mucilage) flux (2nd dissolved component)
        self.seg_fluxes = seg_fluxes.copy()
        self.outer_R_bc_wat_old =  0 # inter-cell water flux in 3d soil model
        self.outer_R_bc_sol_old = 0 # inter-cell solute flux in 3d soil model
        self.outer_R_bc_wat =  outer_R_bc_wat.copy() # inter-cell water flux in 3d soil model
        self.outer_R_bc_sol = outer_R_bc_sol.copy() # inter-cell solute flux in 3d soil model
        
        self.thetaCylOld = self.perirhizalModel.getWaterVolumesCyl(doSum = False, reOrder = True)/self.cylVol # cm3/cm3
        
        # matric potential at the segment-exterior interface, i.e. inner values of the (air or soil) cylindric models 
        self.rsx_init  = comm.bcast(self.perirhizalModel.get_inner_heads(weather=self.perirhizalModel.weatherX), root=0) # store value at beginning time step
        self.rsx_input = self.rsx_init # rsx used as input
        self.rsx_old = self.rsx_input.copy() # rsx used to compute convergence rate
        self.rsx_olds = []#[self.rsx_input.copy()] # rsx used to compute convergence rate
        
        # exchanges between cylinders in same 3d soil voxel
        self.perirhizalModel.flow1d1d_w = np.zeros(self.numSegs)
        self.perirhizalModel.flow1d1d_sol = np.zeros(self.numSegs)
        self.perirhizalModel.flow1d1d_mucil = np.zeros(self.numSegs)
        
        
        gotError, IdComponentWithError =  self.bigErrorIncrease()
        if gotError:
            printData.WarnErrorIncrease(IdComponentWithError, self.perirhizalModel)
            
        ### loop and error values
        self.n_iter = 0 # number of iteration
        self.err = 1.e6 
        self.max_err = self.perirhizalModel.max_err 
        self.perirhizalModel.rhizoMassWError_abs =1.# 
        self.perirhizalModel.rhizoMassCError_abs =1.# 
        self.perirhizalModel.errDiffBCs = np.array([1.])
        self.perirhizalModel.solve_gave_up = False # did we get at least 1 non-convergence error in dumux?
        self.perirhizalModel.diff1d3dCurrant_rel =1e6
        self.perirhizalModel.maxdiff1d3dCurrant_rel =1e6
        
        # store old cumulative error, to get instantenuous (!= cumulative) 1d3d error
        self.perirhizalModel.sumDiff1d3dCW_absOld = self.perirhizalModel.sumDiff1d3dCW_abs 
        self.perirhizalModel.sumDiff1d3dCW_relOld = self.perirhizalModel.sumDiff1d3dCW_rel  
        self.perirhizalModel.maxDiff1d3dCW_absOld = self.perirhizalModel.maxDiff1d3dCW_abs  
        self.perirhizalModel.maxDiff1d3dCW_relOld = self.perirhizalModel.maxDiff1d3dCW_rel  

        
    
    def bigErrorIncrease(self):
        """ check that newRelativeError - oldRelativeError <= 1%
            look at the difference between 1d and 3d model output
            sum of the difference in all the voxels + maximum difference
            no errors: return False, 0
            issue: return True and ids of the components with higher error
        
        """
        perirhizalModel = self.perirhizalModel
        if ( ((np.floor(max(perirhizalModel.sumDiff1d3dCW_rel - perirhizalModel.sumDiff1d3dCW_relOld)) > 1.) \
                or (np.floor(max(perirhizalModel.maxDiff1d3dCW_rel - perirhizalModel.maxDiff1d3dCW_relOld)) > 1.))):
            issueComp = np.where(np.floor((perirhizalModel.sumDiff1d3dCW_rel - perirhizalModel.sumDiff1d3dCW_relOld)) > 1.)
            # it s ok to have higher relative error if the absolute error is very small
            if (perirhizalModel.sumDiff1d3dCW_abs[issueComp] > 1e-13).any():
                return True, issueComp
        return False, 0
    
          
    def distribute3dto1dFlows(self, rs_age_i_dt, dt):
        """
            the inter-cell flow in the 3d soil is distributed between the 1d model contained by the 
            cells according to their water or solute content
        """
        self.getCyldatafor3d1dFlow(rs_age_i_dt, dt)
        outer_R_bc_wat = self.outer_R_bc_wat
        outer_R_bc_sol = self.outer_R_bc_sol
        s = self.s
        perirhizalModel = self.perirhizalModel
        
        # TODO: needed? 
        # ATT: do not move this to the end of the iteration loop
        # if new soil cells get roots after the iteration loop, we need the flux data
        if len(self.emptyCells) > 0:
            outer_R_bc_wat[self.emptyCells] = 0.
            for nc in range(perirhizalModel.numDissolvedSoluteComp):
                outer_R_bc_sol[nc][self.emptyCells] = 0.# only use BC for cells with at least one root
        
        
        if rank == 0:
            if max(abs(outer_R_bc_wat )) > 0:
                assert outer_R_bc_wat.shape == ( s.numberOfCellsTot, )
                try:
                    proposed_outer_fluxes = perirhizalModel.splitSoilVals(soilVals=outer_R_bc_wat / dt,
                                                        seg_values=self.thetaCyl_4splitSoilVals, 
                                                       seg_volume= self.cylVol,dt = dt,
                                                       isWater = True, 
                                                                      verbose = False) #cm3/day
                except:
                    print('thetaCyl_4splitSoilVals',min(self.thetaCyl_4splitSoilVals),max(self.thetaCyl_4splitSoilVals),
                          'self.thetaCylOld',min(self.thetaCylOld),max(self.thetaCylOld))
                    raise Exception                   
                                                 
            else:
                proposed_outer_fluxes = np.full(self.numSegs, 0.)   
                
            rhizoWAfter_eachCyl4splitVals = self.rhizoWAfter_eachCyl.copy()
            rhizoWAfter_eachCyl4splitVals[np.array(perirhizalModel.organTypes) !=2]=1. # to avoind division by 0.

            if max(abs(outer_R_bc_sol[0] )) > 0:
                proposed_outer_sol_fluxes = perirhizalModel.splitSoilVals(soilVals=outer_R_bc_sol[0] / dt, 
                                                seg_values=self.comp1content/rhizoWAfter_eachCyl4splitVals, dt = dt,
                                                seg_volume= rhizoWAfter_eachCyl4splitVals.copy(), isWater = False)#mol/day
            else:
                proposed_outer_sol_fluxes = np.full(self.numSegs, 0.)
                
                
            if max(abs(outer_R_bc_sol[1] )) > 0:
                proposed_outer_mucil_fluxes = perirhizalModel.splitSoilVals(soilVals=outer_R_bc_sol[1] / dt, 
                                              seg_values=self.comp2content/ rhizoWAfter_eachCyl4splitVals, 
                                                                            dt = dt,
                                               seg_volume=rhizoWAfter_eachCyl4splitVals.copy(),
                                                                            isWater = False)
            else:
                proposed_outer_mucil_fluxes = np.full(self.numSegs, 0.)
                
        else:
            proposed_outer_fluxes = None
            proposed_outer_sol_fluxes = None
            proposed_outer_mucil_fluxes = None


        self.proposed_outer_fluxes = comm.bcast(proposed_outer_fluxes, root = 0)
        self.proposed_outer_sol_fluxes = comm.bcast(proposed_outer_sol_fluxes, root = 0)
        self.proposed_outer_mucil_fluxes = comm.bcast(proposed_outer_mucil_fluxes, root = 0)
        
        
    def getCyldatafor3d1dFlow(self, rs_age_i_dt, dt):
        perirhizalModel = self.perirhizalModel
        
        # when there is not transpiration, use data at the beginning of the time step for water flow
        if (perirhizalModel.beforeAtNight and PhloemPhotosynthesis.noTranspiration(perirhizalModel, rs_age_i_dt, dt)):
            self.thetaCyl_4splitSoilVals = self.thetaCylOld
        else:
            # get data before doing the 'reset' => values at the end of the time step
            # self.thetaCyl = perirhizalModel.getWaterVolumesCyl(doSum = False, reOrder = True)/self.cylVol # cm3 water
            # get all potentially available water == water after reset + water taken up by plant
            seg_fluxes_root  = self.seg_fluxes.copy()
            seg_fluxes_root[np.where(np.array(perirhizalModel.organTypes)!=2)]  = 0.
            self.thetaCyl_4splitSoilVals = (self.thetaCylOld * self.cylVol + (seg_fluxes_root + perirhizalModel.flow1d1d_w)* dt )/self.cylVol
            #because of unlimited flow (seg_fluxes_root + perirhizalModel.flow1d1d_w), might get values out of the [theta_r, theta_s] bounds
            self.thetaCyl_4splitSoilVals = np.maximum(np.minimum(
                                                    self.thetaCyl_4splitSoilVals, 
                                                      perirhizalModel.vg_soil.theta_S),
                                                     perirhizalModel.theta_wilting_point)
            self.thetaCyl_4splitSoilVals[np.where(np.array(perirhizalModel.organTypes)!=2)]  = 0.
            
       
        # get data before doing the 'reset' => values at the end of the time step
        self.comp1content = perirhizalModel.getContentCyl(idComp=1, doSum = False, reOrder = True) # [mol] small rhizodeposits
        self.comp2content = perirhizalModel.getContentCyl(idComp=2, doSum = False, reOrder = True) # [mol] mucilage

        # assert (self.thetaCyl_4splitSoilVals >= 0.).all() # because of seg_fluxes, could have thetaCyl_4splitSoilVals < 0
        assert (self.comp1content >= 0.).all()
        assert (self.comp2content >= 0.).all()
        
        
    def do1d1dFlow(self):
        """
            wrapper function for computation of inter-cylinder flow
            use data after doing the reset, otherwise get diverging results
        """
        perirhizalModel = self.perirhizalModel
        soilVals_ = comm.bcast( np.array(self.s.getSolutionHead()),root= 0) #<= data beginning time step
        rhizoWAfter_eachCyl_divide = self.rhizoWAfter_eachCyl.copy()
        rhizoWAfter_eachCyl_divide[np.where(rhizoWAfter_eachCyl_divide==0.)] = 1. # air segments
        if perirhizalModel.do1d1dFlow:
            perirhizalModel.do1d1dFlow_(soilVals = soilVals_,#wat. pot. (cm)
                    seg_valuesW =  self.rhizoWBefore_eachCyl/self.cylVol,#self.rhizoWBefore_eachCyl/self.cylVol,#wat. pot. (cm) <= data beginning time step
                    seg_valuesSol = self.comp1content/rhizoWAfter_eachCyl_divide,# mol/cm3 wat <= data end time step
                    seg_valuesmucil = self.comp2content/rhizoWAfter_eachCyl_divide,# mol/cm3 wat <= data end time step
                    verbose=False)
            
        # when we don t have a dynamic soil for water, set the inter-cylinder flow to 0
        gotDynamicSoil = ((perirhizalModel.spellData['scenario'] == 'none') or 
                          ((perirhizalModel.spellData['scenario'] != 'baseline') and 
                           (perirhizalModel.enteredSpell) and 
                           (not perirhizalModel.leftSpell)))
        if not gotDynamicSoil: # only keep solute flow
            perirhizalModel.flow1d1d_w *= 0.
        
    def storeLimitedFluxes(self):
        """
            store actual BC of 1d models (to compare with prescribed BC)
            also compute actual - limited flux error
        """
        perirhizalModel = self.perirhizalModel
        # inner water flux
        self.seg_fluxes_limited = perirhizalModel.getXcyl(data2share=perirhizalModel.seg_fluxes_limited,idCyll_=None, doSum = False, reOrder = True) 
        # outer water flux 
        self.seg_fluxes_limited_Out = perirhizalModel.getXcyl(data2share=perirhizalModel.seg_fluxes_limited_Out,
                                                              idCyll_=None, doSum = False, reOrder = True) 
        self.seg_fluxes_limited_sol_Out = perirhizalModel.getXcyl(data2share=perirhizalModel.seg_fluxes_limited_sol_Out,
                                                                  idCyll_=None, doSum = False, reOrder = True) 
        self.seg_fluxes_limited_mucil_Out = perirhizalModel.getXcyl(data2share=perirhizalModel.seg_fluxes_limited_mucil_Out,
                                                                    idCyll_=None, doSum = False, reOrder = True) 
        self.seg_fluxes_limited_sol_In = perirhizalModel.getXcyl(data2share=perirhizalModel.seg_fluxes_limited_sol_In,
                                                                 idCyll_=None, doSum = False, reOrder = True)
        self.seg_fluxes_limited_mucil_In = perirhizalModel.getXcyl(data2share=perirhizalModel.seg_fluxes_limited_mucil_In,
                                                                   idCyll_=None, doSum = False, reOrder = True) 

                
        if len(self.airSegsId)>0:                
            try:
                assert (self.seg_fluxes_limited[self.airSegsId] == self.seg_fluxes[self.airSegsId]).all()
            except:
                print('seg_fluxes_limited vs seg_flux', self.seg_fluxes_limited[self.airSegsId] - self.seg_fluxes[self.airSegsId])
                raise Exception

        # get limited vs real boundary fluxes. if we always have a difference, the loop failes        
        perirhizalModel.SinkLim1DS =max( abs((self.seg_fluxes_limited - self.seg_fluxes)/ 
                                            np.where(self.seg_fluxes,self.seg_fluxes,1.))*100.) # at the end of the fixed point iteration, should be ~ 0 
        perirhizalModel.OutLim1DS =max( abs((self.seg_fluxes_limited_Out - self.proposed_outer_fluxes)/ 
                                            np.where(self.proposed_outer_fluxes,
                                            self.proposed_outer_fluxes,1.))*100.) # at the end of the fixed point iteration, should be ~ 0 
        perirhizalModel.InOutBC_Cdiff = []# at the end of the fixed point iteration, should be ~ 0 
        fluxes_limitedAll = [self.seg_fluxes_limited_sol_Out, self.seg_fluxes_limited_mucil_Out,
                                self.seg_fluxes_limited_sol_In, self.seg_fluxes_limited_mucil_In]
        fluxesAll = [self.proposed_outer_sol_fluxes, self.proposed_outer_mucil_fluxes,
                                self.seg_sol_fluxes, self.seg_mucil_fluxes]
        for soluteflux in range(4):
            fluxes_limited = fluxes_limitedAll[soluteflux]
            fluxes = fluxesAll[soluteflux]
            relError = max( abs((fluxes_limited - fluxes)/ 
                                            np.where(fluxes,fluxes,1.))*100.)
            perirhizalModel.InOutBC_Cdiff.append(relError) 
        perirhizalModel.InOutBC_Cdiff = np.array(perirhizalModel.InOutBC_Cdiff)
    
    def compute1dChangesSolute(self, dt):
        """
            compute the net source of solutes per voxel of the 3d soil  
            NOT cuased by inter-cell flow in bulk soil 
            from the content variation simulated by the 1d soil models
            (= plant-soil exchange + biochemical reactions)
            TODO: get sources_sol_from1d directly from dumux instead of
            re-computing it afterwards
        """
        s = self.s
        cellIds = self.cellIds
        # get changes in colute content not caused by flow [mol/day]
        # todo: get it directly from dumux
        sources_sol_from1d = np.full(self.soil_solute1d_perVoxelAfter.shape,0. )
        for nc in range(self.perirhizalModel.numSoluteComp):
            sources_sol_from1d[nc][cellIds] = np.array(
                self.soil_solute1d_perVoxelAfter[nc][cellIds] - self.soil_solute1d_perVoxelBefore[nc][cellIds] - self.outer_R_bc_sol[nc][cellIds]
            )/dt
        sources_sol_from1d = comm.bcast(sources_sol_from1d, root = 0)    
        
        assert sources_sol_from1d.shape == (self.perirhizalModel.numSoluteComp, self.s.numberOfCellsTot)
        
        # store error
        self.s.errSoil_source_sol_abs = sum(sources_sol_from1d.flatten()) - (sum(self.Q_Exud_i) + sum(self.Q_mucil_i))/self.sim_time
        
        if (sum(self.Q_Exud_i) + sum(self.Q_mucil_i))/dt != 0.:
            s.errSoil_source_sol_rel = abs(s.errSoil_source_sol_abs/((sum(self.Q_Exud_i) + sum(self.Q_mucil_i))/self.sim_time)*100)
        else:
            s.errSoil_source_sol_rel = np.nan
        self.sources_sol_from1d = sources_sol_from1d
            
    def compute1dChangesWater(self, dt):
        plantModel = self.plantModel
        perirhizalModel = self.perirhizalModel
        s = self.s
        if (rank == 0):
            net_PWU_ = plantModel.sumSegFluxes(self.seg_fluxes)  #  plant water uptake [cm3/day]  per soil cell
            net_PWU = np.zeros(s.numberOfCellsTot)
            #easier to handle array than dict. maybe change sumSegFluxes() function to choose type of output
            net_PWU[np.array(list(net_PWU_.keys()))] = np.array(list(net_PWU_.values())) 
        else:
            net_PWU = None

        self.net_PWU = comm.bcast(net_PWU, root=0)

        if (rank == 0):
            net_PWU_limited_ = plantModel.sumSegFluxes(self.seg_fluxes_limited)  # [cm3/day]  per soil cell
            net_PWU_limited = np.zeros(s.numberOfCellsTot)
            #easier to handle array than dict. maybe change sumSegFluxes() function to choose type of output
            net_PWU_limited[np.array(list(net_PWU_limited_.keys()))] = np.array(list(net_PWU_limited_.values()))
        else:
            net_PWU_limited = None

        self.net_PWU_limited = comm.bcast(net_PWU_limited, root=0)
        #print('compute1dChangesWater', 'seg_', sum(self.seg_fluxes), sum(net_PWU),sum(self.seg_fluxes)- sum(net_PWU),
        #                    'seg_fluxes_limited', sum(self.seg_fluxes_limited), sum(net_PWU_limited) ,
        #                    sum(self.seg_fluxes_limited)- sum(net_PWU_limited) ,
        #                  'diff',  sum(net_PWU_limited)- sum(net_PWU),sum(self.seg_fluxes)- sum(self.seg_fluxes_limited))
        perirhizalModel.SinkLim3DS =max( abs((self.net_PWU_limited - self.net_PWU)/ 
                                            np.where(self.net_PWU,self.net_PWU,1.))*100.) # at the end of the fixed point iteration, should be ~ 0 
                                            
        
    def compute3DSource(self, dt):
        """
            compute source/sink of water and solutes in the 3d soil according to the changes in the 1d models
        """
        perirhizalModel = self.perirhizalModel
        results_dir = perirhizalModel.results_dir
        s = self.s
        
        self.compute1dChangesSolute(dt)
        self.compute1dChangesWater(dt)
        
        # source/sink in the 3d soil according to changes in the 1d models [(cm3 water or mol C)/day]
        soil_sources_limited = np.concatenate((np.array([self.net_PWU_limited]),self.sources_sol_from1d ))
        

        # when setting source: we limit according to the content at the beginning 
        # of the time step + max(net incoming flow, 0.) 
        # to get the maximum potential value of content that can be
        # in the voxel if needed
        
        for idComp in range(s.numComp):#cm3/day, mol/day

            SSL = soil_sources_limited[idComp].copy()
            
            # compute maximum poentially available content for sink (cm3 water or mol)
            if idComp == 0:
                maxPotentialAvailable = (np.array(self.soil_water3dBefore) +
                                     np.maximum(self.outer_R_bc_wat, self.outer_R_bc_wat*0 ))
            else:
                maxPotentialAvailable = (self.totC3dBefore_eachVoxeleachComp[idComp-1] +
                                         np.maximum(self.outer_R_bc_sol[idComp-1],
                                                     self.outer_R_bc_sol[idComp-1]*0 ) )

            if (max(abs(SSL)) != 0.):
                SSL = np.maximum(SSL, -maxPotentialAvailable/dt)
                
                # how much of the uptake do we need to redistribute?
                toAdd= np.maximum(0., -(maxPotentialAvailable/dt + SSL))
                # redistribute what is missing to 1d segments who have maxPotentialAvailable remaining
                SSL[np.where(toAdd>0.)] += toAdd[np.where(toAdd>0.)] 
                
                # do loop until all the sink has been redistributed
                k_limit_source3d = 0
                epsilon_source3d = 1e-25
                while (not (SSL*dt >= -maxPotentialAvailable).all()) and (k_limit_source3d <= 10):
                    SSL[np.where(SSL*dt < -maxPotentialAvailable)] += epsilon_source3d
                    epsilon_source3d *= 10
                    k_limit_source3d += 1

                try:
                    assert min(maxPotentialAvailable + SSL*dt) >=0.
                except:
                    print(soil_sources_limited[idComp], SSL,dt,  min(maxPotentialAvailable + SSL*dt) )
                    write_file_array("setsourceLim2_"+str(idComp), SSL, directory_ =results_dir, fileType =".csv") 
                    raise Exception
                    
            self.sendSource2dumux(SSL, idComp)
            
           
    def sendSource2dumux(self, SSL, idComp):
        s = self.s
        perirhizalModel=self.perirhizalModel
        results_dir = perirhizalModel.results_dir
        # convert array to dictionnary
        test_values = list(SSL)
        test_keys = np.array([i for i in range(len(test_values))])
        res = {}
        for key in test_keys:
            for value in test_values:
                res[key] = value
                test_values.remove(value)
                break                        

        if not self.perirhizalModel.doMinimumPrint: 
            write_file_array("setsourceLim2_"+str(idComp), SSL, directory_ =results_dir, fileType =".csv") 
        # send to dumux
        s.setSource(res.copy(), eq_idx = idComp)  # [mol/day], in modules/richards.py
        
    def solve3DS(self, dt):
        """
            wrapper around the solving of the 3d soil flow
            TODO: reset maxDt after creating the solving object.
        """
        perirhizalModel=self.perirhizalModel
        s = self.s
        
        k_soil_solve = 0
        redoSolve = True
        maxRelShift = s.MaxRelativeShift
        while redoSolve:
            #s.ddt =min( 1.e-5,s.ddt)#or just reset to 1e-5?
            #s.setMaxTimeStepSize(s.maxDt) # reset each time
            try:
                decreaseMaxRelShift = False
                if rank==0:
                    print("solve 3d soil")
                #s.solve(dt)  # in modules/solverbase.py
                helpfull.run_with_timeout(60.*15.,s.solve,dt) # time out error after Xmn
                if  rank==0:
                    print("solve 3d soil finished")
                
                # if we got solute content < 0, throw exception
                solComp = [s.getSolution(ncom+1) for ncom in range(s.numSoluteComp)]
                whereError = None
                if rank == 0:
                    whereError = [np.where(SC <0.) for SC in solComp]
                    solComp = [min(SC) for SC in solComp]
                solComp = comm.bcast(solComp, root = 0)
                whereError = comm.bcast(whereError, root = 0)
                if min(solComp) <0.:
                    print("min(solComp) <0.", rank, solComp, whereError)
                    decreaseMaxRelShift = True
                    raise Exception
                    
                redoSolve = False
                ## solving succeded, reset solving parameters (in case  they were changed)
                # newton parameters are re-read at each 'solve()' calls
                s.setParameter("Newton.MaxRelativeShift", str(s.MaxRelativeShift))# reset value
                s.setParameter("Newton.EnableResidualCriterion", "false") # sometimes helps, sometimes makes things worse
                s.setParameter("Newton.EnableAbsoluteResidualCriterion", "false")
                s.setParameter("Newton.SatisfyResidualAndShiftCriterion", "false")

                s.setParameter("Newton.MaxSteps", "18")
                s.setParameter("Newton.MaxTimeStepDivisions", "10")
                s.createNewtonSolver() # re-create Newton solver to implement the new newton parameters
                
            except Exception as err:
                s.setParameter("Newton.EnableResidualCriterion", "false") # sometimes helps, sometimes makes things worse
                s.setParameter("Newton.EnableAbsoluteResidualCriterion", "false")
                s.setParameter("Newton.SatisfyResidualAndShiftCriterion", "false")
                
                print(rank, f"Unexpected {err=}, {type(err)=}", 'k_soil_solve',k_soil_solve)
                
                if k_soil_solve > 6: # to avoid going to infinit loop
                    raise Exception
                if k_soil_solve == 0: 
                    # for the 1st fail, simply increase number of steps allowed
                    s.setParameter("Newton.MaxSteps", "200")
                    s.setParameter("Newton.MaxTimeStepDivisions", "100")
                elif k_soil_solve == 1: 
                    # 2nd fail: try making the computation more precise
                    print(rank,
                          'soil.solve() failed. making the computation more precise')
                    s.setParameter("Newton.EnableResidualCriterion", "true") # sometimes helps, sometimes makes things worse
                    s.setParameter("Newton.EnableAbsoluteResidualCriterion", "true")
                    s.setParameter("Newton.SatisfyResidualAndShiftCriterion", "true")
                    s.setParameter("Newton.MaxRelativeShift", str(s.MaxRelativeShift/10.))# reset value
                else:
                    # try decreasing MaxRelShift (if got solute content < 0)
                    # try increasing maxrelshift (if got a non-convergence error)
                    if decreaseMaxRelShift:
                        change = 0.1
                    else:
                        change = 10
                    print(rank,
                          'soil.solve() failed. NewtonMaxRelativeShift from',
                          maxRelShift,'to',maxRelShift*change)
                    maxRelShift *= change
                    # newton parameters are re-read at each 'solve()' calls
                    s.setParameter("Newton.MaxRelativeShift", str(maxRelShift))
                s.reset() # reset solution vector
                s.createNewtonSolver() # re-create Newton solver to implement the new newton parameters
                k_soil_solve += 1
                
        
    def storeOldMassData1d(self):
        """
            store mass data of the 1d models beginning time step
            to evaluate later mole balance error
            @see  massBalanceError1d
        """
        perirhizalModel = self.perirhizalModel
        s = self.s
        self.rhizoWBefore_eachCyl = perirhizalModel.getWaterVolumesCyl(doSum = False, reOrder = True) #cm3 water per 1d model
        self.rhizoWBefore = sum(self.rhizoWBefore_eachCyl) 
        self.rhizoTotCBefore_eachCyl = perirhizalModel.getTotCContentAll(doSum = False, reOrder = True) # total carbon content per 1d model
        self.rhizoTotCBefore = sum(self.rhizoTotCBefore_eachCyl) 

        # return component content per voxel according to 1d model data 
        # used to compute 3d sink
        self.soil_solute1d_perVoxelBefore = np.array( [np.array(perirhizalModel.getC_rhizo(
                                                                        idComp = idc + 1, konz = False
                                                                    )) for idc in range(perirhizalModel.numSoluteComp)])
        
    def storeNewMassData1d(self):
        """
            store mass data of the 1d models end time step
            to evaluate later mole balance error
            @see  massBalanceError1d
        """ 
        perirhizalModel = self.perirhizalModel
        s = self.s
        self.rhizoWAfter_eachCyl = perirhizalModel.getWaterVolumesCyl(doSum = False, reOrder = True) #cm3 water per 1d model
        self.rhizoWAfter = sum(self.rhizoWBefore_eachCyl) 
        self.rhizoTotCAfter_eachCyl = perirhizalModel.getTotCContentAll(doSum = False, reOrder = True) # total carbon content per 1d model
        self.rhizoTotCAfter = sum(self.rhizoTotCAfter_eachCyl) 
        # return component content per voxel according to 1d model data 
        # used to compute 3d sink
        self.soil_solute1d_perVoxelAfter = np.array( [np.array(perirhizalModel.getC_rhizo(
                                                                              idComp = idc + 1, konz = False
                                                                             )) for idc in range(perirhizalModel.numSoluteComp)])
        
        
        try: # todo: take out?
            assert min(self.soil_solute1d_perVoxelAfter.flatten()) >=0
        except:
            print("min(self.soil_solute1d_perVoxelAfter)",
                  min(self.soil_solute1d_perVoxelAfter.flatten()),
                  [min(nss) for nss in self.soil_solute1d_perVoxelAfter])
            raise Exception
            
        # evaluate convergence for water and solute content 1d models
        self.errW1ds = np.linalg.norm(self.rhizoWAfter_eachCyl - self.rhizoWAfter_eachCyl_old)
        self.errC1ds = np.linalg.norm(self.rhizoTotCAfter_eachCyl - self.rhizoTotCAfter_eachCyl_old)
        self.rhizoWAfter_eachCyl_old = self.rhizoWAfter_eachCyl.copy()
        self.rhizoTotCAfter_eachCyl_old = self.rhizoTotCAfter_eachCyl.copy()

            
    def massBalanceError1d(self, dt):
        """
            get mass balance error of 1d models
        """
        perirhizalModel = self.perirhizalModel
        rhizoSegsId = perirhizalModel.rhizoSegsId # plant segments with a rhizosphere model
        airSegsId = self.airSegsId
        
        ############ solutes
        # get error according to the 'limited' (== realised ) boundary fluxes
        # should be always ~ 0 
        errorsEachC = self.rhizoTotCAfter_eachCyl - ( 
            self.rhizoTotCBefore_eachCyl + (self.seg_fluxes_limited_sol_In+ self.seg_fluxes_limited_mucil_In+\
                                                              self.seg_fluxes_limited_sol_Out+ self.seg_fluxes_limited_mucil_Out)*dt)
        # store absolute total error for limited flow
        perirhizalModel.rhizoMassCError_absLim = sum(abs(errorsEachC[rhizoSegsId]))
        
        
        # get error according to the proposed (==prescribed) flux
        # need to be ~ 0 when leaving fixed point iteration 
        errorsEachC = self.rhizoTotCAfter_eachCyl - ( self.rhizoTotCBefore_eachCyl + 
                                                     (self.seg_sol_fluxes+ 
                                                      self.proposed_outer_sol_fluxes+ 
                                                      perirhizalModel.flow1d1d_sol +
                                                      self.seg_mucil_fluxes+ 
                                                      self.proposed_outer_mucil_fluxes+ 
                                                      perirhizalModel.flow1d1d_mucil )*dt)
        
        perirhizalModel.rhizoMassCError_abs  = sum(abs(errorsEachC[rhizoSegsId]))
        
        # store relative total error 
        if sum(self.rhizoTotCAfter_eachCyl) != 0:
            perirhizalModel.rhizoMassCError_relLim = abs(perirhizalModel.rhizoMassCError_absLim/sum(self.rhizoTotCAfter_eachCyl)*100)
            perirhizalModel.rhizoMassCError_rel = abs(perirhizalModel.rhizoMassCError_abs/sum(self.rhizoTotCAfter_eachCyl)*100)
        else:
            perirhizalModel.rhizoMassCError_relLim = np.nan
            perirhizalModel.rhizoMassCError_rel = np.nan
            

        
        ############ water
        # get error according to the 'limited' (== realised ) boundary fluxes
        # should be always ~ 0 
        errorsEachW = self.rhizoWAfter_eachCyl - ( 
            self.rhizoWBefore_eachCyl + (self.seg_fluxes_limited + self.seg_fluxes_limited_Out)*dt)
        
        
        # store absolute total error for limited flow
        perirhizalModel.rhizoMassWError_absLim = sum(abs(errorsEachW[rhizoSegsId]))
        perirhizalModel.errorsEachWLim = errorsEachW

        # get error according to the proposed (==prescribed) flux
        # need to be ~ 0 when leaving fixed point iteration 
        perirhizalModel.errorsEachW = self.rhizoWAfter_eachCyl - ( self.rhizoWBefore_eachCyl + (self.seg_fluxes+ self.proposed_outer_fluxes+ perirhizalModel.flow1d1d_w)*dt)
        perirhizalModel.rhizoMassWError_abs  = sum(abs(perirhizalModel.errorsEachW[rhizoSegsId]))
        
        # store relative total error 
        perirhizalModel.rhizoMassWError_relLim = abs(perirhizalModel.rhizoMassWError_absLim/sum(self.rhizoWAfter_eachCyl)*100)
        perirhizalModel.rhizoMassWError_rel = abs(perirhizalModel.rhizoMassWError_abs/sum(self.rhizoWAfter_eachCyl)*100)
                
        # to be sure that the value is the same in all the threads 
        # (I'm not sure when/if thread 0 is the only one with the important infor)
        # perirhizalModel.rhizoMassWError_abs = comm.bcast(perirhizalModel.rhizoMassWError_abs,root= 0)
        # perirhizalModel.rhizoMassCError_abs = comm.bcast(perirhizalModel.rhizoMassCError_abs,root= 0)
        
        if rank == 0:
            print(f'relative error balance soil 1d (%)?\n\t\tfrom PWU: {perirhizalModel.rhizoMassWError_rel:.2e}, from PWU-limited: {perirhizalModel.rhizoMassWError_relLim:.2e}, from PCU: {perirhizalModel.rhizoMassCError_rel:.2e}, from PCU-limited: {perirhizalModel.rhizoMassCError_relLim:.2e}'
                )

        
    def storeOldMassData3d(self):
        """
            store mass data of the 3d model at beginning time step
            to evaluate later mole balance error
            @see  massBalanceError3d
        """
        s = self.s
        water_content = comm.bcast( np.array(s.getWaterContent()),root= 0)  # theta per cell [1]
        self.soil_water3dBefore = np.multiply(water_content, self.cell_volumes)  # water per cell [cm3]
        
        # 3d soil solute content per solute type and voxel
        self.totC3dBefore_eachVoxeleachComp = comm.bcast(s.getTotCContent_each(), root = 0) 
        self.totC3dBefore = np.sum(self.totC3dBefore_eachVoxeleachComp)
        assert isinstance(self.totC3dBefore ,float) # just test once: we did the complete sum?
        
    def storeNewMassData3d(self):
        """
            store mass data of the 3d model at end time step
            to evaluate later mole balance error
            @see  massBalanceError3d
        """
        s = self.s
        water_content = comm.bcast( np.array(s.getWaterContent()),root= 0)  # theta per cell [1]
        self.soil_water3dAfter = np.multiply(water_content, self.cell_volumes)  # water per cell [cm3]
        
        # 3d soil solute content per solute type and voxel
        self.totC3dAfter_eachVoxeleachComp = comm.bcast(s.getTotCContent_each(), root = 0) 
        self.totC3dAfter = np.sum(self.totC3dAfter_eachVoxeleachComp)
        #assert isinstance(self.totC3dAfter ,float) # just test once: we did the complete sum?
        
        
        # convergence water and solute content 3d models
        self.errW3ds = np.linalg.norm(self.soil_water3dAfter - self.soil_water3dAfter_old)
        self.errC3ds = np.linalg.norm(self.totC3dAfter_eachVoxeleachComp - self.totC3dAfter_eachVoxeleachComp_old)
        self.soil_water3dAfter_old = self.soil_water3dAfter.copy()
        self.totC3dAfter_eachVoxeleachComp_old = self.totC3dAfter_eachVoxeleachComp.copy()

            
    def massBalanceError3d(self, dt):
        s = self.s
        # according to plant data # I need to add the water flux to be sure
        s.bulkMassErrorWater_absLim = sum(abs(self.soil_water3dAfter - self.soil_water3dBefore  -self.sources_wat_from3d-self.outer_R_bc_wat))
        
        s.bulkMassErrorWater_abs = sum(abs(self.soil_water3dAfter -  self.soil_water3dBefore-self.net_PWU*dt-self.outer_R_bc_wat))
        s.bulkMassErrorWater_relLim = abs(s.bulkMassErrorWater_absLim /sum(self.soil_water3dAfter) )*100
        s.bulkMassErrorWater_rel = abs(s.bulkMassErrorWater_abs /sum(self.soil_water3dAfter) )*100
        s.bulkMassErrorWaterSink_abs = self.net_PWU_limited*dt- self.sources_wat_from3d # sink sent to dumux == sink implemented?


        s.bulkMassCErrorPlant_abs = abs(self.totC3dAfter - ( self.totC3dBefore + sum(self.Q_Exud_i) + sum(self.Q_mucil_i)))
        if self.totC3dAfter > 0:
            s.bulkMassCErrorPlant_rel = abs(s.bulkMassCErrorPlant_abs/self.totC3dAfter*100)
        else:
            s.bulkMassCErrorPlant_rel = np.nan

        # according to 1d data
        s.bulkMassCError1ds_abs = abs(self.totC3dAfter - ( self.totC3dBefore + sum(self.sources_sol_from1d.flatten())*dt))
        if self.totC3dAfter > 0:
            s.bulkMassCError1ds_rel = abs(s.bulkMassCError1ds_abs/self.totC3dAfter*100)
        else:
            s.bulkMassCError1ds_rel = np.nan
            


        if rank == 0:
            print(f"relative error balance soil 3d (%)?\n\t\tfrom PWU {s.bulkMassErrorWater_rel:.2e}, from PWU-limited {s.bulkMassErrorWater_relLim:.2e},from PCU {s.bulkMassCErrorPlant_rel:.2e}, from 1ds C-change {s.bulkMassCError1ds_rel:.2e}"
                )
                
        ### for each voxel
        s.bulkMassErrorWaterAll_abs = abs(self.soil_water3dAfter - (self.soil_water3dBefore  + self.sources_wat_from3d + self.outer_R_bc_wat))
        s.bulkMassErrorWaterAll_rel = abs(s.bulkMassErrorWaterAll_abs /self.soil_water3dAfter )*100







        ## inner mass balance error 

        s.bulkMassCError1dsAll_abs = abs(self.totC3dAfter_eachVoxeleachComp - (self.totC3dBefore_eachVoxeleachComp + 
                                                    self.sources_sol_from3d + self.outer_R_bc_sol))


        totC3dAfter_eachVoxeleachCompTemp = np.where(self.totC3dAfter_eachVoxeleachComp != 0,
                                                  self.totC3dAfter_eachVoxeleachComp,
                                                  1)
        s.bulkMassCError1dsAll_rel = abs(s.bulkMassCError1dsAll_abs*100/totC3dAfter_eachVoxeleachCompTemp)

        try:
            assert s.bulkMassCError1dsAll_rel.shape == (s.numSoluteComp, s.numberOfCellsTot )
        except:
            print(s.bulkMassCError1dsAll_rel.shape , 
                  (s.numSoluteComp, s.numberOfCellsTot ))
            raise Exception




    def computeConvergence(self):
        """
            long funciton, compute different types of errors:
            non-convergence, mass balance errors, prescribed vs realised BCs
        """
        perirhizalModel = self.perirhizalModel
        s = self.s
        plantModel = self.plantModel
        
        
        # convergence plant wat. pot
        rx = comm.bcast(np.array(plantModel.psiXyl), root = 0) 
        rx_divide = np.where(rx !=0, rx, 1.)
        errRxPlant = max(abs((rx - self.rx_old)/rx_divide)*100.)
        self.rx_old = rx.copy()
        
        # convergence wat. pot. at root-soil interface            
        rsx =comm.bcast( perirhizalModel.get_inner_heads(weather=perirhizalModel.weatherX), root = 0)  # works for all threads
        rsx_divide = np.where(abs(rsx) < abs(self.rsx_input),
                                  rsx,self.rsx_input)
        perirhizalModel.errWrsiRealInputs =  abs((rsx - self.rsx_input)/rsx_divide)*100.
        perirhizalModel.errWrsiRealInputs[abs(rsx - self.rsx_input)<= 5]= 0.
        errWrsiRealInputf = max(perirhizalModel.errWrsiRealInputs)
        if rank == 0:
            print("errWrsiRealInput:",errWrsiRealInputf)
        perirhizalModel.errWrsiRealInput = errWrsiRealInputf
        
        rsx_divide = np.where(rsx!=0.,rsx,1.)
        perirhizalModel.errWrsis =  abs((rsx - self.rsx_old)/rsx_divide)*100.#max()
        perirhizalModel.errWrsis[abs(rsx - self.rsx_old)<= 5]= 0.
        errWrsi = max(perirhizalModel.errWrsis)
        self.rsx_old = rsx;
        self.rsx_olds.append(self.rsx_old)
        perirhizalModel.errWrsi = errWrsi

        # convergence BC flow 1d models
        diffBCS1dsFluxIn =   np.array(self.seg_fluxes)  - self.seg_fluxes_old   #only for water as plant exud is outside of loop
        self.diffBCS1dsFluxOut =   self.proposed_outer_fluxes  - self.proposed_outer_fluxes_old 
        
                         
        diffBCS1dsFluxOut_sol =   self.proposed_outer_sol_fluxes  - self.proposed_outer_sol_fluxes_old 
        diffBCS1dsFluxOut_mucil =   self.proposed_outer_mucil_fluxes  - self.proposed_outer_mucil_fluxes_old

        errBCS1dsFluxIn =max(abs((
            diffBCS1dsFluxIn/ np.where(np.array(self.seg_fluxes),
                                       np.array(self.seg_fluxes),1.))*100))
        errBCS1dsFluxOut = max(abs((
            self.diffBCS1dsFluxOut/ np.where(self.proposed_outer_fluxes,
                                        self.proposed_outer_fluxes,1.))*100))
        
        errBCS1dsFluxOut_sol = max(abs((
            diffBCS1dsFluxOut_sol/ np.where(self.proposed_outer_sol_fluxes,
                                            self.proposed_outer_sol_fluxes,1.))*100))
        errBCS1dsFluxOut_mucil = max(abs((
            diffBCS1dsFluxOut_mucil/ np.where(self.proposed_outer_mucil_fluxes,
                                              self.proposed_outer_mucil_fluxes,1.))*100))
                                              
        self.seg_fluxes_old = np.array(self.seg_fluxes).copy()        
        self.proposed_outer_fluxes_old = self.proposed_outer_fluxes.copy()
        self.proposed_outer_sol_fluxes_old =self.proposed_outer_sol_fluxes.copy()
        self.proposed_outer_mucil_fluxes_old =self.proposed_outer_mucil_fluxes.copy()
        
        # convergence inter-cell flow for 3d soil
        diffouter_R_bc_wat =   self.outer_R_bc_wat  - self.outer_R_bc_wat_old 
        diffouter_R_bc_sol =   self.outer_R_bc_sol  - self.outer_R_bc_sol_old 

        errOuter_R_bc_wat = max(abs((diffouter_R_bc_wat/ np.where(self.outer_R_bc_wat,
                                                                  self.outer_R_bc_wat,1.))*100))
        errOuter_R_bc_sol = max(abs((diffouter_R_bc_sol[:1].reshape(-1)/ np.where(self.outer_R_bc_sol[:1].reshape(-1),
                                                                                               self.outer_R_bc_sol[:1].reshape(-1),1.))*100))
        
        self.outer_R_bc_wat_old = self.outer_R_bc_wat.copy()
        self.outer_R_bc_sol_old = self.outer_R_bc_sol.copy()
        
        
        #### store error values

        perirhizalModel.maxDiff1d3dCW_abs =np.array( comm.bcast(perirhizalModel.maxDiff1d3dCW_abs,root= 0))

        # only look at the relative error if the absolute error is high enough
        compErrorAboveLim = np.where(perirhizalModel.sumDiff1d3dCW_abs > 1e-13) 

        # to not depend on cumulative error
        perirhizalModel.diff1d3dCurrant = max(np.append((perirhizalModel.sumDiff1d3dCW_abs - perirhizalModel.sumDiff1d3dCW_absOld)[compErrorAboveLim],0.)) 
        perirhizalModel.diff1d3dCurrant_rel =max(np.append((perirhizalModel.sumDiff1d3dCW_rel - perirhizalModel.sumDiff1d3dCW_relOld)[compErrorAboveLim],0.))

        perirhizalModel.maxdiff1d3dCurrant = max(np.append((perirhizalModel.maxDiff1d3dCW_abs - perirhizalModel.maxDiff1d3dCW_absOld)[compErrorAboveLim],0.)) 
        perirhizalModel.maxdiff1d3dCurrant_rel =max(np.append((perirhizalModel.maxDiff1d3dCW_rel - perirhizalModel.maxDiff1d3dCW_relOld)[compErrorAboveLim],0.))

        perirhizalModel.solve_gave_up = (np.array(comm.bcast(comm.gather(perirhizalModel.solve_gave_up ,root = 0),root = 0))).any()
        
        # one metric to decide if we stay in the iteration loop or not
        
        perirhizalModel.err = comm.bcast(max(errRxPlant,self.errW1ds, self.errW3ds,self.errC1ds, self.errC3ds,
                                             errWrsi,errWrsiRealInputf,
                                s.bulkMassErrorWater_rel, 
                               s.bulkMassCErrorPlant_rel, s.bulkMassCError1ds_rel,
                               perirhizalModel.rhizoMassWError_rel,
                               perirhizalModel.rhizoMassCError_rel,
                               #perirhizalModel.SinkLim3DS,
                        #perirhizalModel.SinkLim1DS
                        ),root= 0)
        # one array to do printing
        perirhizalModel.errs =np.array([
                        # non-convergence metrics
                        errRxPlant, self.errW1ds, self.errW3ds,self.errC1ds, self.errC3ds, errWrsi,errWrsiRealInputf,
                        errBCS1dsFluxIn, errBCS1dsFluxOut,errBCS1dsFluxOut_sol, errBCS1dsFluxOut_mucil,
                        errOuter_R_bc_wat, errOuter_R_bc_sol,
                        # realised vs prescribed fluxes and sinks
                        perirhizalModel.SinkLim3DS,
                        perirhizalModel.SinkLim1DS,perirhizalModel.OutLim1DS, 
                        perirhizalModel.InOutBC_Cdiff[0],perirhizalModel.InOutBC_Cdiff[1],perirhizalModel.InOutBC_Cdiff[2],perirhizalModel.InOutBC_Cdiff[3],
                       # 1d-3d differences/errors
                        max(perirhizalModel.sumDiff1d3dCW_abs),max(perirhizalModel.sumDiff1d3dCW_rel), perirhizalModel.diff1d3dCurrant,perirhizalModel.diff1d3dCurrant_rel, 
                        max(perirhizalModel.maxDiff1d3dCW_abs), max(perirhizalModel.maxDiff1d3dCW_rel), perirhizalModel.maxdiff1d3dCurrant,perirhizalModel.maxdiff1d3dCurrant_rel, 
                        # mass balance error 3d model
                        s.bulkMassErrorWater_abs,s.bulkMassErrorWater_rel,
                        s.bulkMassErrorWater_absLim,s.bulkMassErrorWater_relLim,
                        s.bulkMassCErrorPlant_abs, s.bulkMassCError1ds_abs,
                        s.bulkMassCErrorPlant_rel, s.bulkMassCError1ds_rel,  
                       # mass balance error 1d models
                       perirhizalModel.rhizoMassWError_relLim, perirhizalModel.rhizoMassCError_relLim,
                       perirhizalModel.rhizoMassWError_rel, perirhizalModel.rhizoMassCError_rel,
                       # summary metric
                       perirhizalModel.err ])

        