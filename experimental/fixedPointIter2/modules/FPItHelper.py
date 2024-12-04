
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



def storeOldMassData1d(perirhizalModel):
    """
        store mass data of the 1d models beginning time step
        to evaluate later mole balance error
        @see  massBalanceError1d
    """
    perirhizalModel.rhizoWBefore_eachCyl = perirhizalModel.getWaterVolumesCyl(doSum = False, reOrder = True) #cm3 water per 1d model
    # each solute for each cyl.
    perirhizalModel.soil_solute1d_eachCylechSolBefore = perirhizalModel.getTotCContentAll(doSum = False, reOrder = True)
    
    if rank ==0:

        cellIds= perirhizalModel.getCellIds()
        # total carbon content per 1d model
        perirhizalModel.rhizoTotCBefore_eachCyl = perirhizalModel.soil_solute1d_eachCylechSolBefore.sum(axis=1) 
        perirhizalModel.rhizoWBefore = sum(perirhizalModel.rhizoWBefore_eachCyl) 
        
        perirhizalModel.rhizoTotCBefore = sum(perirhizalModel.rhizoTotCBefore_eachCyl) 

        # return component content per voxel according to 1d model data 
        # used to compute 3d sink
        perirhizalModel.soil_solute1d_perVoxelBefore = np.array([
            perirhizalModel.soil_solute1d_eachCylechSolBefore[perirhizalModel.getIdCyllMPI(cellId)[0]].sum(axis=0) for cellId in cellIds
        ]).T
        
        if perirhizalModel.debugMode:
            write_file_array("fpit_soil_solute1d_eachCylechSolBefore1", perirhizalModel.soil_solute1d_eachCylechSolBefore[1], 
                             directory_ =perirhizalModel.results_dir, fileType = '.txt') 

            write_file_array("fpit_rhizoWBefore_eachCyl", perirhizalModel.rhizoWBefore_eachCyl, 
                                 directory_ =perirhizalModel.results_dir, fileType = '.txt') 
    

def storeNewMassData1d(perirhizalModel):
    """
        store mass data of the 1d models end time step
        to evaluate later mole balance error
        @see  massBalanceError1d
    """ 
    perirhizalModel.rhizoWAfter_eachCyl = perirhizalModel.getWaterVolumesCyl(doSum = False, reOrder = True) #cm3 water per 1d model

    # each solute for each cyl.
    perirhizalModel.soil_solute1d_eachCylechSolAfter = perirhizalModel.getTotCContentAll(doSum = False, reOrder = True) # total carbon content per 1d model
    
    if rank ==0:
        cellIds= perirhizalModel.getCellIds()
       
        assert ((perirhizalModel.soil_solute1d_eachCylechSolAfter.shape == (len(perirhizalModel.eidx_all),perirhizalModel.numSoluteComp))|(perirhizalModel.soil_solute1d_eachCylechSolAfter.shape[0] == len(perirhizalModel.eidx_all) == 0)) 

        # total carbon content per 1d model
        perirhizalModel.rhizoTotCAfter_eachCyl = perirhizalModel.soil_solute1d_eachCylechSolAfter.sum(axis=1)

        perirhizalModel.rhizoWAfter = sum(perirhizalModel.rhizoTotCAfter_eachCyl) 
        perirhizalModel.rhizoTotCAfter = sum(perirhizalModel.rhizoTotCAfter_eachCyl) 
        
        # return component content per voxel according to 1d model data 
        # used to compute 3d sink
        perirhizalModel.soil_solute1d_perVoxelAfter = np.array([
            perirhizalModel.soil_solute1d_eachCylechSolAfter[perirhizalModel.getIdCyllMPI(cellId)[0]].sum(axis=0) for cellId in cellIds
        ]).T
        
        assert perirhizalModel.soil_solute1d_perVoxelAfter.shape == (perirhizalModel.numSoluteComp, len(cellIds))

    if rank ==0:
        try: # todo: take out?
            assert min(perirhizalModel.soil_solute1d_perVoxelAfter.flatten()) >=0
        except:
            print("min(perirhizalModel.soil_solute1d_perVoxelAfter)",
                  min(perirhizalModel.soil_solute1d_perVoxelAfter.flatten()),
                  [min(nss) for nss in perirhizalModel.soil_solute1d_perVoxelAfter])
            raise Exception

        # evaluate convergence for water and solute content 1d models
        if len(perirhizalModel.rhizoWAfter_eachCyl_old) == len(perirhizalModel.rhizoWAfter_eachCyl): # not 1st loop
            perirhizalModel.errW1ds = np.linalg.norm(perirhizalModel.rhizoWAfter_eachCyl - perirhizalModel.rhizoWAfter_eachCyl_old)
            perirhizalModel.errC1ds = np.linalg.norm(perirhizalModel.rhizoTotCAfter_eachCyl - perirhizalModel.rhizoTotCAfter_eachCyl_old)
            
        perirhizalModel.rhizoWAfter_eachCyl_old = perirhizalModel.rhizoWAfter_eachCyl.copy()
        perirhizalModel.rhizoTotCAfter_eachCyl_old = perirhizalModel.rhizoTotCAfter_eachCyl.copy()

        if perirhizalModel.debugMode:
            write_file_array("fpit_soil_solute1d_perVoxelAfter1", perirhizalModel.soil_solute1d_perVoxelAfter[1], 
                             directory_ =perirhizalModel.results_dir, fileType = '.txt') 

            write_file_array("fpit_rhizoWAfter_eachCyl", perirhizalModel.rhizoWAfter_eachCyl, 
                                 directory_ =perirhizalModel.results_dir, fileType = '.txt') 
            print("perirhizalModel.rhizoWAfter_eachCyl",perirhizalModel.rhizoWAfter_eachCyl)



def storeOldMassData3d(s,perirhizalModel):
    """
        store mass data of the 3d model at beginning time step
        to evaluate later mole balance error
        @see  massBalanceError3d
    """
    water_content = s.getWaterContent()  # theta per cell [1]


    # 3d soil solute content per solute type and voxel
    perirhizalModel.totC3dBefore_eachVoxeleachComp = s.getTotCContent_each() 

    if rank ==0:
        perirhizalModel.soil_water3dBefore = np.multiply(water_content, s.CellVolumes)  # water per cell [cm3]
        perirhizalModel.totC3dBefore = np.sum(perirhizalModel.totC3dBefore_eachVoxeleachComp)
        assert isinstance(perirhizalModel.totC3dBefore ,float) # just test once: we did the complete sum?
        
        if perirhizalModel.debugMode:
            write_file_array("fpit_soil_water3dBefore", perirhizalModel.soil_water3dBefore, 
                             directory_ =perirhizalModel.results_dir, fileType = '.txt') 
            print("perirhizalModel.soil_water3dBefore",perirhizalModel.soil_water3dBefore)

def storeNewMassData3d(s,perirhizalModel):
    """
        store mass data of the 3d model at end time step
        to evaluate later mole balance error
        @see  massBalanceError3d
    """
    water_content = s.getWaterContent()  # theta per cell [1]
    # 3d soil solute content per solute type and voxel,  (numSoluteComp,numberOfCellsTot)
    perirhizalModel.totC3dAfter_eachVoxeleachComp = s.getTotCContent_each()

    if rank==0:
        perirhizalModel.soil_water3dAfter = np.multiply(water_content, s.CellVolumes)  # water per cell [cm3]

        perirhizalModel.totC3dAfter = np.sum(perirhizalModel.totC3dAfter_eachVoxeleachComp)
        

        # convergence water and solute content 3d models
        perirhizalModel.errW3ds = np.linalg.norm(perirhizalModel.soil_water3dAfter - perirhizalModel.soil_water3dAfter_old)
        perirhizalModel.errC3ds = np.linalg.norm(perirhizalModel.totC3dAfter_eachVoxeleachComp - perirhizalModel.totC3dAfter_eachVoxeleachComp_old)
        perirhizalModel.soil_water3dAfter_old = perirhizalModel.soil_water3dAfter.copy()
        perirhizalModel.totC3dAfter_eachVoxeleachComp_old = perirhizalModel.totC3dAfter_eachVoxeleachComp.copy()

        if perirhizalModel.debugMode:
            write_file_array("fpit_soil_water3dAfter", perirhizalModel.soil_water3dAfter, 
                             directory_ =perirhizalModel.results_dir, fileType = '.txt') 
            print("perirhizalModel.soil_water3dAfter",perirhizalModel.soil_water3dAfter)




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
        self.cell_volumes = s.CellVolumes 
        self.cellIds = perirhizalModel.cellWithRoots # only id of cells with roots
        
        
        # plant-soil solute flow, defined outside of iteration loop
        self.Q_Exud_i = Q_Exud_i
        self.Q_mucil_i = Q_mucil_i
        self.seg_sol_fluxes = Q_Exud_i /sim_time# mol/day for segments
        self.seg_mucil_fluxes = Q_mucil_i/sim_time
        self.sim_time = sim_time
        self.dt = dt
        self.initializeOldData(seg_fluxes, outer_R_bc_wat, outer_R_bc_sol)
        self.err = 1000
        self.errs= 1000
        
    def initializeOldData(self, seg_fluxes, outer_R_bc_wat, outer_R_bc_sol):
        """ create old data object to evaluate convergence """
        # data
        self.rx_old = 0 # plant water potential [cm]
        
        self.seg_fluxes_old = 0 # plant-soil water exchanges
        self.proposed_outer_fluxes_old = 0 # 3d to 1d water flux
        self.proposed_outer_sol_fluxes_old = 0 # 3d to 1d small solute flux (1st dissolved component)
        self.proposed_outer_mucil_fluxes_old = 0 # 3d to 1d large solute (==mucilage) flux (2nd dissolved component)
        self.seg_fluxes = seg_fluxes #.copy()
        self.outer_R_bc_wat_old =  0 # inter-cell water flux in 3d soil model
        self.outer_R_bc_sol_old = 0 # inter-cell solute flux in 3d soil model
        self.outer_R_bc_wat =  outer_R_bc_wat #.copy() # inter-cell water flux in 3d soil model
        self.outer_R_bc_sol = outer_R_bc_sol #.copy() # inter-cell solute flux in 3d soil model
        
        self.thetaCylOld = self.perirhizalModel.getWaterVolumesCyl(doSum = False, reOrder = True, verbose=False)/self.cylVol # cm3/cm3
        
        
        #self.comp1contentOld = perirhizalModel.getContentCyl(idComp=1, doSum = False, reOrder = True) # [mol] small rhizodeposits
        #self.comp2contentOld = perirhizalModel.getContentCyl(idComp=2, doSum = False, reOrder = True) # [mol] mucilage
        
        # matric potential at the segment-exterior interface, i.e. inner values of the (air or soil) cylindric models 
        self.rsx_init  = self.perirhizalModel.get_inner_heads(weather=self.perirhizalModel.weatherX) # store value at beginning time step
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
        self.err2 = 1.e6
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
        if rank ==0:
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
        
        if self.n_iter == 0:
            rhizoWAfter_eachCyl4splitVals = perirhizalModel.getWaterVolumesCyl(doSum = False, reOrder = True) #cm3 water per 1d model
    # each solute for each cyl.
        
        if rank == 0:
            if len(self.emptyCells) > 0:
                outer_R_bc_wat[self.emptyCells] = 0.
                for nc in range(perirhizalModel.numDissolvedSoluteComp):
                    outer_R_bc_sol[nc][self.emptyCells] = 0.# only use BC for cells with at least one root
        
        
            if max(abs(outer_R_bc_wat )) > 0:
                
                try:
                    
                    # thetaCyl_4splitSoilVals: potentially available water in the cylinder
                    proposed_outer_fluxes = perirhizalModel.splitSoilVals(soilVals=outer_R_bc_wat / dt,
                                                        seg_values=self.thetaCyl_4splitSoilVals, # water 'potentially available'
                                                       seg_volume= self.cylVol,dt = dt,
                                                       isWater = True, 
                                                                      verbose = False) #cm3/day
                except:
                    print('distribute3dto1dFlows: outer_R_bc_wat',outer_R_bc_wat.shape)
                    print('thetaCyl_4',self.thetaCyl_4splitSoilVals.shape, self.cylVol.shape,
                          min(self.thetaCyl_4splitSoilVals),max(self.thetaCyl_4splitSoilVals),
                          'self.thetaCylOld',min(self.thetaCylOld),max(self.thetaCylOld),
                         'self.thetaCyl_4splitSoilVals',self.thetaCyl_4splitSoilVals)
                    raise Exception                   
                                                 
            else:
                proposed_outer_fluxes = np.full(self.numSegs, 0.) 
                
            if self.n_iter > 0:    
                # rhizoWAfter_eachCyl4splitVals: value to go from C to C/cm3 water for the solutes
                rhizoWAfter_eachCyl4splitVals = perirhizalModel.rhizoWAfter_eachCyl.copy()
            
                
            if len(perirhizalModel.airSegs) > 0:
                rhizoWAfter_eachCyl4splitVals[perirhizalModel.airSegs]=1. # to avoind division by 0.

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

        if (rank == 0)  and (perirhizalModel.debugMode):
            
            print("outer_R_bc_wat",outer_R_bc_wat)
            #print("self.thetaCyl_4splitSoilVals",self.thetaCyl_4splitSoilVals)
            print("rhizoWAfter_eachCyl4splitVals",rhizoWAfter_eachCyl4splitVals)
            print("proposed_outer_fluxes",proposed_outer_fluxes)
            

        self.proposed_outer_fluxes = proposed_outer_fluxes #comm.bcast(, root = 0)
        self.proposed_outer_sol_fluxes = proposed_outer_sol_fluxes#comm.bcast(, root = 0)
        self.proposed_outer_mucil_fluxes = proposed_outer_mucil_fluxes#comm.bcast(, root = 0)
        
        
    def getCyldatafor3d1dFlow(self, rs_age_i_dt, dt):
        perirhizalModel = self.perirhizalModel

        
        # get data before doing the 'reset' => values at the end of the time step
        self.comp1content = perirhizalModel.getContentCyl(idComp=1, doSum = False, reOrder = True) # [mol] small rhizodeposits
        self.comp2content = perirhizalModel.getContentCyl(idComp=2, doSum = False, reOrder = True) # [mol] mucilage
    
        
        if rank ==0:
            # when there is not transpiration, use data at the beginning of the time step for water flow
            if (perirhizalModel.beforeAtNight and PhloemPhotosynthesis.noTranspiration(perirhizalModel, rs_age_i_dt, dt)):
                self.thetaCyl_4splitSoilVals = self.thetaCylOld
            else:
                # get data before doing the 'reset' => values at the end of the time step
                # self.thetaCyl = perirhizalModel.getWaterVolumesCyl(doSum = False, reOrder = True)/self.cylVol # cm3 water
                # get all potentially available water == water after reset + water taken up by plant
                seg_fluxes_root  = self.seg_fluxes
                if len(perirhizalModel.airSegs) > 0:
                    seg_fluxes_root[perirhizalModel.airSegs]  = 0.
                self.thetaCyl_4splitSoilVals = (self.thetaCylOld * self.cylVol + (seg_fluxes_root + perirhizalModel.flow1d1d_w)* dt )/self.cylVol
                #because of unlimited flow (seg_fluxes_root + perirhizalModel.flow1d1d_w), might get values out of the [theta_r, theta_s] bounds
                self.thetaCyl_4splitSoilVals = np.maximum(np.minimum(
                                                        self.thetaCyl_4splitSoilVals, 
                                                          perirhizalModel.vg_soil.theta_S),
                                                         perirhizalModel.theta_wilting_point)
                if len(perirhizalModel.airSegs) > 0:
                    self.thetaCyl_4splitSoilVals[perirhizalModel.airSegs]  = 0.
                    

            if (len(perirhizalModel.airSegs) > 0):
                assert (self.comp1content[perirhizalModel.airSegs] == 0.).all() #no solutes in the space around the air segments
                assert (self.comp2content[perirhizalModel.airSegs] == 0.).all() #no solutes in the space around the air segments
            assert len(self.comp1content) == len(perirhizalModel.eidx_all)    

            # assert (self.thetaCyl_4splitSoilVals >= 0.).all() # because of seg_fluxes, could have thetaCyl_4splitSoilVals < 0
            assert (self.comp1content.astype(np.float64).flatten().min() >= 0.)
            assert ( self.comp2content.astype(np.float64).flatten().min()>= 0.)

            # self.comp1content = perirhizalModel.soil_solute1d_eachCylechSolAfter[:,0]
            
            # [mol] small rhizodeposits
            # self.comp2content = perirhizalModel.soil_solute1d_eachCylechSolAfter[:,1]  # [mol] mucilage

            if perirhizalModel.debugMode:
                print("self.cylVol",self.cylVol)
                print("seg_fluxes_root",seg_fluxes_root)
                print("flow1d1d_w",perirhizalModel.flow1d1d_w)
                write_file_array("fpit_thetaCyl_4splitSoilVals", self.thetaCyl_4splitSoilVals, 
                             directory_ =perirhizalModel.results_dir, fileType = '.txt') 
            
        
        
    def do1d1dFlow(self):
        """
            wrapper function for computation of inter-cylinder flow
            use data after doing the reset, otherwise get diverging results
        """
        perirhizalModel = self.perirhizalModel
        soilVals_ = np.array(self.s.getSolutionHead()) #<= data beginning time step
        rhizoWAfter_eachCyl_divide = perirhizalModel.rhizoWAfter_eachCyl.copy()
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

        if rank ==0:
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
        perirhizalModel = self.perirhizalModel
        cellIds = self.cellIds
        if rank == 0:
            # get changes in colute content not caused by flow [mol/day]
            # todo: get it directly from dumux
            sources_sol_from1d = np.full( (self.perirhizalModel.numSoluteComp, self.s.numberOfCellsTot),0. )
            for nc in range(self.perirhizalModel.numSoluteComp):
                sources_sol_from1d[nc][cellIds] = np.array(
                    perirhizalModel.soil_solute1d_perVoxelAfter[nc] - perirhizalModel.soil_solute1d_perVoxelBefore[nc] - self.outer_R_bc_sol[nc][cellIds]
                )/dt
            #sources_sol_from1d = sources_sol_from1d  
            
            assert sources_sol_from1d.shape == (self.perirhizalModel.numSoluteComp, self.s.numberOfCellsTot)
            
            # store error
            self.s.errSoil_source_sol_abs = sum(sources_sol_from1d.flatten()) - (sum(self.Q_Exud_i) + sum(self.Q_mucil_i))/self.sim_time

            if perirhizalModel.debugMode:
                print("self.s.errSoil_source_sol_abs",self.s.errSoil_source_sol_abs,"(sum(self.Q_Exud_i) + sum(self.Q_mucil_i))",
                      (sum(self.Q_Exud_i) + sum(self.Q_mucil_i))," sum(sources_sol_from1d.flatten())*dt", sum(sources_sol_from1d.flatten())*dt)
            
            if (sum(self.Q_Exud_i) + sum(self.Q_mucil_i))/self.sim_time != 0.:
                s.errSoil_source_sol_rel = abs(s.errSoil_source_sol_abs/((sum(self.Q_Exud_i) + sum(self.Q_mucil_i))/self.sim_time)*100)
            else:
                s.errSoil_source_sol_rel = np.nan
            self.sources_sol_from1d = sources_sol_from1d

            
    def compute1dChangesWater(self, dt):
        if (rank == 0):
            plantModel = self.plantModel
            perirhizalModel = self.perirhizalModel
            s = self.s
            net_PWU_ = plantModel.sumSegFluxes(self.seg_fluxes)  #  plant water uptake [cm3/day]  per soil cell
            net_PWU = np.zeros(s.numberOfCellsTot)
            #easier to handle array than dict. maybe change sumSegFluxes() function to choose type of output
            net_PWU[np.array(list(net_PWU_.keys()))] = np.array(list(net_PWU_.values())) 


            self.net_PWU = net_PWU

            net_PWU_limited_ = plantModel.sumSegFluxes(self.seg_fluxes_limited)  # [cm3/day]  per soil cell
            net_PWU_limited = np.zeros(s.numberOfCellsTot)
            #easier to handle array than dict. maybe change sumSegFluxes() function to choose type of output
            net_PWU_limited[np.array(list(net_PWU_limited_.keys()))] = np.array(list(net_PWU_limited_.values()))


            self.net_PWU_limited = net_PWU_limited
            #print('compute1dChangesWater', 'seg_', sum(self.seg_fluxes), sum(net_PWU),sum(self.seg_fluxes)- sum(net_PWU),
            #                    'seg_fluxes_limited', sum(self.seg_fluxes_limited), sum(net_PWU_limited) ,
            #                    sum(self.seg_fluxes_limited)- sum(net_PWU_limited) ,
            #                  'diff',  sum(net_PWU_limited)- sum(net_PWU),sum(self.seg_fluxes)- sum(self.seg_fluxes_limited))
            perirhizalModel.SinkLim3DS =max( abs((self.net_PWU_limited - self.net_PWU)/ 
                                                np.where(self.net_PWU,self.net_PWU,1.))*100.) # at the end of the fixed point iteration, should be ~ 0 

            if perirhizalModel.debugMode:
                write_file_array("fpit_net_PWU", net_PWU, 
                                     directory_ =perirhizalModel.results_dir, fileType = '.txt') 
                write_file_array("fpit_net_PWU_limited", net_PWU_limited, 
                                     directory_ =perirhizalModel.results_dir, fileType = '.txt') 
                                
        
    def compute3DSource(self, dt):
        """
            compute source/sink of water and solutes in the 3d soil according to the changes in the 1d models
        """
        perirhizalModel = self.perirhizalModel
        results_dir = perirhizalModel.results_dir
        s = self.s
        
        self.compute1dChangesSolute(dt)
        self.compute1dChangesWater(dt)
        
        if rank == 0:
            # source/sink in the 3d soil according to changes in the 1d models [(cm3 water or mol C)/day]
            soil_sources_limited = np.concatenate((np.array([self.net_PWU_limited]),self.sources_sol_from1d ))
        

        # when setting source: we limit according to the content at the beginning 
        # of the time step + max(net incoming flow, 0.) 
        # to get the maximum potential value of content that can be
        # in the voxel if needed
        
        for idComp in range(s.numComp):#cm3/day, mol/day
            if rank == 0:
                SSL = soil_sources_limited[idComp].copy()
                
                # compute maximum poentially available content for sink (cm3 water or mol)
                if idComp == 0:
                    maxPotentialAvailable = (np.array(perirhizalModel.soil_water3dBefore) +
                                         np.maximum(self.outer_R_bc_wat, self.outer_R_bc_wat*0 ))
                else:
                    maxPotentialAvailable = (perirhizalModel.totC3dBefore_eachVoxeleachComp[idComp-1] +
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
                        
                if idComp < 2 and perirhizalModel.debugMode:
                    write_file_array("setsourceLim3_" + str(idComp), SSL,
                                     directory_=results_dir, fileType=".csv")
            else:
                SSL = None
                                 
            self.sendSource2dumux(SSL, idComp)
            
           
    def sendSource2dumux(self, SSL, idComp):
        s = self.s
        perirhizalModel=self.perirhizalModel
        results_dir = perirhizalModel.results_dir
        
        SSL = comm.bcast(SSL, root=0)
        
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

        # todo: probably to adapt
        assert dt >= 10/(24*3600)
        s.ddt = min(max(s.ddt, 1/(24*3600)), dt/10.)
        #s.maxDt = max(min(s.maxDt, dt/4.), s.ddt)
        assert dt > s.ddt
        assert s.maxDt > s.ddt
        #print('troubleshoot data',
        #      'time',s.ddt, s.maxDt, dt)
        #print('params',s.getParameter("Newton.MaxRelativeShift"),
        #      s.getParameter("Newton.MaxSteps"))
        while redoSolve:
            #s.ddt =min( 1.e-5,s.ddt)#or just reset to 1e-5?
            #s.setMaxTimeStepSize(s.maxDt) # reset each time
            
            try:
                decreaseMaxRelShift = False
                if rank==0:
                    print("solve 3d soil")
                #s.solve(dt)  # in modules/solverbase.py

                helpfull.run_with_timeout(60.*5,s.solve,dt) # time out error after Xmn
                if rank==0:
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

                s.setParameter("Newton.MaxSteps", "100")
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
           
    def computeConvergence2(self):
        perirhizalModel = self.perirhizalModel
        s = self.s
        plantModel = self.plantModel
        
        # convergence wat. pot. at root-soil interface
        rsx = perirhizalModel.get_inner_heads(weather=perirhizalModel.weatherX)

        if rank == 0:
            # convergence plant wat. pot
            rx = np.array(plantModel.psiXyl)
            rx_divide = np.where(rx != 0, rx, 1.)
            self.errRxPlant = max(abs((rx - self.rx_old) / rx_divide) * 100.)
            self.rx_old = rx.copy()

            
            rsx_divide = np.where(abs(rsx) < abs(self.rsx_input),
                                  rsx, self.rsx_input)
            perirhizalModel.errWrsiRealInputs = abs(
                (rsx - self.rsx_input) / rsx_divide) * 100.
            perirhizalModel.errWrsiRealInputs[abs(rsx - self.rsx_input) <= 5] = 0.
            perirhizalModel.errWrsiRealInput = max(perirhizalModel.errWrsiRealInputs)
            
            rsx_divide = np.where(rsx != 0., rsx, 1.)
            perirhizalModel.errWrsis = abs((rsx - self.rsx_old) / rsx_divide) * 100.  # max()
            perirhizalModel.errWrsis[abs(rsx - self.rsx_old) <= 5] = 0.
            perirhizalModel.errWrsi =  max(perirhizalModel.errWrsis)
            self.rsx_old = rsx;
            self.rsx_olds.append(self.rsx_old)


            # convergence BC flow 1d models
            diffBCS1dsFluxIn = np.array(
                self.seg_fluxes) - self.seg_fluxes_old  # only for water as plant exud is outside of loop

            self.errBCS1dsFluxIn =max(abs((
                diffBCS1dsFluxIn/ np.where(np.array(self.seg_fluxes),
                                           np.array(self.seg_fluxes),1.))*100))


            perirhizalModel.err2 = max(self.errRxPlant,perirhizalModel.errW1ds, perirhizalModel.errC1ds,
                                                 perirhizalModel.errWrsi,perirhizalModel.errWrsiRealInput,
                                   perirhizalModel.rhizoMassWError_rel,
                                   perirhizalModel.rhizoMassCError_rel)
                            

    def computeConvergence(self):
        """
            long funciton, compute different types of errors:
            non-convergence, mass balance errors, prescribed vs realised BCs
        """
        perirhizalModel = self.perirhizalModel
        s = self.s
        

        perirhizalModel.solve_gave_up = (np.array(comm.gather(perirhizalModel.solve_gave_up ,root = 0))).any()
        
        if rank ==0:
            self.diffBCS1dsFluxOut =   self.proposed_outer_fluxes  - self.proposed_outer_fluxes_old 
            
                             
            diffBCS1dsFluxOut_sol =   self.proposed_outer_sol_fluxes  - self.proposed_outer_sol_fluxes_old 
            diffBCS1dsFluxOut_mucil =   self.proposed_outer_mucil_fluxes  - self.proposed_outer_mucil_fluxes_old

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

            perirhizalModel.maxDiff1d3dCW_abs =np.array( perirhizalModel.maxDiff1d3dCW_abs)

            # only look at the relative error if the absolute error is high enough
            compErrorAboveLim = np.where(perirhizalModel.sumDiff1d3dCW_abs > 1e-13) 

            # to not depend on cumulative error
            perirhizalModel.diff1d3dCurrant = max(np.append((perirhizalModel.sumDiff1d3dCW_abs - perirhizalModel.sumDiff1d3dCW_absOld)[compErrorAboveLim],0.)) 
            perirhizalModel.diff1d3dCurrant_rel =max(np.append((perirhizalModel.sumDiff1d3dCW_rel - perirhizalModel.sumDiff1d3dCW_relOld)[compErrorAboveLim],0.))

            perirhizalModel.maxdiff1d3dCurrant = max(np.append((perirhizalModel.maxDiff1d3dCW_abs - perirhizalModel.maxDiff1d3dCW_absOld)[compErrorAboveLim],0.)) 
            perirhizalModel.maxdiff1d3dCurrant_rel =max(np.append((perirhizalModel.maxDiff1d3dCW_rel - perirhizalModel.maxDiff1d3dCW_relOld)[compErrorAboveLim],0.))

            
            # one metric to decide if we stay in the iteration loop or not
            
            perirhizalModel.err = max(perirhizalModel.err2,
                                    s.bulkMassErrorWater_rel, 
                                   s.bulkMassCErrorPlant_rel, s.bulkMassCError1ds_rel
                            )
            # one array to do printing
            perirhizalModel.errs =np.array([
                            # non-convergence metrics
                            self.errRxPlant, perirhizalModel.errW1ds, perirhizalModel.errW3ds,perirhizalModel.errC1ds,
                perirhizalModel.errC3ds, perirhizalModel.errWrsi,
                perirhizalModel.errWrsiRealInput,
                            self.errBCS1dsFluxIn, errBCS1dsFluxOut,errBCS1dsFluxOut_sol,
                errBCS1dsFluxOut_mucil,
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


    def massBalanceError1d(self,dt):
        """
            get mass balance error of 1d models
        """

        if rank == 0:
            perirhizalModel = self.perirhizalModel
            rhizoSegsId = perirhizalModel.rhizoSegsId # plant segments with a rhizosphere model
            airSegsId = self.airSegsId

            ############ solutes
            # get error according to the 'limited' (== realised ) boundary fluxes
            # should be always ~ 0 
            errorsEachC = perirhizalModel.rhizoTotCAfter_eachCyl - ( 
                perirhizalModel.rhizoTotCBefore_eachCyl + (self.seg_fluxes_limited_sol_In+ self.seg_fluxes_limited_mucil_In+\
                                                                  self.seg_fluxes_limited_sol_Out+ self.seg_fluxes_limited_mucil_Out)*dt)
            # store absolute total error for limited flow
            perirhizalModel.rhizoMassCError_absLim = sum(abs(errorsEachC[rhizoSegsId]))


            # get error according to the proposed (==prescribed) flux
            # need to be ~ 0 when leaving fixed point iteration
            if not perirhizalModel.doSoluteUptake:
                errorsEachC = perirhizalModel.rhizoTotCAfter_eachCyl - ( perirhizalModel.rhizoTotCBefore_eachCyl +
                                                         (self.seg_sol_fluxes+ 
                                                          self.proposed_outer_sol_fluxes+ 
                                                          perirhizalModel.flow1d1d_sol +
                                                          self.seg_mucil_fluxes+ 
                                                          self.proposed_outer_mucil_fluxes+ 
                                                          perirhizalModel.flow1d1d_mucil )*dt)

            perirhizalModel.rhizoMassCError_abs  = sum(abs(errorsEachC[rhizoSegsId]))

            # store relative total error 
            if sum(perirhizalModel.rhizoTotCAfter_eachCyl) != 0:
                perirhizalModel.rhizoMassCError_relLim = abs(perirhizalModel.rhizoMassCError_absLim/sum(perirhizalModel.rhizoTotCAfter_eachCyl)*100)
                perirhizalModel.rhizoMassCError_rel = abs(perirhizalModel.rhizoMassCError_abs/sum(perirhizalModel.rhizoTotCAfter_eachCyl)*100)
            else:
                perirhizalModel.rhizoMassCError_relLim = np.nan
                perirhizalModel.rhizoMassCError_rel = np.nan



            ############ water
            # get error according to the 'limited' (== realised ) boundary fluxes
            # should be always ~ 0 
            errorsEachW = perirhizalModel.rhizoWAfter_eachCyl - ( 
                perirhizalModel.rhizoWBefore_eachCyl + (self.seg_fluxes_limited + self.seg_fluxes_limited_Out)*dt)


            # store absolute total error for limited flow
            perirhizalModel.rhizoMassWError_absLim = sum(abs(errorsEachW[rhizoSegsId]))
            perirhizalModel.errorsEachWLim = errorsEachW

            # get error according to the proposed (==prescribed) flux
            # need to be ~ 0 when leaving fixed point iteration 
            perirhizalModel.errorsEachW = perirhizalModel.rhizoWAfter_eachCyl - ( perirhizalModel.rhizoWBefore_eachCyl + (self.seg_fluxes+ self.proposed_outer_fluxes+ perirhizalModel.flow1d1d_w)*dt)
            perirhizalModel.rhizoMassWError_abs  = sum(abs(perirhizalModel.errorsEachW[rhizoSegsId]))

            # store relative total error 
            perirhizalModel.rhizoMassWError_relLim = abs(perirhizalModel.rhizoMassWError_absLim/sum(perirhizalModel.rhizoWAfter_eachCyl)*100)
            perirhizalModel.rhizoMassWError_rel = abs(perirhizalModel.rhizoMassWError_abs/sum(perirhizalModel.rhizoWAfter_eachCyl)*100)


            print(f'\t\trelative error balance soil 1d (%)?\n\t\t\t\tfrom PWU: {perirhizalModel.rhizoMassWError_rel:.2e}, from PWU-limited: {perirhizalModel.rhizoMassWError_relLim:.2e}, from PCU: {perirhizalModel.rhizoMassCError_rel:.2e}, from PCU-limited: {perirhizalModel.rhizoMassCError_relLim:.2e}'
                )



    def massBalanceError3d(self,dt):
        if rank == 0:
            s = self.s
            perirhizalModel = self.perirhizalModel
            # according to plant data # I need to add the water flux to be sure
            s.bulkMassErrorWater_absLim = sum(abs(perirhizalModel.soil_water3dAfter - perirhizalModel.soil_water3dBefore  -self.sources_wat_from3d-self.outer_R_bc_wat))

            s.bulkMassErrorWater_abs = sum(abs(perirhizalModel.soil_water3dAfter -  perirhizalModel.soil_water3dBefore-self.net_PWU*dt-self.outer_R_bc_wat))
            
            
            assert s.bulkMassErrorWater_abs - s.bulkMassErrorWater_absLim >= -1e-13
            
            s.bulkMassErrorWater_relLim = abs(s.bulkMassErrorWater_absLim /sum(perirhizalModel.soil_water3dAfter) )*100
            s.bulkMassErrorWater_rel = abs(s.bulkMassErrorWater_abs /sum(perirhizalModel.soil_water3dAfter) )*100
            s.bulkMassErrorWaterSink_abs = self.net_PWU_limited*dt- self.sources_wat_from3d # sink sent to dumux == sink implemented?


            s.bulkMassCErrorPlant_abs = abs(perirhizalModel.totC3dAfter - ( perirhizalModel.totC3dBefore + sum(self.Q_Exud_i) + sum(self.Q_mucil_i)))
            if (perirhizalModel.totC3dAfter > 0) and (not perirhizalModel.doSoluteUptake):
                s.bulkMassCErrorPlant_rel = abs(s.bulkMassCErrorPlant_abs/perirhizalModel.totC3dAfter*100)
            else:
                s.bulkMassCErrorPlant_rel = np.nan

            # according to 1d data
            s.bulkMassCError1ds_abs = abs(perirhizalModel.totC3dAfter - ( perirhizalModel.totC3dBefore + sum(self.sources_sol_from1d.flatten())*dt))
            if perirhizalModel.totC3dAfter > 0:
                s.bulkMassCError1ds_rel = abs(s.bulkMassCError1ds_abs/perirhizalModel.totC3dAfter*100)
            else:
                s.bulkMassCError1ds_rel = np.nan



            print(f"relative error balance soil 3d (%)?\n\t\tfrom PWU {s.bulkMassErrorWater_rel:.2e}, from PWU-limited {s.bulkMassErrorWater_relLim:.2e},from PCU {s.bulkMassCErrorPlant_rel:.2e}, from 1ds C-change {s.bulkMassCError1ds_rel:.2e}"
                )
            if perirhizalModel.debugMode:
                print("s.bulkMassCErrorPlant_abs",s.bulkMassCErrorPlant_abs," s.bulkMassCError1ds_abs", s.bulkMassCError1ds_abs,"sum(self.Q_Exud_i) + sum(self.Q_mucil_i)",sum(self.Q_Exud_i) + sum(self.Q_mucil_i)," sum(self.sources_sol_from1d.flatten())*dt)", sum(self.sources_sol_from1d.flatten())*dt)

            ### for each voxel
            s.bulkMassErrorWaterAll_abs = abs(perirhizalModel.soil_water3dAfter - (perirhizalModel.soil_water3dBefore  + self.sources_wat_from3d + self.outer_R_bc_wat))
            s.bulkMassErrorWaterAll_rel = abs(s.bulkMassErrorWaterAll_abs /perirhizalModel.soil_water3dAfter )*100



            ## inner mass balance error 

            s.bulkMassCError1dsAll_abs = abs(perirhizalModel.totC3dAfter_eachVoxeleachComp - (perirhizalModel.totC3dBefore_eachVoxeleachComp + 
                                                        self.sources_sol_from3d + self.outer_R_bc_sol))


            totC3dAfter_eachVoxeleachCompTemp = np.where(perirhizalModel.totC3dAfter_eachVoxeleachComp != 0,
                                                      perirhizalModel.totC3dAfter_eachVoxeleachComp,
                                                      1)
            s.bulkMassCError1dsAll_rel = abs(s.bulkMassCError1dsAll_abs*100/totC3dAfter_eachVoxeleachCompTemp)

            try:
                assert s.bulkMassCError1dsAll_rel.shape == (s.numSoluteComp, s.numberOfCellsTot )
            except:
                print(s.bulkMassCError1dsAll_rel.shape , 
                      (s.numSoluteComp, s.numberOfCellsTot ))
                raise Exception

