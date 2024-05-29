
import plantbox as pb
import functional.xylem_flux as xylem_flux
import sys
from functional.xylem_flux import XylemFluxPython
from rosi_richards10c_cyl import RichardsNCCylFoam  # C++ part (Dumux binding)
from richards_no_mpi import RichardsNoMPIWrapper  # Python part of cylindrcial model (a single cylindrical model is not allowed to run in parallel)
from fv.fv_grid import *
import fv.fv_richards as rich  # local pure Python cylindrical models
import functional.van_genuchten as vg
from scenario_setup import write_file_array, setDefault, write_file_float

import numpy as np
import matplotlib.pyplot as plt
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank(); max_rank = comm.Get_size(); size = comm.Get_size()
import psutil
from air_modelsPlant import AirSegment
from scipy import sparse
import scipy.sparse.linalg as LA

from scipy.interpolate import PchipInterpolator,  CubicSpline
from helpfull import StdoutRedirector

class RhizoMappedSegments(pb.MappedPlant):#XylemFluxPython):#
    """
        Adds 1-dimensional rhizospere models to each root segment of a MappedSegments (or later MappedPlant)    
        
        modes:        
        "dumux_w"             
        "dumux_3c"        
        "dumux_10c"        
    """

    # TODO copy mapped segments constructors (!)...

    def __init__(self,  soilModel,
                usemoles, seedNum=None,limErr1d3dAbs = 1e-11):
        """ @param file_name is either a pb.MappedRootSystem, pb.MappedSegments, or a string containing a rsml filename"""
        if not seedNum is None:
            super().__init__(seednum = seedNum)
        else:
            super().__init__()
            
        # simulation characteristics            
        self.results_dir = soilModel.results_dir # directory for output data
        self.mode = "dumux_10c" # mode of the simulation represented
        self.useMoles = usemoles
        self.numFluidComp = soilModel.numFluidComp
        self.numDissolvedSoluteComp = soilModel.numDissolvedSoluteComp
        self.numSoluteComp = soilModel.numSoluteComp
        self.numComp = soilModel.numComp
        self.debugMode = True
        
        
        # constants
        self.molarMassWat = soilModel.molarMassWat # [g/mol]
        self.densityWat_m3 =soilModel.densityWat_m3 #[g/m3]
        self.molarDensityWat_m3 =  self.densityWat_m3 /self.molarMassWat # [mol/m3] = [g/m3] /  [g/mol]  
        
        # changes with growing plant (to update)
        # plant shape        
        self.wilting_point = soilModel.wilting_point
        self.theta_wilting_point = vg.water_content( self.wilting_point, soilModel.vg_soil)
        self.seg2cell_old = {}
        self.cell2seg_old = {}
        self.seg_length = np.array([]) 
        self.seg_length_old = np.array([])
        
        # 3d soil 
        self.soilModel = soilModel
        self.soil = self.soilModel.soil
        self.vg_soil = self.soilModel.vg_soil
        self.cellWithRoots = [] # id of the (3d) bulk soil with at least one root
        self.sizeSoilCell = comm.bcast(self.soilModel.getCellVolumes(), root = 0) #cm3
        self.solidDensity = self.soilModel.solidDensity#2700 # [kg/m^3 solid]
        self.solidMolarMass = self.soilModel.solidMolarMass#60.08e-3 # [kg/mol]           
        self.solidMolDensity = self.solidDensity/self.solidMolarMass # [mol / m3 solid] =[kg/m^3 solid] / [kg/mol] 
        self.bulkDensity_m3 = self.solidMolDensity*(1.- self.soilModel.vg_soil.theta_S) #porosity == theta_s # [mol / m3 scv] = [mol / m3 solid] * [m3 solid /m3 space]
        
        
        # all 1d models
        self.eidx_all = [] # all 1d segs id 
        self.eidx = np.array([], dtype = np.int64) # id of 1d segments on that thread
        self.cyls = [] # local cylinders
        self.outer_radii = None # outer radius of the 1d models (cm)
        self.rhizoVol = None # volume of 1d models (cm3)
        self.dx2 = [] # distance between inner cell center and surface of plant segment (cm)
        # 1d soil models
        self.NC = 10 # number of faces of the 1d meshes, == dof + 1
        self.logbase = 0.5 # logbase for scaling of the 1d mesh cells
        self.cylSoilidx_all = [] # all soil perirhizal cylinders id ( = self.eidx_all - airsegs_all)
        self.cylSoilidx = [] # perirhizal cylinders id on local thread without airSegments ( = self.eidx - airsegs)
        self.cylsSoil = [] # 1d soil models carried by thread
        self.nsSoilOld = 0 # total number of 1d soil segs at the last time step
        self.toAddSoil = np.array([]) # number of 1d soil segs  to add for each threads
        self.repartitionSoilOld = np.array([0 for i in range( max_rank)]) # old division of the 1d soil models between the threads
        # 1d air domains == dummy perirhizal models
        self.airSegs = np.array([])# all air perirhizal cylinders id
        self.nsAirOld = 0 # total number of 1d air segs at the last time step
        self.toAddAir = np.array([]) # number of 1d air segs  to add for each threads
        self.repartitionAirOld = np.array([0 for i in range( max_rank)])  # old division of the 1d air models between the threads
        
        # error rates and fails
        self.limErr1d3dAbs = limErr1d3dAbs
        self.maxDiff1d3dCW_absOld = np.full(self.numComp, 0.)
        self.maxDiff1d3dCW_relOld = np.full(self.numComp, 0.)
        self.maxdiff1d3dCurrant_rel = np.Inf
        self.maxdiff1d3dCurrant = np.Inf
        self.sumDiff1d3dCW_absOld = np.full(self.numComp, 0.)
        self.sumDiff1d3dCW_relOld = np.full(self.numComp, 0.)
        self.diff1d3dCurrant_rel = 0.
        self.diff1d3dCurrant = np.Inf
        self.solve_gave_up = False
        
        
        # 1d1d flow
        self.do1d1dFlow = True
        self.flow1d1d_w = np.zeros(1)
        self.flow1d1d_sol = np.zeros(1)
        self.flow1d1d_mucil = np.zeros(1)
        

    def set_phloem_flux(self, plantModel):
        """ save the plant model"""
        self.plantModel = plantModel
        
    def reset(self):
        for cyl in self.cyls:
            cyl.reset()
        
    def save(self):
        for cyl in self.cyls:
            cyl.save()
            
    def resetManual(self):
        for cyl in self.cyls:
            cyl.resetManual()
            
    def saveManual(self):
        for cyl in self.cyls:
            cyl.saveManual()
            
    
    def getVolumes(self,vertices, length):
        """ volume of each cell of a hollow cylinder mesh """
        return   np.pi * (vertices[1:] ** 2 - vertices[:-1]**2) * length
    
    def phaseDensity(self, compId):# mol / m3
        return self.soilModel.phaseDensity( isDissolved = (compId <= self.numDissolvedSoluteComp))
        
    def getCellIds(self):# only send cells which have roots in them
        return self.cellWithRoots
            
    def getDistribution(self, ns:int):      
        """ divide ns segments between the threads """
        dcyl = int(np.floor(ns /  max_rank))
        repartition = np.array([dcyl for i in range( max_rank)])
        residual = ns - max_rank * dcyl 
        repartition[:residual] += 1  
        return repartition
        
    def getAirSegsId(self):
        """ get id of abovegronud root and shoot segments """
        aboveGround = np.array([])
        if not (self.cell2seg.get(-1) is None):
            aboveGround = self.cell2seg.get(-1)                
        self.airSegs = np.array(list(set(np.concatenate((aboveGround,
                                                    np.where(np.array(self.organTypes) != 2)[0])) )))
        self.airSegs.sort()#to go from local segIdx to global segIdx easely
    
    def getNewCyldistribution(self):
        """ get new ids and shape of 1d models and define distribute between the threads """
        self.outer_radii = np.array(self.segOuterRadii(type = 0)) 
        
        
        self.getAirSegsId()
        if len (self.airSegs) > 0:
            self.outer_radii[self.airSegs ] = np.array(self.radii)[self.airSegs ]*1.1 #dummy outer radii                
        
        ## check value outer radiis        
        if self.debugMode :
            if isinstance(self.outer_radii_old, type(self.outer_radii)):
                assert (self.outer_radii_old >= self.outer_radii[:len(self.outer_radii_old)]).all()
            
        self.seg_length = np.array(self.segLength())#node might have shifted: new length for pre-existing segments          
        repartitionSoil = self.getDistribution(len(self.seg_length) - len(self.airSegs )) # distribute soil segments between threads
        self.toAddSoil = repartitionSoil - self.repartitionSoilOld # see how many new segments in each threads
        
        repartitionAir = self.getDistribution(len(self.airSegs ))
        self.toAddAir = repartitionAir - self.repartitionAirOld
        
            
        self.repartitionSoilOld = repartitionSoil
        self.repartitionAirOld = repartitionAir
        
        if self.debugMode :
            assert (min(self.toAddSoil) >= 0.) and (min(self.toAddAir) >=0.) 
            assert sum(repartitionSoil) == len(self.seg_length) - len(self.airSegs )
            assert sum(self.toAddSoil) == (len(self.seg_length) - len(self.airSegs ) - self.nsSoilOld )
            assert sum(repartitionAir) ==  len(self.airSegs )
            assert sum(self.toAddAir) == ( len(self.airSegs ) - self.nsAirOld )
    
    
    
                
    def checkVolumeAndRadii(self, finishedUpdate):
        """ make sure that 1d cyl volume == voxel volume for each voxel"""
        cellIds = self.getCellIds()
        for cellId in cellIds:
            vol_total = self.sizeSoilCell[cellId]  # cm3.
            
            idCylsAll, idCyls =  self.getIdCyllMPI(cellId)
            if(len(idCylsAll) >0):
            
                # from dumux volume
                vol_rhizo =self.getVolumesCyl(idCyls,doSum = True)
                assert (abs((vol_total - vol_rhizo)) <= 1e-10)
                
            if not finishedUpdate:
                idSegs = np.array(self.cell2seg[cellId])#all segments in the cell
                idCylsAll = np.array([ids for ids in idSegs if self.organTypes[ids]==2 ]) # 
            if(len(idCylsAll) >0):
                # from inner and outer radii
                lengths_I = np.array(self.seg_length)[idCylsAll]
                radii_in = np.array(self.radii)[idCylsAll]
                radii_out = np.array(self.outer_radii)[idCylsAll]
                vol_rootRhizo = radii_out * radii_out * np.pi * lengths_I
                vol_root = radii_in * radii_in * np.pi * lengths_I 
                try:
                    assert abs(((vol_total - sum(vol_rootRhizo - vol_root))/vol_total)*100) < 1e-12
                except:
                    print("checkVolumeAndRadii ", idCylsAll,
                            cellId,vol_total , vol_rootRhizo,vol_root,vol_total - sum(vol_rootRhizo - vol_root) )
                    print("self.getVolumesCyl(idCyls,doSum = True)",self.getVolumesCyl(idCyls,doSum = False))
                    print('self.rhizoVol',self.rhizoVol[idCylsAll])
                    raise Exception
    
                    
    def check1d3dDiff(self, 
                        diff1d3dCW_abs_lim = np.Inf, # maximal axcepted absolute arror
                              verbose_ = False,
                              diff1d3dCW_rel_lim =np.full(10, np.Inf) # maximal accepted relative error
                              ):
        """ check difference between soil content according to 1d and 3d soil """
        
        if diff1d3dCW_abs_lim is None:
            diff1d3dCW_abs_lim = self.limErr1d3dAbs
        
        
        cellIds = self.getCellIds()
        self.allDiff1d3dCW_rel = np.full((self.numComp,max(cellIds)+1), 0.)
        self.allDiff1d3dCW_abs = np.full((self.numComp,max(cellIds)+1), 0.)
        self.contentIn3d = np.full((self.numComp,max(cellIds)+1), 0.)
        self.contentIn1d = np.full((self.numComp,max(cellIds)+1), 0.)
        
        self.sumDiff1d3dCW_abs = np.full(self.numComp, 0.)
        self.sumDiff1d3dCW_rel = np.full(self.numComp, 0.)
        self.maxDiff1d3dCW_abs = np.full(self.numComp, 0.)
        self.maxDiff1d3dCW_rel = np.full(self.numComp, 0.)
        
        
        wat_total_ = comm.bcast( self.soilModel.getWaterVolumes(), root=0)# [cm3].
        
        self.all1d3dDiff =  np.full(len(wat_total_), 0.)

        mol_total_ =comm.bcast( [self.soilModel.getContent(idComp, idComp <= 2)  for  idComp in range(1,self.numComp) ], root=0)# solute content [mol].
        
        
        for cellId in cellIds:
        
            wat_total = wat_total_[cellId] # [cm3].
            
            idCylsAll, idCyls = self.getIdCyllMPI(cellId)
            
            if len(idCylsAll)> 0:
                localIdCyls = self.getLocalIdCyls(idCyls)
                wat_rhizo = self.getWaterVolumesCyl(idCyls, doSum = True) #cm3 
                
                
                if rank == 0:
                    # water
                    diff1d3dCW_abs = abs(wat_total - wat_rhizo)
                    self.all1d3dDiff[cellId] = diff1d3dCW_abs
                    
                    diff1d3dCW_rel = abs(diff1d3dCW_abs/(wat_total ) *100)
                    self.allDiff1d3dCW_abs[0][cellId] = diff1d3dCW_abs
                    self.allDiff1d3dCW_rel[0][cellId] = diff1d3dCW_rel
                    self.sumDiff1d3dCW_abs[0] += diff1d3dCW_abs
                    
                    self.maxDiff1d3dCW_abs[0] = max(self.maxDiff1d3dCW_abs[0], diff1d3dCW_abs)
                    if diff1d3dCW_abs > 1e-13:
                        self.maxDiff1d3dCW_rel[0] = max(self.maxDiff1d3dCW_rel[0], diff1d3dCW_rel)
                        
                    if (diff1d3dCW_abs > diff1d3dCW_abs_lim) and (diff1d3dCW_rel > diff1d3dCW_rel_lim[0]): 
                        print("check1d3dDifferror")
                        print(cellId, 0, 'wat_total',wat_total,'wat_rhizo', wat_rhizo, )
                        print("wat_rhizo_",wat_rhizo)
                        print("Diff1d3dW",diff1d3dCW_abs, diff1d3dCW_rel, self.sumDiff1d3dCW_abs[0],
                             self.maxDiff1d3dCW_abs[0], self.maxDiff1d3dCW_rel[0],
                            'diff1d3dCW_abs,rel_lim',diff1d3dCW_abs_lim,diff1d3dCW_rel_lim[0])
                        print('ots',np.array(self.organTypes)[idCylsAll])
                        print('radii',np.array(self.radii)[idCylsAll],np.array(self.outer_radii)[idCylsAll])
                        raise Exception #only used to check that error is the same after self.update() 
                    
                for idComp in range(1, self.numComp): 
                    self.check_Ci_1d3dDiff(cellId, idComp, mol_total_[idComp -1][cellId], idCyls, diff1d3dCW_rel_lim)
                    
            for idComp in range(0, self.numComp): 
                divideTemp = sum(self.contentIn1d[idComp])
                if divideTemp == 0:
                    divideTemp = 1.
                self.sumDiff1d3dCW_rel[idComp] = self.sumDiff1d3dCW_abs[idComp]/divideTemp
        
        self.maxDiff1d3dCW_abs = comm.bcast(self.maxDiff1d3dCW_abs,root= 0) # used later as max axceptable error           
        
    def check_Ci_1d3dDiff(self, cellId, idComp, mol_total, idCyls, diff1d3dCW_rel_lim):    
        """ check difference between s1d and 3d soil for one specific solute
            @param: cellId index of the cell
            @param: idComp id of the solute component
            @param: content in the 3d soil voxel
            @param: idCyls: id of cylinder
            @param: diff1d3dCW_rel_lim max acceptable error
        """
        isDissolved = (idComp <= self.numDissolvedSoluteComp)
                    
        mol_rhizo = self.getContentCyl(idCyls, idComp, doSum = True)
        if rank == 0:
            if (mol_total    == 0) :
                diff1d3dCW_abs = abs(mol_total)
                if diff1d3dCW_abs!=0:
                    diff1d3dCW_rel = 100.
                else:
                    diff1d3dCW_rel = 0.
                    
                if ((mol_rhizo > 1e-16)):
                    print('error for component', rank, cellId)
                    print(idComp, mol_total, mol_rhizo)
                    raise Exception
            else:
                diff1d3dCW_abs = abs(mol_total  - mol_rhizo)
                diff1d3dCW_rel = abs(diff1d3dCW_abs/(mol_total ) *100)
                
            self.allDiff1d3dCW_abs[idComp][cellId] = diff1d3dCW_abs
            self.allDiff1d3dCW_rel[idComp][cellId] = diff1d3dCW_rel
            self.contentIn3d[idComp][cellId] = mol_total
            self.contentIn1d[idComp][cellId] = mol_rhizo
            
            self.sumDiff1d3dCW_abs[idComp] += diff1d3dCW_abs
            
            self.maxDiff1d3dCW_abs[idComp] = max(self.maxDiff1d3dCW_abs[idComp], diff1d3dCW_abs)
            if diff1d3dCW_abs > 1e-13:
                self.maxDiff1d3dCW_rel[idComp] = max(self.maxDiff1d3dCW_rel[idComp], diff1d3dCW_rel)
            if (diff1d3dCW_rel > diff1d3dCW_rel_lim[idComp]) and (diff1d3dCW_rel > 0.1) and (diff1d3dCW_abs > 1e-13) :
                print("check compBis", idComp,cellId, rank, 
                      'mol_rhizo',mol_rhizo,
                      'mol_total',mol_total, 
                      "error",diff1d3dCW_abs, diff1d3dCW_rel,diff1d3dCW_rel_lim[idComp],
                     'diff1d3dCW_abs_lim',diff1d3dCW_abs_lim)
                print('for water:', self.maxDiff1d3dCW_abs[0], 
                  self.maxDiff1d3dCW_rel[0], diff1d3dCW_rel_lim[0])
                raise Exception
    
    def broadcastPlantShape(self):
        """
            1) save plant shape at last time step
            2) get new relevant plant shape data on thread 0
            3) broadcst plant shape data across the threads
        """
        #backups
        self.outer_radii_old = self.outer_radii
        self.seg_length_old = self.seg_length
        self.seg2cell_old = self.seg2cell
        self.cell2seg_old = self.cell2seg
        if rank == 0:
            self.getNewCyldistribution()
          
        # bcast data
        self.nodesPos = comm.bcast(np.array(self.plantModel.get_nodes()), root = 0) 
        self.airSegs = comm.bcast(self.airSegs , root = 0) 
        self.toAddSoil = comm.bcast(self.toAddSoil, root = 0) 
        self.toAddAir = comm.bcast(self.toAddAir, root = 0) 
        self.seg_length = comm.bcast(self.seg_length, root = 0)  
        self.outer_radii = np.array( comm.bcast(self.outer_radii, root = 0))
        self.radii =np.array( comm.bcast(self.radii, root = 0)    )
        self.seg2cell = comm.bcast(self.seg2cell, root = 0)    
        self.cell2seg = comm.bcast(self.cell2seg, root = 0)         
        self.organTypes =np.array( comm.bcast(self.organTypes, root = 0) ) 
        self.rhizoVol = (self.outer_radii**2 - np.array(self.radii)**2)*np.pi*self.seg_length
        
        # get gloab indx of the new segments
        newLidxSoil = np.array([i for i in range(self.toAddSoil[rank])]) + sum(self.toAddSoil[:rank])+self.nsSoilOld 
        newLidxAir = np.array([i for i in range(self.toAddAir[rank])]) + sum(self.toAddAir[:rank])+self.nsAirOld
        self.rhizoSegsId = np.array([i for i in range(len(self.organTypes)) if i not in self.airSegs ])
        newEidxSoil = np.array([], dtype = np.int64)
        if len(newLidxSoil) > 0:
            newEidxSoil = self.rhizoSegsId[newLidxSoil]
        newEidxAir = np.array([], dtype = np.int64)
        if len(newLidxAir) > 0:
            newEidxAir = self.airSegs [newLidxAir]
        self.newEidx = np.concatenate((newEidxSoil, newEidxAir), dtype = np.int64)
        
        
        if self.debugMode :
            if len(newEidxSoil) > 0:
                assert (np.array(self.organTypes)[newEidxSoil] == 2).all()
        # id of 3d soil cells with at least one root
        self.cellWithRoots = np.unique([self.seg2cell[cylid] for cylid in self.rhizoSegsId]).reshape(-1)

        
        
    def update(self):
        """ 
            creates new 1d models or update the sape of existing 1d models
        """
        self.printData('updateBefore')
        try:
            maxlim1d3d = max(self.maxDiff1d3dCW_abs[0]*10,self.limErr1d3dAbs)
        except:
            maxlim1d3d = self.limErr1d3dAbs
        
        self.broadcastPlantShape() # compute and share plant data from thread 0 
        
        self.eidx = np.concatenate((self.eidx,np.array(self.newEidx, dtype = np.int64)), dtype = np.int64) # 1d model global id on this thread
        
        
        if self.debugMode :
            self.checkVolumeAndRadii(finishedUpdate=False) 
            self.check1d3dDiff( diff1d3dCW_abs_lim = maxlim1d3d) # check mass balance before updating size
            
        cellIds = self.getCellIds()
        
        wat_total = comm.bcast(self.soilModel.getWaterContent() , root = 0) # m3/m3 #self.soilWatVol_old
        mol_total = comm.bcast([self.soilModel.getContent(ncomp,  (ncomp <= self.numDissolvedSoluteComp)) for ncomp in range(1, self.numComp)], root = 0) 

        
        for i, cyl in enumerate(self.cyls):
            if not isinstance(cyl, AirSegment) :
                #update cylinder which have shrunk
                self.updateOld(i, cyl,smaller = True)
                
        try:
            
            # compute the leftover volume after cylinder shrinkage
            volLeftOver = np.array([self.get_Vol_leftoverI(i) for i in range(len(self.sizeSoilCell))])# cm3 scv 
            # compute the leftover water after cylinder shrinkage  
            watVol_leftover = np.array([self.get_watVol_leftoverI(i, wat_total) for i in range(len(self.sizeSoilCell))])#cm3 wat

            # cm3 => cm and cm3/cm3 
            WatPsi_leftover, theta_leftover = self.get_vol2theta_leftover(watVol_leftover[cellIds], 
                                                               volLeftOver[cellIds],cellIds )
            
            c_content_leftover = dict([(i,np.array([self.getC_content_leftoverI(i, ncomp, mol_total[ncomp -1 ]) for ncomp in range(1, self.numComp )])) for i in cellIds])# mol    

            #molar content of the phase ()water or soil the element is in # m3 wat * (mol wat/m3 wat) or m3 scv * (mol scv/m3 scv)
            phaseMol = dict([(i, np.array([theta_leftover[i]*volLeftOver[i]/1e6* self.molarDensityWat_m3 if ncomp <= self.numDissolvedSoluteComp else volLeftOver[i]/1e6*self.bulkDensity_m3 for ncomp in range(1, self.numComp)])) for i in cellIds])

            # in mol/mol wat or mol/mol bulk soil
            molFr_leftover, conc_leftover = self.getMolFrAndConcLeftover(c_content_leftover, # mol C
                                                              phaseMol,# mol wat or mol mineral soil
                                                              volLeftOver,# mol/molscv or mol/cm3 scv
                                                              cellIds)
                                                              
            
            if self.debugMode :
                try:
                    assert molFr_leftover[cellIds[0]].shape == (self.numSoluteComp,)
                except:
                    print('molFr_leftover[cellIds[0]].shape != (self.numSoluteComp,)',c_content_leftover,  molFr_leftover)
                    raise Exception
                write_file_array('volLeftOver',volLeftOver, 
                                 directory_ = self.results_dir )
                write_file_array('watVol_leftover',watVol_leftover, 
                                 directory_ = self.results_dir )
                write_file_float('theta_leftover',theta_leftover, 
                                 directory_ = self.results_dir )
                write_file_float('c_content_leftover',c_content_leftover, 
                                 directory_ = self.results_dir )
                write_file_array('molFr_leftover',molFr_leftover, 
                                 directory_ = self.results_dir )
                write_file_array('conc_leftover',conc_leftover, 
                                 directory_ = self.results_dir )
                                 

            for i, cyl in enumerate(self.cyls):
                if not isinstance(cyl, AirSegment) :
                    try:
                        self.updateOld(i, cyl,smaller = False, 
                                       thetaLeftOver = theta_leftover[self.seg2cell[cyl.gId]],
                                       konzLeftOver = conc_leftover[self.seg2cell[cyl.gId]])
                    except:
                        print('error when creating smaller cylB',i, cyl.gId,
                              theta_leftover[self.seg2cell[cyl.gId]],conc_leftover[self.seg2cell[cyl.gId]] )
                        raise Exception
                    
            for gId in self.newEidx:#only initialize the new eidx
                self.initialize_(gId,WatPsi_leftover,molFr_leftover)
        except:
            self.printData('updateAfter')
            raise Exception
            
        self.cylSoilidx = np.array([gId for lId, gId in enumerate(self.eidx) if (not isinstance(self.cyls[lId], AirSegment))])
        self.cylsSoil = np.array([cyl for cyl in self.cyls if (not isinstance(cyl, AirSegment))])

        # DO NOT SORT THEM
        self.eidx_all = self.allgatherv(np.array(self.eidx), data2share_type_default = np.int64)# all epreviously existsing segs
        self.cylSoilidx_all = self.allgatherv(np.array(self.cylSoilidx), data2share_type_default = np.int64)# all epreviously existsing segs
        
        
        if self.debugMode :
            self.printData('updateAfter')
            
            self.checkVolumeAndRadii(finishedUpdate=True)
            
            
            self.check1d3dDiff(diff1d3dCW_abs_lim = max(self.maxDiff1d3dCW_abs[0]*2,
                                                                  self.limErr1d3dAbs),
                                         diff1d3dCW_rel_lim = self.maxDiff1d3dCW_rel*2 )
            
        
        
        self.nsSoilOld = sum(self.repartitionSoilOld )
        self.nsSoilOld = comm.bcast(self.nsSoilOld, root = 0)  
        self.nsAirOld = sum(self.repartitionAirOld )
        self.nsAirOld = comm.bcast(self.nsAirOld, root = 0)  
        
    
    def printData(self, title):
        '''to see how the water and volume gets divided between the volumes'''
        print("save data update", title)
        cellIds = self.getCellIds()
        for cellId in cellIds:
            idCylsAll, idCyls = self.getIdCyllMPI(cellId, doSum = False)
            if len(idCylsAll)> 0:
                wat_rhizo_ = self.getWaterVolumesCyl(idCyls, doSum = False, reOrder=True) #cm3
                vol_rhizo_ = self.getVolumesCyl(idCyls, doSum = False, reOrder=True) #cm3
                write_file_float(title+str(cellId)+'waterSum',sum(wat_rhizo_), 
                                 directory_ = self.results_dir +"printData/")
                write_file_array(title+str(cellId)+'water',wat_rhizo_, 
                                 directory_ = self.results_dir +"printData/", fileType = '.csv' )
                write_file_array(title+str(cellId)+'vol',vol_rhizo_, 
                                 directory_ = self.results_dir +"printData/", fileType = '.csv' )
                write_file_array(title+str(cellId)+'theta',wat_rhizo_/vol_rhizo_, 
                                 directory_ = self.results_dir +"printData/", fileType = '.csv' )
                write_file_array(title+str(cellId)+'id',idCylsAll, 
                                 directory_ = self.results_dir +"printData/", fileType = '.csv' )
                for ncomp in range(1, self.numComp):
                    cont = self.getContentCyl(idCyls, idComp=ncomp, doSum = False, reOrder=True)
                    write_file_array(title+str(cellId)+'cont'+str(ncomp),cont, 
                                 directory_ = self.results_dir +"printData/", fileType = '.csv' )

    
    
    def getMolFrAndConcLeftover(self,c_content_leftover, #mol
                        phaseMol,     #mol water or mol scv
                        volLeftOver,#cm3 scv
                        cellIds):   
        """ go from left over content to molar fraction (dumux unit) and concentration """
        assert phaseMol[cellIds[0]].shape == (self.numSoluteComp,)#
        assert c_content_leftover[cellIds[0]].shape == (self.numSoluteComp,)
        CC_leftover = np.full((len(cellIds), self.numSoluteComp),np.nan)
        konz_leftover = np.zeros((len(cellIds), self.numSoluteComp))
        for i, cid in enumerate(cellIds):
            pm = phaseMol[cid]
            c_content = c_content_leftover[cid]
            CC_leftover[i][np.where(pm != 0)] = c_content[np.where(pm != 0)]/pm[np.where(pm != 0)]
            
            if volLeftOver[cid] != 0:
                konz_leftover[i,:] = c_content/volLeftOver[cid]
            
        
        assert CC_leftover.shape == (len(cellIds), self.numSoluteComp)
        
        CC_leftover = dict([(cellIds[i],CC_leftover[i]) for i in range(len(cellIds)) ]) # molar fraction
        konz_leftover = dict([(cellIds[i],konz_leftover[i]) for i in range(len(cellIds)) ]) # 

        return CC_leftover, konz_leftover #mol/mol, mol/cm3 scv
        
    def get_vol2theta_leftover(self,watVol_leftover, #cm3, len(cellIds)
                        volLeftOver,     #cm3, len(cellIds)
                        cellIds):   
        verbose =  False#(rank == 0) #(idCell== 916) and 
        theta_leftOver = np.array([watVol_leftover[idvol]/vol if vol > 0 else np.nan for idvol, vol in enumerate(volLeftOver) ])
        
        lowTheta = np.where(theta_leftOver <= self.vg_soil.theta_R )
        if len(volLeftOver[lowTheta]) > 0:
            try:
                if verbose:
                    print('get_vol2theta_leftover, some theta too low', cellIds[lowTheta],theta_leftOver[lowTheta],
                          min((theta_leftOver - self.vg_soil.theta_R)/self.vg_soil.theta_R)*100)
                assert (
                    np.logical_or(
                    (((theta_leftOver[lowTheta] - self.vg_soil.theta_R)/self.vg_soil.theta_R)*100 > -1.),
                        ((-theta_leftOver[lowTheta] + self.vg_soil.theta_R)*volLeftOver[lowTheta]<=self.maxDiff1d3dCW_abs[0] *10)
                )
                ).all()
                theta_leftOver[lowTheta] = self.vg_soil.theta_R + 1e-14
            except:
                print('min((theta_leftOver - self.vg_soil.theta_R)/self.vg_soil.theta_R)*100 < -1.', volLeftOver,theta_leftOver )
                print(volLeftOver[lowTheta],
                      theta_leftOver[lowTheta],self.vg_soil.theta_R, 
                      'lowTheta',lowTheta, 'cellIds',cellIds,self.maxDiff1d3dCW_abs[0])
                raise Exception
        highTheta = np.where(theta_leftOver > self.vg_soil.theta_S )
        if len(volLeftOver[highTheta]) > 0:
        
            try:
            
                if verbose:
                    print('get_vol2theta_leftover, some theta too high', volLeftOver[highTheta],
                      watVol_leftover[highTheta], theta_leftOver[highTheta],highTheta,
                      self.maxDiff1d3dCW_abs[0],self.vg_soil.theta_S)
                assert (np.logical_or(
                    ((theta_leftOver[highTheta] - self.vg_soil.theta_S)/self.vg_soil.theta_S)*100 < 1.,
                    (theta_leftOver[highTheta] - self.vg_soil.theta_S)*volLeftOver[highTheta]<=self.maxDiff1d3dCW_abs[0] * 10
                     )).all()
                theta_leftOver[highTheta] = self.vg_soil.theta_S 
            except:
                print('theta too high',volLeftOver[highTheta],
                      watVol_leftover[highTheta], theta_leftOver[highTheta],highTheta,
                      self.maxDiff1d3dCW_abs[0],self.vg_soil.theta_S)
                raise Exception
                
        WatPsi_leftover = np.full(theta_leftOver.shape, np.nan)
        if len(theta_leftOver) > 0:
            try:
                idtocompute = np.where(~np.isnan(theta_leftOver))
                WatPsi_leftover[idtocompute] = np.array([vg.pressure_head( tlo, self.vg_soil) for tlo in theta_leftOver[idtocompute]])
            except:
                print(theta_leftOver[lowTheta],volLeftOver[lowTheta])
                print(theta_leftOver[highTheta],volLeftOver[highTheta])
                print(theta_leftOver )
                raise Exception
        
        emptyCell = np.where(np.isnan(theta_leftOver) )
        theta_leftOver[emptyCell] = 0.
        
        return dict(zip(cellIds, WatPsi_leftover)), dict(zip(cellIds, theta_leftOver))
        
    
    

    def allgatherv(self,data2share,keepShape = False, data2share_type_default= float):
        return self.soilModel.allgatherv(data2share,keepShape = keepShape, data2share_type_default = data2share_type_default)
        
    
    def getIdCyllMPI(self,cellId, getidCylsAll = True, doSum =False):    
        idSegs = self.cell2seg[cellId]#all segments in the cell
        idCyls =np.array( list(set(idSegs).intersection(self.cylSoilidx))) #all the segments which already have cylinder
        idCyls.sort()
        if getidCylsAll:
            idCylsAll = np.array( list(set(self.cylSoilidx_all).intersection(idSegs))) # to get the newest segments at the end 
            idCylsAll.sort()
            if doSum:
                idCylsAll = len(idCylsAll)
            return idCylsAll, idCyls
        else:
            return idCyls
            
    def getXcyl(self,data2share,idCyll_=None, doSum = True, reOrder = True):
        if not doSum:
            data2share = self.allgatherv(data2share)
            if reOrder:
                if idCyll_ is None:
                    idCyll =self.eidx_all
                else:
                    idCyll = self.allgatherv(idCyll_, data2share_type_default = np.int64) #  
                try:
                    assert data2share.shape == idCyll.shape
                except:
                    print(data2share.shape, idCyll.shape)
                    raise Exception
                if len(idCyll) > 0:
                    try:
                        sorted_order = np.argsort(idCyll)
                        data2share2 = data2share[sorted_order]
                    except:
                        print('issue data2share2[idCyll] = data2share', idCyll, data2share, self.eidx)
                        raise Exception
                else:
                    assert len(data2share) == 0
            else:
                data2share2 = data2share
        else:# faster than allgather
            data2share2 = comm.allreduce( sum(data2share), op=MPI.SUM)#self.allgatherv(sum(data2share))#= sum(data2share)
        return data2share2
        
    def getLocalIdCyls(self,idCyls=None):# idCyls needs to be sorted
        """ get local from global cylinder id  for the thread """
        if idCyls is None:
            idCyls = self.eidx
        try:
            outout = np.array([ np.where(self.eidx == i)[0] for i in idCyls if  np.where(self.eidx == i)[0] < len(self.cyls)])
        except:
            print('error getLocalIdCyls', self.eidx, idCyls)
            raise Exception
        return outout.flatten()#np.where(idCyls in self.eidx)[0] 
        
        
    def getContentCyl(self,idCyls=None, idComp=1, doSum = True, reOrder = True):#mol
        isDissolved = (idComp <= self.numDissolvedSoluteComp)
        localIdCyls =   self.getLocalIdCyls(idCyls)             
        mol_rhizo = np.array([sum(self.cyls[i].getContentCyl( idComp, isDissolved)) for i in localIdCyls ]) #cm3
        return self.getXcyl(data2share=mol_rhizo, idCyll_ = idCyls,doSum=doSum, reOrder=reOrder)
        
    def getTotCContentAll(self,idCyls=None, doSum = True, reOrder = True):#mol
        localIdCyls =   self.getLocalIdCyls(idCyls)                      
        mol_rhizo = np.array([sum(self.getTotCContent(self.cyls[i] )) for i in localIdCyls ]) #mol
        return self.getXcyl(data2share=mol_rhizo,idCyll_ = idCyls, doSum=doSum, reOrder=reOrder)
        
    def getThetaCyl(self,idCyls=None, doSum = True, reOrder = True, verbose = False):#cm3
        localIdCyls =   self.getLocalIdCyls(idCyls)           
        wat_rhizo = np.array([sum(self.cyls[i].getWaterVolumesCyl()) for i in localIdCyls ]) #cm3
        V_rhizo = np.array([sum(self.cyls[i].getCellVolumes()) for i in localIdCyls ]) #cm3
        theta= wat_rhizo/V_rhizo#att for air segs
        if verbose and (len(wat_rhizo) > 0):
            print('getWaterVolumesCyl',wat_rhizo,idCyls)#
        return self.getXcyl(data2share=theta, idCyll_ = idCyls,doSum=doSum, reOrder=reOrder)
    
    def getWaterVolumesCyl(self,idCyls=None, doSum = True, reOrder = True, verbose = False):#cm3
        localIdCyls =   self.getLocalIdCyls(idCyls)           
        wat_rhizo = np.array([sum(self.cyls[i].getWaterVolumesCyl()) for i in localIdCyls ]) #cm3
        if verbose and (len(wat_rhizo) > 0):
            print('getWaterVolumesCyl',wat_rhizo,idCyls)
        return self.getXcyl(data2share=wat_rhizo, idCyll_ = idCyls,doSum=doSum, reOrder=reOrder)
        
    def getVolumesCyl(self,idCyls=None, doSum = True, reOrder = True):#cm3
        localIdCyls =   self.getLocalIdCyls(idCyls)                     
        V_rhizo = np.array([sum(self.cyls[i].getCellVolumes()) for i in localIdCyls ]) #cm3
        return self.getXcyl(data2share=V_rhizo,idCyll_ = idCyls, doSum=doSum, reOrder=reOrder)
        
    def getKrw(self,idCyls=None):#[-] unitless
        doSum = False
        reOrder = True
        localIdCyls =   self.getLocalIdCyls(idCyls)                     
        krws = np.array([self.cyls[i].getKrw()[0] for i in localIdCyls ]) 
        return self.getXcyl(data2share=krws ,idCyll_ = idCyls, doSum=doSum, reOrder=reOrder)
        
    def getDeltaR(self,idCyls=None):#[-] unitless, TODO: replace with self.dx2?
        doSum = False
        reOrder = True
        localIdCyls =   self.getLocalIdCyls(idCyls)      
        cc0 = np.array([self.cyls[i].getCellCenters()[0] for i in localIdCyls ]) 
        p0 = np.array([ self.cyls[i].getPoints()[0]  for i in localIdCyls ], dtype=object).flatten() 
        deltaR = cc0 - p0
        return self.getXcyl(data2share=deltaR ,idCyll_ = idCyls, doSum=doSum, reOrder=reOrder)
    
    def getC_rhizo(self, idComp = 1, konz = True): # if konz:  mol/m3 wat or mol/m3 scv, else: mol
        """ return component concentration or content, only in voxel with 1D-domains"""
        
        isDissolved = (idComp <= self.numDissolvedSoluteComp)
        
        contentOrKonz = np.full(len(self.sizeSoilCell), 0.)
        cellIds = self.getCellIds()
            
        for ncell, idCell in enumerate(cellIds):
            idSegs = self.cell2seg[idCell]#all segments in the cell
            idCylsAll, idCyls= self.getIdCyllMPI(idCell, doSum = True)
            if idCylsAll > 0:
                mol_rhizo = self.getContentCyl(idCyls, idComp)
                if konz:
                    if isDissolved:
                        V_rhizo = self.getWaterVolumesCyl(idCyls)/1e6
                    else:
                        V_rhizo = self.getVolumesCyl(idCyls)/1e6
                    V_tot =  V_rhizo 
                else:
                    V_tot = 1
                    V_rhizo = np.nan
                res_CC = ( mol_rhizo)/(V_tot)#mol/m3 or mol
                
                if res_CC < 0:
                    print("res_CC < 0:",rank , idComp,  mol_rhizo,V_tot ,V_rhizo, flush=True)
                    print(idCyls, idSegs, idCell, flush=True)
                    raise Exception
                
                contentOrKonz[idCell] =  res_CC   
        return contentOrKonz
        
        
    def updateConcentration(self,totContent,changeRatio, gradient, theta_old_, volumes,isWater, verbose_ = False):
        """ give the new theta or solute mole fraction to get specific total content and gradient
            volumes: either bulk soil or water volumes
        """
        verbose = verbose_#self.mpiVerboseInner# and (size > 1)) or verbose:
            
        if verbose:
            print(rank,'isWater',isWater, 'updateConcentration error', 'totContent',totContent, 
                  'gradient',gradient, 'theta_old',theta_old_, 'volumes',volumes,'changeRatio',changeRatio ) 
        
        # compute concentration using matrices
        matrix_size = self.NC -1
        sub_diag_values = -1.
        main_diag_values = 1.
        matrix = np.diag(np.full(matrix_size-1,sub_diag_values), k=1) + np.diag(np.full(matrix_size,main_diag_values), k=0) 
        matrix[-1,] = volumes
        aB = - gradient 
        aB = np.append(aB,totContent)
        SPmatrix = sparse.csc_matrix(sparse.coo_matrix(matrix))
        val_new = LA.spsolve(SPmatrix, aB, use_umfpack = True) #either that or mol fraction
        if verbose:
            print(rank, 'val_new',val_new) 
        try:
            assert abs(sum(val_new *volumes ) - totContent) < 1e-13 #check tot content ok
            assert (abs(((val_new[1:] - val_new[:-1] )) - gradient) < 1e-13).all() #check gradient ok / (cellCenters[1:] - cellCenters[:-1])
        except:
            print(rank,'isWater',isWater, 'updateConcentration error', 'val_new',val_new,'totContent',totContent, 
                  'gradient',gradient, 'theta_old',theta_old_, 'volumes',volumes,'changeRatio',changeRatio )
            raise Exception
        if isWater:
            maxVal = self.vg_soil.theta_S # keep that or set lower higher bound?
            minVal = self.theta_wilting_point # self.vg_soil.theta_R# do not use thetar => lead to -Inf pressure head
        else: # solutes
            maxVal = np.Inf
            minVal = 0.
        try:
            assert ( (minVal - totContent/sum(volumes)) < 1e-14) 
            assert ( (totContent/sum(volumes) - maxVal) < 1e-14) 
        except:
            print('(totContent/sum(volumes) < minVal) or (totContent/sum(volumes) > maxVal)',totContent,volumes,totContent/sum(volumes),'limits',minVal,maxVal,'diffs', minVal - totContent/sum(volumes),(totContent/sum(volumes) - maxVal ))
            raise Exception
        if verbose:
            print('BEFORE possible changes. new content',val_new* volumes, 'totContent',totContent,
                  'sum new content',sum(val_new* volumes),
                  'sum(val_new* volumes)-totContent', sum(val_new* volumes)-totContent,
                 'max(val_new) - maxVal ',max(val_new - maxVal))
        n_iter = 0
        
        try:
            assert (val_new > minVal).all()
            assert (val_new <= maxVal).all()
            assert abs(sum(val_new* volumes)-totContent) < 1e-14
        except:
            
            # update concentration manually if mass balance is incorect
            
            if verbose:
                print(rank, 'updateConcentration: need to adapt val_new', val_new, 'to', minVal, maxVal,totContent, volumes  )
            n_iter = 0
            while  ((max(val_new - maxVal) >  0) or  (  max(minVal - val_new) > 0)) and (n_iter <= 5):
                n_iter += 1
                if (val_new - maxVal > 0).any():
                    if (max(val_new - maxVal) < 1e-14):
                        val_new = np.minimum(val_new, maxVal)
                    else:
                        tempContent = val_new.copy()
                        idtemp = np.where(val_new >  maxVal)
                        extraElement = sum((val_new[idtemp] - maxVal) * volumes[idtemp])
                        val_new[idtemp] = maxVal

                        n_iter_ = 0
                        while  (extraElement >  1e-20) and (n_iter_ <= 5):
                            n_iter_ += 1
                            if verbose:
                                print(rank, 'updateConcentration: take out extra element content', extraElement )
                            canAdd = (maxVal - val_new)* volumes 
                            assert sum(canAdd) >= extraElement 
                            idtemp = np.where(canAdd > 0 )
                            totCtemp = sum(val_new*volumes) 
                            val_new[idtemp] = np.minimum(maxVal*volumes[idtemp],(val_new[idtemp]*volumes[idtemp]+extraElement/len(idtemp))/volumes[idtemp])
                            extraElement -=  (sum(val_new*volumes) - totCtemp) #extraElement/len(idtemp)
                        if verbose:
                            print(rank, 'updateConcentration: finished taking out extra element content', extraElement , val_new, 
                              sum(val_new * volumes), sum(tempContent * volumes),sum(val_new * volumes) - sum(tempContent * volumes))
                        try:
                            assert abs(sum(val_new * volumes) - sum(tempContent * volumes)) < 1e-14
                        except:
                            print('abs(sum(val_new * volumes) - sum(tempContent * volumes)) > 1e-14', abs(sum(val_new * volumes) - sum(tempContent * volumes)),
                                  sum(val_new * volumes), sum(tempContent * volumes), val_new * volumes,tempContent * volumes, 'totContent',totContent)
                            raise Exception
                        
                if (minVal - val_new  > 0).any():
                    if (  max(minVal -val_new) <  1e-14 ):
                        if verbose:
                            print('max(minVal -val_new) <  1e-14', max(minVal -val_new),np.maximum(val_new, minVal))
                        val_new = np.maximum(val_new, minVal)
                    else:
                        tempContent = val_new.copy()
                        if verbose:
                            print('sum(tempContent* volumes)-totContent', sum(val_new* volumes)-totContent, sum(tempContent* volumes)-totContent)
                        idtemp = np.where(val_new <  minVal )
                        if verbose:
                            print('np.where(val_new <  minVal)',val_new,idtemp)
                        missingElement = sum((minVal - val_new[idtemp]) * volumes[idtemp])

                        val_new[idtemp] = minVal
                        if verbose:
                            print('val_new[idtemp] = minVal',val_new)

                        n_iter_ = 0
                        while  (missingElement > 1e-20) and (n_iter_ <= 5):
                                n_iter_ += 1
                                if verbose:
                                    print(rank, 'updateConcentration: add missing water', missingElement,(minVal - val_new[idtemp]) * volumes[idtemp],(minVal - val_new) * volumes)
                                canTake = (val_new - minVal )* volumes 
                                if verbose:
                                    print('canTake',canTake , (val_new - minVal )* volumes )
                                assert sum(canTake) >= missingElement 
                                idtemp = np.where(canTake > 0 )
                                if verbose:
                                    print('np.where(canTake > 0 )',idtemp,val_new[idtemp])
                                totCtemp = sum(val_new*volumes) 
                                val_new[idtemp] = np.maximum(minVal,(val_new[idtemp]*volumes[idtemp] - missingElement/len(volumes[idtemp]))/volumes[idtemp])# # might be too high in some places
                                missingElement -= ( totCtemp - sum(val_new*volumes) )

                        if verbose:
                            print(rank, 'updateConcentration: finished adding missing element content', missingElement , val_new, 
                              sum(val_new * volumes), sum(tempContent * volumes),sum(val_new * volumes) - sum(tempContent * volumes))
                            print('sum(tempContent* volumes)-totContent', sum(val_new* volumes)-totContent, sum(tempContent* volumes)-totContent)
                            print('val_new',val_new,'tempContent',tempContent,
                              'diff good',abs(sum(val_new * volumes) - sum(tempContent * volumes)),'diff bad',abs(sum(val_new) - sum(tempContent)) )
                        try:
                            assert abs(sum(val_new * volumes) - sum(tempContent * volumes)) < 1e-14
                        except:
                            print('abs(sum(val_new * volumes) - sum(tempContent * volumes)) < 1e-14', abs(sum(val_new * volumes) - sum(tempContent * volumes)),
                                  sum(val_new * volumes), sum(tempContent * volumes), val_new * volumes,tempContent * volumes)
                            raise Exception

        
        if verbose:
            print('AFTER possible changes. new content',val_new* volumes, 'totContent',totContent,
                  'sum new content',sum(val_new* volumes),
                  'change ratio error', sum(val_new* volumes)-totContent,'n_iter',n_iter)
        # final check  
        assert (val_new >= minVal).all()
        assert (val_new <= maxVal).all()
        assert abs(sum(val_new* volumes)-totContent) < 1e-14
        return val_new
    
    def interAndExtraPolation_(self,pt,pointsOld_, chip, spl):
        if (pt >= min(pointsOld_)) and (pt <= max(pointsOld_)): #interpolation
            return chip(pt)
        else:#extrapolation
            return spl(pt)
        
    def interAndExtraPolation(self,pointsNew,pointsOld, valOld):
        chip = PchipInterpolator(pointsOld, valOld)#for interpolation
        spl =  CubicSpline(pointsOld, valOld, bc_type='not-a-knot') #extrapolation  
        return np.array([self.interAndExtraPolation_(pt, pointsOld, chip, spl) for pt in pointsNew])

            
    def updateOld(self, lId,  cyl,smaller = False,
                  thetaLeftOver = 0., #cm3/cm3 scv
                  konzLeftOver = None,#mol/cm3 scv 
                  verbose = False):
        """ update distribution of cylinder if it s volume has changed
            @param: lId: local cylinder id
            @param: cyl: cylinder object
            @
        """
        if konzLeftOver is None:
            konzLeftOver = np.zeros(self.numSoluteComp)
        alwaysUpdateWat = True
        
        gId = self.eidx[lId]
        assert gId == cyl.gId
        cellId = self.seg2cell[gId]
        verbose = self.mpiVerboseInner #True#(cellId== 916)#(cellId == 332) #or (cellId == 257)
        
        segsId =  np.array(self.cell2seg[cellId])
        segsId = np.array([ids for ids in segsId if self.organTypes[ids]==2 ])# only root organs
        try:
            hasNewSegs =(np.array([ bool(rootId not in self.cylSoilidx_all) for rootId in segsId]).any())
        except:
            print(rank, 'error with ahsnewsegs', lId, gId,  np.array(self.cell2seg[cellId]), np.array(self.organTypes)[np.array(self.cell2seg[cellId])])
            print(rank, segsId,self.cylSoilidx, self.cylSoilidx_all)
            raise Exception
        oldPoints = np.array(cyl.getPoints()).flatten() # cm
        lb = self.logbase            
        a_in = self.radii[gId] # cm
        a_out = self.outer_radii[gId]  # cm

        if not ((self.seg2cell[gId] >= 0) and (self.organTypes[gId] == 2)):
            points = np.array([a_in,a_out ])
            try:
                assert len(points) == len(oldPoints)
            except:
                print('wrong number of points', points, oldPoints, gId, self.seg2cell[gId], self.organTypes[gId])
                raise Exception
        else:
            points = np.logspace(np.log(a_in) / np.log(lb), np.log(a_out) / np.log(lb), self.NC, base = lb) # cm

            volOld = cyl.getCellVolumes()
            ## new shape: 
            volNew = self.getVolumes(points, self.seg_length[gId] ) # cm^3
            deltaVol = sum(volNew) - sum(volOld) 
            
            if verbose:
                print('updateOld, gId', gId,'rank',rank, self.seg2cell[gId], 'deltaVol',sum(volNew),
                      sum(volOld),deltaVol,'leftover',thetaLeftOver, #cm3/cm3 scv
                  konzLeftOver,'smaller',smaller,
                      (deltaVol >0) and (not smaller),(deltaVol <0) and smaller,
                     (((deltaVol >0) and (not smaller)) or ((deltaVol <0) and smaller)) and (abs(deltaVol) > 1e-12) )
            
            
            if (((deltaVol >0) and (not smaller)) or ((deltaVol <0) and smaller)) and (abs(deltaVol) > 1e-20):
                try:
                    assert (thetaLeftOver > -self.maxDiff1d3dCW_abs[0]) or (deltaVol <= 0.)
                except:
                    print('error in updateOld', gId, cellId, deltaVol, smaller,'leftover', thetaLeftOver,konzLeftOver)
                    raise Exception
                #max(abs(points - oldPoints)) > 1e-16:#new vertex distribution
                if verbose:
                    print(rank, 'updateOld', cellId, lId, gId)#, 'hasNewSegs',hasNewSegs )
                ## old shape:
                centersOld = np.array(cyl.getCellCenters()).flatten();
                assert len(centersOld)==9      # cm 
                centersNew = (points[1:] + points[:-1])/2  # cm
                ## change ratio
                if hasNewSegs or alwaysUpdateWat:# currently, only new segments can gain content from old segmetns.
                    #changeRatio_ =(sum(volNew) - sum(volOld))/vol3d[cellId]
                    changeRatio = min(sum(volNew)/sum(volOld), 1.)# we migth have volNew > volOld if the gain by L increase is higher than loss via r_out decrease
                else:
                    #changeRatio_ =0.
                    changeRatio = 1.
                if False:
                    try:
                        #assert ((changeRatio <= 1.) and (changeRatio > 0.))
                        assert changeRatio > 0.
                    except:
                        print('volNew', volNew,'volOld',volOld)
                        print("points",oldPoints, points)
                        print('lengths', self.seg_length_old[gId],self.seg_length[gId])
                        print('radii', self.radii[gId], self.outer_radii[gId],self.outer_radii_old[gId])
                        raise Exception
                if verbose:
                    print('\t',gId,'volNew', volNew,'volOld',volOld, 'deltaVol',deltaVol)
                    print('\t',gId,"points",oldPoints, points,'thetaLeftOver,deltaVol',thetaLeftOver,deltaVol)
                    print('\t',gId,'lengths', self.seg_length_old[gId],self.seg_length[gId])
                    print('\t',gId,'radii', self.radii[gId], self.outer_radii[gId],self.outer_radii_old[gId])
                ##  water:
                theta_old = cyl.getWaterContent() # cm3/cm3
                gradientOld = (theta_old[1:] - theta_old[:-1])#/(centersOld[1:] - centersOld[:-1])           
                gradientNew = self.interAndExtraPolation(points[1:-1],oldPoints[1:-1], gradientOld)
                
                wOld = sum(theta_old*volOld)
                maxVal = self.vg_soil.theta_S  # keep that or set lower higher bound?
                minVal = self.theta_wilting_point
                #changeRatioW_ = changeRatio_ * contentW[cellId]/wOld + 1.
                try:
                    #changeRatioW = max(min(changeRatioW_, sum(maxVal*volNew)/wOld),sum(minVal*volNew)/wOld)
                    changeRatioW = max(min(changeRatio, sum(maxVal*volNew)/wOld),sum(minVal*volNew)/wOld)
                    
                    assert ((changeRatioW <= 1.) and (changeRatioW > 0.))
                    assert (changeRatioW > 0.)
                except:
                    print('volNew', volNew,'volOld',volOld)
                    print("points",oldPoints, points)
                    print('lengths', self.seg_length_old[gId],self.seg_length[gId])
                    print('radii', self.radii[gId], self.outer_radii[gId],self.outer_radii_old[gId])
                    print('theta', theta_old)#,'changeRatioW_',changeRatioW_,contentW.shape)
                    raise Exception
                try:
                    assert wOld*changeRatioW > 0.
                    assert max(thetaLeftOver*deltaVol,0.) >= 0.
                except:
                    print('error new vals in updateOld',wOld,changeRatioW , thetaLeftOver,deltaVol)
                    raise Exception
                try:
                    try:
                        theta_new = self.updateConcentration(totContent = wOld*changeRatioW + max(thetaLeftOver*deltaVol,0.),
                                                             changeRatio=changeRatioW, 
                                                         gradient =gradientNew, theta_old_ = theta_old, volumes = volNew,isWater = True,
                                                            verbose_ = verbose)
                    except:
                        print('1st faile of updateConcentration for water, try again, try again with verbose = True. rank', rank)
                        theta_new = self.updateConcentration(totContent = wOld*changeRatioW + max(thetaLeftOver*deltaVol,0.),
                                                             changeRatio=changeRatioW, 
                                                         gradient =gradientNew, theta_old_ = theta_old, volumes = volNew,isWater = True,
                                                            verbose_ = True)
                    #newHead = np.full(theta_new.shape, np.nan)
                    #idtocompute = np.where(~np.isnan(theta_leftOver))
                    #newHead[idtocompute] = np.array([vg.pressure_head(nt, self.vg_soil) for nt in theta_new[idtocompute]])# cm
                    newHead = np.array([vg.pressure_head(nt, self.vg_soil) for nt in theta_new])# cm
                except:
                    print('issue updateConcentration, gradientNew', wOld,changeRatioW,gradientNew,points[1:-1],
                          oldPoints[1:-1], 'gradientOld',gradientOld,'volNew',volNew )
                    raise Exception
                if verbose:            
                    print('\t',gId,"x_old",theta_old* volOld,volOld )
                    print('\t',gId,"xnew", theta_new* volNew,volNew )
                    print('\t',gId,"newHead", newHead )
                    print('\t',gId,"points",oldPoints, points,centersNew)
                    print('\t',gId,"gradient",gradientOld, gradientNew)
                    print('\t',gId,"theta",theta_new, theta_old)
                    print('\t',gId,"vg param", self.vg_soil.theta_R, self.vg_soil.theta_S,self.theta_wilting_point)   
                    print('\t',gId,"changeRatio",changeRatio, changeRatioW)    
                
                    
                ## new contents:  
                molFrOld =np.array( [np.array(cyl.getSolution(nC+1)) for nC in range(self.numSoluteComp)])   #mol/mol 
                volWatNew = theta_new *volNew
                molFrNew = []
                for nComp in range(1, self.numComp):
                    if (molFrOld[nComp -1] != 0.).any():
                        isDissolved = (nComp <= self.numDissolvedSoluteComp)
                        if isDissolved: # mol phase = [cm3 phase] * [m3/cm3] * [mol phase /m3 phase] 
                            molarPhaseOld = theta_old*volOld/1e6 * self.phaseDensity(nComp ) 
                            molarPhaseNew = volWatNew/1e6 * self.phaseDensity(nComp ) # 
                        else:
                            molarPhaseOld = volOld/1e6 * self.phaseDensity(nComp )
                            molarPhaseNew = volNew/1e6 * self.phaseDensity(nComp ) 
                        gradientOld = (molFrOld[nComp -1][1:] - molFrOld[nComp -1][:-1])#/(centersOld[1:] - centersOld[:-1])   
                        gradientNew = self.interAndExtraPolation(points[1:-1],oldPoints[1:-1], gradientOld)
                        cOld = sum(molFrOld[nComp -1] *   molarPhaseOld  )    
                        
                        try:
                            assert abs(cOld - sum(cyl.getContentCyl(nComp, isDissolved )))< 1e-13
                        except:
                            print('contentoldA',cyl.getContentCyl(nComp, isDissolved), 
                                    sum(cyl.getContentCyl(nComp, isDissolved)))
                            print('contentoldB',molFrOld[nComp -1] *   molarPhaseOld, cOld,'diff', 
                                  cOld - sum(cyl.getContentCyl(nComp, isDissolved )))
                            print('watcontA',cyl.getWaterContent(), cyl.getCellVolumes() / 1e6 )
                            print('watcontB',theta_old,volOld/1e6, self.phaseDensity(nComp ) , cyl.phaseDensity(nComp ))
                            print('molFrOld',molFrOld[nComp -1])
                            print('molarPhaseOld', molarPhaseOld,np.multiply(cyl.getCellVolumes() / 1e6 , 
                                                    cyl.getWaterContent()) *cyl.molarDensityWat_m3)
                            raise Exception
                        try:
                            #changeRatio = changeRatio_ * contentC[nComp-1][cellId] /cOld + 1
                            #assert changeRatio > 0
                            try:
                                molFrNew.append(
                                    self.updateConcentration(totContent = cOld*changeRatio+ max(konzLeftOver[nComp -1]*deltaVol,0.),
                                                                         changeRatio=changeRatio,gradient =gradientNew, 
                                        theta_old_ = molFrOld[nComp -1], volumes = molarPhaseNew,isWater = False, verbose_ = False))                             
                            except:
                                print('1st faile of updateConcentration for comp no',nComp,', try again, try again with verbose = True. rank', rank,
                                     'max(konzLeftOver[nComp -1]*deltaVol,0.)',max(konzLeftOver[nComp -1]*deltaVol,0.))
                                molFrNew.append(
                                    self.updateConcentration(totContent = cOld*changeRatio + max(konzLeftOver[nComp -1]*deltaVol,0.),
                                                                         changeRatio=changeRatio,gradient =gradientNew, 
                                        theta_old_ = molFrOld[nComp -1], volumes = molarPhaseNew,isWater = False, verbose_ = True))  
                        except:
                            print('updateConcentration failed','cOld',cOld,#'contentC[nComp-1][cellId] ',contentC[nComp-1][cellId] ,
                                  'changeRatio',changeRatio,'gradientNew',gradientNew,'molarPhaseNew',
                                 molarPhaseNew,'molFrOld[nComp -1]',molFrOld[nComp -1],'molFrOld',molFrOld[nComp -1] ,
                                  'molarPhaseOld',   molarPhaseOld,'sum(volNew)',sum(volNew),
                                  'sum(volOld)',sum(volOld),'vol3d[cellId]',vol3d[cellId] )
                            raise Exception
                        if verbose:            
                            print('\t',gId,'nComp', nComp, "cOld",cOld,molFrOld[nComp -1],molarPhaseOld )
                            print('\t',gId,"molFrNew", molFrNew[nComp -1], 
                                        "cNew", molFrNew[nComp -1]* molarPhaseNew, 
                                        sum(molFrNew[nComp -1]* molarPhaseNew) )
                            print('\t',gId,"gradient", gradientOld, gradientNew,'ratio', 
                            sum(molFrNew[nComp -1]* molarPhaseNew)/sum(molFrOld[nComp -1] *   molarPhaseOld),
                                 'changeRatio',changeRatio,#'changeRatio_',changeRatio_,
                                  'cOld',cOld) #'contentWC[nComp][cellId]',contentWC[nComp][cellId],
                            print('\t',gId,"error",nComp,  
                                  abs(sum(molFrNew[nComp -1]* molarPhaseNew)/sum(molFrOld[nComp -1] *   molarPhaseOld) -changeRatio),
                                    abs(sum(molFrNew[nComp -1]* molarPhaseNew)/sum(molFrOld[nComp -1] *   molarPhaseOld) -changeRatio) < 1e-13)
                        try:
                            # assert (abs(sum(molFrNew[nComp -1]* molarPhaseNew)/sum(molFrOld[nComp -1] *   molarPhaseOld) -changeRatio) < 1e-13) or (abs(sum(molFrNew[nComp -1]* molarPhaseNew)<1e-18))
                            assert (abs(sum(molFrNew[nComp -1]* molarPhaseNew)- cOld*changeRatio - max(konzLeftOver[nComp -1]*deltaVol,0.)) < 1e-14) or (abs(sum(molFrNew[nComp -1]* molarPhaseNew)<1e-14)) or (abs(abs(sum(molFrNew[nComp -1]* molarPhaseNew)/ (cOld*changeRatio  + max(konzLeftOver[nComp -1]*deltaVol,0.)))*100 -100) < 1e-5)
                        except:
                            print('\t',rank,gId,"error",nComp, 'totContent', cOld*changeRatio,'changeRatio',changeRatio,
                                   'added',max(konzLeftOver[nComp -1]*deltaVol,0.),
                                 'new content', molFrNew[nComp -1]* molarPhaseNew, 'old content',molFrOld[nComp -1] *   molarPhaseOld,
                                  'sum new content',sum(molFrNew[nComp -1]* molarPhaseNew),#sum(molFrNew[nComp -1]* molarPhaseOld),
                                  'sum old content',sum(molFrOld[nComp -1] *   molarPhaseOld)  + max(konzLeftOver[nComp -1]*deltaVol,0.),
                                   'change ratio error',sum(molFrNew[nComp -1]* molarPhaseNew)- cOld*changeRatio - max(konzLeftOver[nComp -1]*deltaVol,0.))
                            raise Exception
                    else:
                        molFrNew.append(molFrOld[nComp -1])
                        
                self.cyls[lId] = self.initialize_dumux_nc_( gId, 
                                                            x = newHead,# cm
                                                            cAll = molFrNew, # mol/mol water or mol/mol scv
                                                            Cells = centersNew) # cm

                        
                
                
            
    
    def get_Vol_leftoverI(self, idCell):# cm3
        """ difference between volume of cylinder and volume of 3d soil voxel  """
        verbose = (rank == 0)  and self.mpiVerboseInner #(idCell== 916) and (rank == 0)#
        idCylsAll = 0
        idCyls =[]
        newVol = 0
        vol_rhizo = 0
        if idCell in self.cell2seg:
            idCylsAll, idCyls =   self.getIdCyllMPI(idCell, True, doSum = True)  

            if idCylsAll > 0: # we have old segments
                vol_rhizo= self.getVolumesCyl(idCyls) #m3
        newVol = self.sizeSoilCell[idCell] - vol_rhizo
        if verbose:
            print('get_Vol_leftoverI',rank, idCell,newVol, self.sizeSoilCell[idCell], vol_rhizo,  len(idCyls), idCylsAll)
        try:
            assert newVol > -1e-10
        except:
            print('issue in get_Vol_leftoverI',rank, idCylsAll, 
                  len(np.array(self.outer_radii)), np.array(self.outer_radii))
            print( idCell,newVol,self.sizeSoilCell[idCell], vol_rhizo)
            raise Exception
            
        return max(newVol,0.)
        
    def get_watVol_leftoverI(self, idCell, wat_total_):# cm3
        """ difference between water volume of cylinder and volume of 3d soil voxel  """
        wat_total = wat_total_[idCell] * self.sizeSoilCell[idCell] # [cm3] water contentw
        wat_rhizo = 0
        idSegs = []
        idCyls = []
        idCylsAll = 0
        if idCell in self.cell2seg:
            idSegs = self.cell2seg[idCell]#all segments in the cell
            idSegs=np.array( [ids for ids in idSegs if (self.organTypes[ids] == 2)])#only root segments
            idCylsAll, idCyls =   self.getIdCyllMPI(idCell, True, doSum = True)   # id of segments which already have a 1d model
            if idCylsAll > 0: # we have old segments
                wat_rhizo= self.getWaterVolumesCyl(idCyls) #cm3
        watVol_res = wat_total - wat_rhizo #cm3
        if watVol_res < 0.:            
            try:
                # no new segments (watVol_res should be 0) and error >= max water error -1e-13
                assert (len(idSegs) == idCylsAll) and (watVol_res > (-1e-13-self.maxDiff1d3dCW_abs[0])) # rounding error probably
                watVol_res = 0.
            except:
                print("getWatPsi_leftoverI")
                print(idCell, wat_total,wat_rhizo, wat_total - wat_rhizo,self.maxDiff1d3dCW_abs[0])
                print( len(idSegs), len(idCyls), idCylsAll)# how many old and new cylinders?
                raise Exception
                
        return watVol_res
        
    
    def getC_content_leftoverI(self, idCell, idComp, mol_total):# mol
        """ difference between C content of cylinder and of 3d soil voxel  
            @param: idCell: id of the cell
            @param: idComp: solute component id
            @param: mol_total: C content of component idComp according to 3d soil models
        
        """
        mol_total = mol_total[idCell] # mol
        mol_rhizo = 0
        mol_rhizo_ = 0
        idCyls = []
        idCylsAll = 0
        if idCell in self.cell2seg:
            idSegs = self.cell2seg[idCell]#all segments in the cell
            idSegs=np.array( [ids for ids in idSegs if (self.organTypes[ids] == 2)])#only root segments
            idCylsAll, idCyls =  self.getIdCyllMPI(idCell, True, doSum = True)   
            
            if idCylsAll > 0:#we have old segments
                mol_rhizo = self.getContentCyl(idCyls, idComp=idComp, doSum = True)#adapted for MPI
        res_CC = mol_total - mol_rhizo
        
        if res_CC < 0.:
            if (res_CC >(-1e-13-self.maxDiff1d3dCW_abs[idComp]*10)):# and (len(idSegs) == len(idCylsAll))): # rounding error probably 
                res_CC = 0.
            else:
                print("getC_content_leftoverI", idCell)
                print("res_CC = ",res_CC," < ",(-self.maxDiff1d3dCW_abs[idComp]*10),"or",len(idSegs),"=/=",idCylsAll,
                ", idComp:",idComp,'mol_total',mol_total ,
                'mol_rhizo', mol_rhizo,'self.maxDiff1d3dCW_abs',self.maxDiff1d3dCW_abs ) 
                print('lens', len(idSegs),idCylsAll, len(idCyls))# how many old and new cylinders?
                raise Exception
                
        return res_CC

                
    def initialize_(self,gId,x,cc ):
        """ create new segments for soil element content data
            @param gId: global segment index
            @param x: mean water potential per 3d soil cell (cm)
            @param cc: solute molar fraction per 3d soil cell (mol/mol)
        
        """
        if ((self.seg2cell[gId] >= 0) and (self.organTypes[gId] == 2)):
            try:
                self.cyls.append(self.initialize_dumux_nc_(gId, x[self.seg2cell[gId]], cc[self.seg2cell[gId]]))
            except:
                print('error at initialization',gId,self.seg2cell[gId],
                      self.organTypes[gId],
                    self.radii[gId], self.outer_radii[gId],self.seg_length[gId],
                    x[self.seg2cell[gId]],
                      cc[self.seg2cell[gId]],self.seg2cell[gId] )
                raise Exception
        else:
            # shoot or root aboveground
            a_in = self.radii[gId]#get Perimeter instead? not used for now anyway
            a_out = self.outer_radii[gId]
            self.cyls.append(AirSegment(a_in, a_out,self.seg_length[gId])) #psi_air computed by the photosynthesis module.
    

    def initialize_dumux_nc_(self, gId, x,                                                # cm
                                    cAll = [0.1 / 18,                                   # mol/mol wat
                                         10/ 18,                                        # mol/mol wat
                                         0.011 * 1e6/ (2700/ 60.08e-3* (1. - 0.43)),    # mol/mol scv
                                         0.05 * 1e6/(2700/ 60.08e-3* (1. - 0.43)),      # mol/mol scv    
                                         0.011 * 1e6/(2700/ 60.08e-3* (1. - 0.43)),     # mol/mol scv
                                         0.05 * 1e6/(2700/ 60.08e-3* (1. - 0.43)),      # mol/mol scv
                                         0./(2700/ 60.08e-3* (1. - 0.43)),              # mol/mol scv
                                         0./(2700/ 60.08e-3* (1. - 0.43))],             # mol/mol scv
                                         Cells = []):                                   # cm
        verbose = self.mpiVerboseInner # (self.seg2cell[gId]==916)
        a_in = self.radii[gId]
        a_out = self.outer_radii[gId]
        lId =int( np.where(self.eidx == gId)[0])
        
        if a_in < a_out:
            with StdoutRedirector(self.results_dir + 'stdcout_cpp.txt'):
                cyl = RichardsNoMPIWrapper(RichardsNCCylFoam(), self.useMoles)  # only works for RichardsCylFoam compiled without MPI
                cyl.initialize(verbose = False)
            cyl.setVGParameters([self.soil])
            lb = self.logbase
            
            points = np.logspace(np.log(a_in) / np.log(lb), np.log(a_out) / np.log(lb), 
                                 self.NC, base = lb)
            
            with StdoutRedirector(self.results_dir + 'stdcout_cpp.txt'):
                cyl.createGrid1d(points)# cm
                
            if len(self.dx2) > lId:
                self.dx2[lId] = 0.5 * (points[1] - points[0]) #when updating
            else:
                self.dx2.append(0.5 * (points[1] - points[0]))
                
            
            cyl.setParameter("Soil.BC.dzScaling", "1")
            cyl.seg_length = self.seg_length[gId]
            cyl.setParameter( "Soil.css1Function", str(self.soilModel.css1Function))
            cyl.setParameter("Problem.verbose", "0")
            cyl.setParameter("Problem.reactionExclusive", "0")    
            cyl.setParameter("Problem.segLength", str(self.seg_length[gId]))   # cm
            cyl.setParameter( "Soil.Grid.Cells",str( self.NC-1)) # -1 to go from vertices to cell (dof)
            if verbose:
                print("Soil.IC.P", cyl.dumux_str(x), Cells)
            cyl.setParameter("Soil.IC.P", cyl.dumux_str(x))# cm
            
            #default: no flux
            cyl.setInnerBC("fluxCyl", 0.)  # [cm/day] #Y 0 pressure?
            cyl.setOuterBC("fluxCyl", 0.)
            
            cyl.setParameter("Soil.MolarMass", str(self.soilModel.solidMolarMass))
            cyl.setParameter("Soil.solidDensity", str(self.soilModel.solidDensity))
            cyl.setParameter("Flux.UpwindWeight", "1")#very important because we get high solute gradient.
            cyl.setParameter("Soil.betaC", str(self.soilModel.betaC ))
            cyl.setParameter("Soil.betaO", str(self.soilModel.betaO))
            cyl.setParameter("Soil.C_S_W_thresC", str(self.soilModel.C_S_W_thresC )) #mol/cm3
            cyl.setParameter("Soil.C_S_W_thresO", str(self.soilModel.C_S_W_thresO )) #mol/cm3
            cyl.setParameter("Soil.k_decay", str(self.soilModel.k_decay))
            cyl.setParameter("Soil.k_decay2", str(self.soilModel.k_decay2 ))
            cyl.setParameter("Soil.k_DC", str(self.soilModel.k_DC  )) # 1/d
            cyl.setParameter("Soil.k_DO", str(self.soilModel.k_DO  )) # 1/d
            cyl.setParameter("Soil.k_growthC", str(self.soilModel.k_growthC))
            cyl.setParameter("Soil.k_growthO", str(self.soilModel.k_growthO))
            cyl.setParameter("Soil.K_L", str(self.soilModel.K_L))#[mol/cm3]
            cyl.setParameter("Soil.k_phi", str(self.soilModel.k_phi ))
            cyl.setParameter("Soil.k_RC", str(self.soilModel.k_RC))
            cyl.setParameter("Soil.k_RO", str(self.soilModel.k_RO ))

            cyl.setParameter("Soil.k_SC", str(self.soilModel.k_SC )) #cm^3/mol/d
            cyl.setParameter("Soil.k_SO", str(self.soilModel.k_SO )) #cm^3/mol/d
            
            cyl.setParameter("Soil.m_maxC", str(self.soilModel.m_maxC  ))# 1/d
            cyl.setParameter("Soil.m_maxO", str(self.soilModel.m_maxO  ))# 1/d
            cyl.setParameter("Soil.micro_maxC", str(self.soilModel.micro_maxC ))# 1/d
            cyl.setParameter("Soil.micro_maxO", str(self.soilModel.micro_maxO ))# 1/d
            cyl.setParameter("Soil.v_maxL", str(self.soilModel.v_maxL))#[d-1]

            cyl.setParameter("Soil.k_sorp", str(self.soilModel.k_sorp)) # mol / cm3 or mol
            cyl.setParameter("Soil.f_sorp", str(self.soilModel.f_sorp)) #[-]

            cyl.setParameter("Soil.kads", str(self.soilModel.kads)) #[cm3/mol/d]
            cyl.setParameter("Soil.kdes", str(self.soilModel.kdes)) #[1/d]            
            cyl.setParameter("Soil.CSSmax", str(self.soilModel.CSSmax)) #[mol/cm3 scv zone 1] or mol
            cyl.setParameter("Soil.alpha", str(self.soilModel.alpha)) #[1/d]


            cyl.setParameter("Soil.C_aOLim", str(self.soilModel.C_aOLim)) #[molC/cm3 scv]
            cyl.setParameter("Soil.C_aCLim", str(self.soilModel.C_aCLim)) #[molC/cm3 scv]
            cyl.setParameter("1.Component.LiquidDiffusionCoefficient", str(self.soilModel.Ds)) #m^2/s

            cyl.setParameter("2.Component.LiquidDiffusionCoefficient", str(self.soilModel.Dl)) #m^2/s

            
            for j in range( 1, self.numComp):
                cyl.setParameter( "Soil.BC.Bot.C"+str(j)+"Type", str(3))
                cyl.setParameter( "Soil.BC.Top.C"+str(j)+"Type", str(3))
                cyl.setParameter( "Soil.BC.Bot.C"+str(j)+"Value", str(0)) 
                cyl.setParameter( "Soil.BC.Top.C"+str(j)+"Value", str(0 )) 
                
            for j in range( 1, self.numComp):       
                cyl.setParameter("Soil.IC.C"+str(j), cyl.dumux_str(cAll[j-1]) ) 

            if len(Cells) > 0:#in case we update a cylinder, to allow dumux to do interpolation
                assert(len(cAll[j-1])==len(Cells))
                CellsStr = cyl.dumux_str(Cells/100)#cm -> m
                cyl.setParameter("Soil.IC.Z",CellsStr)# m
                if len(Cells)!= len(x):
                    print("Cells, x",Cells, x, len(Cells), len(x))
                    raise Exception
                for j in range( 1, self.numComp):
                    cyl.setParameter("Soil.IC.C"+str(j)+"Z",CellsStr) # m
                    if len(Cells)!= len( cAll[j-1]):
                        print("Cells,  cAll[j-1]",Cells,  cAll[j-1], 
                                len(Cells), len(cAll[j-1]), j)
                        raise Exception
                        
            cyl.setParameter("Newton.MaxRelativeShift",str(self.soilModel.MaxRelativeShift))
            cyl.setParameter("Newton.Verbosity", "0") 
            cyl.initializeProblem(maxDt = 250/(3600*24))
            
            cyl.setCriticalPressure(self.soilModel.wilting_point)  # cm pressure head
            cyl.bulkDensity_m3 = self.soilModel.bulkDensity_m3
            cyl.solidDensity =self.soilModel.solidDensity 
            cyl.solidMolarMass =self.soilModel.solidMolarMass
            cyl.solidMolDensity =self.soilModel.solidMolDensity         
            cyl.k_sorp =self.soilModel.k_sorp               
            cyl.CSSmax = self.soilModel.CSSmax               
            cyl.f_sorp =self.soilModel.f_sorp    
            cyl.css1Function = self.soilModel.css1Function
            cyl.ddt = 1.e-5
            cyl.gId = gId    
            cyl.l = self.seg_length[gId]    
            cyl.a_in = a_in    
            cyl.a_out = a_out
            ThetaCyl = cyl.getWaterContent()
            setDefault(cyl)
            try:
                assert (ThetaCyl >= self.vg_soil.theta_R).all()
                assert (ThetaCyl <= self.vg_soil.theta_S).all()
            except:
                print('issue thetaCyl',rank,ThetaCyl, self.vg_soil.theta_R, self.vg_soil.theta_S )
                raise Exception
            pHeadcyl = cyl.getSolutionHead()
            
            if verbose:
                print('end initialize_',gId,self.seg2cell[gId],'wat vol?',sum(cyl.getWaterVolumesCyl()),self.seg_length[gId],
                  pHeadcyl , x,pHeadcyl - x,'Cells',Cells)
            try:
                assert len(pHeadcyl) == (self.NC - 1)
                x_divide = np.where(x!=0,x,1)
                assert (np.logical_or( (abs((pHeadcyl - x)/x_divide)*100 < 1e-5) , 
                                       (abs(pHeadcyl - x) < 1e-9) )).all()
            except:
                print('error: issue with cylinder creations', rank)
                print('len(pHeadcyl) == (self.NC - 1)?', len(pHeadcyl), (self.NC - 1))
                print('(abs((pHeadcyl - x)/x_divide)*100 > 1e-5).all()?',abs((pHeadcyl - x)/x_divide)*100, pHeadcyl,x,'x_divide',x_divide)
                raise Exception
            
            return cyl
        else:
            print("RhizoMappedSegments.initialize_dumux_: Warning, segment {:g} might not be in domain, radii [{:g}, {:g}] cm".format(i, a_in, a_out))
            return []



    def get_inner_heads(self, shift = 0, weather:dict={}):#        
        """ matric potential at the root surface interface [cm]"""
        rsx = np.zeros((len(self.cyls),))
        for i, cyl in enumerate(self.cyls):  # run cylindrical models
            if isinstance(cyl, AirSegment):
                try:
                    rsx[i] = self.plantModel.getPsiAir(weather["ea"]/weather["es"], weather["TairC"])  # [cm]
                except:
                    print('issue weather', rank, weather["ea"], weather["es"], weather["TairC"], weather)
                    raise Exception
            else:
                rsx[i] = cyl.getInnerHead(shift)  # [cm] (in richards.py, then richards_cyl.hh)
        inner_heads= self._map(self.allgatherv( rsx)) # gathers and maps correctly
        return inner_heads

    def get_inner_solutes(self, compId = 1):
        """ matric potential at the root surface interface [mol/cm3]"""
        rsx = np.full(len(self.cyls),0.)
        isDissolved = (compId <= self.numDissolvedSoluteComp)
        for i, cyl in enumerate(self.cyls):  # run cylindrical models
            rsx[i] = cyl.getInnerSolutes( compId, isDissolved)  # [cm]
        inner_solutes= self._map(self.allgatherv(rsx)).flatten()  # gathers and maps correctly
        return inner_solutes
    
    
    def get_inner_concentrations(self):  # TODO
        raise Exception("do not use get_inner_concentrations yet")
        """ solute concentration at the root surface interface [g / cm3]"""
        rsx = np.zeros((len(self.cyls),))
        if self.mode.startswith("dumux"):
            for i, cyl in enumerate(self.cyls):  # run cylindrical models
                 rsx[i] = cyl.getInnerHead()  # [cm] ?!
            #pass  # currently zero flux !!!!!!
        else:
            raise Exception("RhizoMappedSegments.get_inner_concentrations: Warning, mode {:s} unknown".format(self.mode))
        inner_concentrations = self._map(self.allgatherv(rsx))  # gathers and maps correctly
        return inner_concentrations


    def get_dx2(self):
        """ TODO doc me AND only for mode="dumux" yet (set in initialize)"""
        dx2 = self._map(self.allgatherv(self.dx2))
        return dx2

    def get_inner_fluxes(self): # todo
        """ fluxes [cm3/day] at the root surface, i.e. inner boundary  """
        raise Exception
        

    def get_concentrationOrContent(self, compId, concentration = True):  
        """ total content or mean concentration of solute or water  """
        isAirSegs = self.airSegs 
        if self.mode.startswith("dumux"):
            if compId == 0: # water
                V_rhizo = self.getVolumesCyl(doSum = False, reOrder = True)#[cm3]
                contentRhizo = self.getWaterVolumesCyl(doSum = False, reOrder = True) # cm3
            else:
                isDissolved = (compId <= self.numDissolvedSoluteComp)
                if isDissolved: 
                    V_rhizo = self.getWaterVolumesCyl(doSum = False, reOrder = True) # cm3     
                    V_rhizo[isAirSegs] = 1 # avoid RuntimeWarning: invalid value encountered in true_divide                    
                else:
                    V_rhizo = self.getVolumesCyl(doSum = False, reOrder = True) #[cm3]                
                contentRhizo = self.getContentCyl(idComp=compId, doSum = False, reOrder = True)# mol
            assert (contentRhizo[isAirSegs] == 0.).all() #no solutes in the space around the air segments
            if not concentration:
                return contentRhizo
            concentration = contentRhizo/V_rhizo
            concentration[isAirSegs] = np.nan
            return concentration # mol/cm3
        else:
            raise Exception("RhizoMappedSegments.get_inner_concentrations: Warning, mode {:s} unknown".format(self.mode))
        concentration = self._map(self.allgatherv(rsx))  
        return concentration



    def solve(self, dt, n_iter, *argv):#dt, seg_rx, proposed_outer_fluxes,kex,proposed_outer_sol_fluxes
        """ set bc for cyl and solves it """
        
        #inner = bot = plant
        #outer = top = soil
        proposed_inner_fluxes = argv[0]
        proposed_outer_fluxes = argv[1]
        proposed_inner_fluxes_sol = argv[2]
        proposed_outer_fluxes_sol = argv[3]
        proposed_inner_fluxes_mucil = argv[4]
        proposed_outer_fluxes_mucil = argv[5]
        
        self.diffW = 0.
        self.seg_fluxes_limited = np.full(len(self.cyls), np.nan) # store for post processing
        self.seg_fluxes_limited_Out = np.full(len(self.cyls), np.nan) # store for post processing
        self.seg_fluxes_limited_sol_Out = np.full(len(self.cyls), np.nan) # store for post processing
        self.seg_fluxes_limited_mucil_Out = np.full(len(self.cyls), np.nan) # store for post processing
        self.seg_fluxes_limited_sol_In = np.full(len(self.cyls), np.nan) # store for post processing
        self.seg_fluxes_limited_mucil_In = np.full(len(self.cyls), np.nan) # store for post processing
        
        for lId, cyl in enumerate(self.cyls):  # run cylindrical models
            verbose = False
            gId = self.eidx[lId]  # for one process gId == lId
            if verbose:
                print(rank, "cyl no ",lId+1,"/",len(self.cyls),'gId')
            if isinstance(cyl, AirSegment):  
                cyl.setInnerFluxCyl(proposed_inner_fluxes[gId])
                self.seg_fluxes_limited[lId] = proposed_inner_fluxes[gId] # [cm3/day]
                self.seg_fluxes_limited_Out[lId] = proposed_outer_fluxes[gId] # [cm3/day]
                
                self.seg_fluxes_limited_sol_Out[lId] = proposed_outer_fluxes_sol[gId] # [mol/day]
                self.seg_fluxes_limited_sol_In[lId] = proposed_inner_fluxes_sol[gId] # [mol/day]
                
                self.seg_fluxes_limited_mucil_Out[lId] = proposed_outer_fluxes_mucil[gId] # [mol/day]
                self.seg_fluxes_limited_mucil_In[lId] = proposed_inner_fluxes_mucil[gId] # [mol/day]
                
            else:
                
                l = cyl.segLength#self.seg_length[gId]
                QflowIn = proposed_inner_fluxes[gId] 
                qIn = QflowIn/ (2 * np.pi * self.radii[gId] * l) # [cm3/day] -> [cm /day]
                QflowOut = proposed_outer_fluxes[gId] 
                if self.do1d1dFlow:
                    QflowOut += self.flow1d1d_w[gId]
                cyl.setInnerFluxCyl(qIn)  
                
                qOut = cyl.distributeSource(QflowOut, 0,self.numDissolvedSoluteComp)
                if QflowOut!=0.:
                    try:
                        assert (abs(sum(qOut) - QflowOut) < 1e-13) or (abs((sum(qOut) - QflowOut)/QflowOut) < 1e-13)
                    except:
                        print('distributeSourceEnd_error',sum(qOut) , QflowOut,sum(qOut) - QflowOut)
                        raise Exception
                        

                    
                     
                if not isinstance(proposed_inner_fluxes_sol,float):  # [mol/day]
                    botVal = proposed_inner_fluxes_sol[gId]
                else:
                    botVal = proposed_inner_fluxes_sol
                if not isinstance(proposed_outer_fluxes_sol,float):  # [mol/day]
                    topVal = proposed_outer_fluxes_sol[gId]
                else:
                    topVal = proposed_outer_fluxes_sol
                if not isinstance(proposed_inner_fluxes_mucil,float):  # [mol/day]
                    botVal_mucil = proposed_inner_fluxes_mucil[gId]
                else:
                    botVal_mucil = proposed_inner_fluxes_mucil
                    
                if not isinstance(proposed_outer_fluxes_mucil,float):  # [mol/day]
                    topVal_mucil = proposed_outer_fluxes_mucil[gId]
                else:
                    topVal_mucil = proposed_outer_fluxes_mucil
                    
                if self.do1d1dFlow:
                    topVal += self.flow1d1d_sol[gId]
                    topVal_mucil += self.flow1d1d_mucil[gId]
                    
                typeBC = np.full(self.numSoluteComp,3)
                
                topVals = np.array([ topVal, topVal_mucil])
                valueTopBC = np.full(self.numSoluteComp,0.)
                valueBotBC = np.full(self.numSoluteComp,0.)

                valueBotBC[0] = botVal / (2 * np.pi * self.radii[gId] * l) # [mol/day] -> [mol/cm2 /day]
                valueBotBC[1] = botVal_mucil / (2 * np.pi * self.radii[gId] * l) # [mol/day] -> [mol/cm2 /day]
                    
                    
                cyl.setSoluteBotBC(typeBC, valueBotBC)

                valueTopBC = np.zeros((self.numDissolvedSoluteComp,cyl.numberOfCellsTot))
                for nc in range(self.numDissolvedSoluteComp):
                    if (cyl.getContent(nc+1, isDissolved = (nc < self.numDissolvedSoluteComp))  == 0.).all():
                        cellIndx4Distrib = 0
                    else:
                        cellIndx4Distrib = None
                    valueTopBC[nc] = cyl.distributeSource(topVals[nc],nc+1 , 
                                                   self.numDissolvedSoluteComp, cellIndx4Distrib)

                
                
                buWBefore_ = cyl.getWaterVolumesCyl()
                buWBefore = sum( buWBefore_ )
                WaterContentOld = cyl.getWaterContent()
                SolutionHeadOld = cyl.getSolutionHead()
                buCCBefore_ = cyl.getContentCyl(1, True)
                buTotCBefore = self.getTotCContent(cyl)
                solsBefore = [np.array(cyl.getSolution(ooo + 1)).flatten() for ooo in range(self.numSoluteComp)]
                
                QflowIn_temp = QflowIn
                QflowOut_temp = QflowOut
                QCflowIn_temp =np.array([botVal,botVal_mucil])
                QCflowOut_temp =np.array([topVal,topVal_mucil])
                QflowIn_limited = 0
                QflowOut_limited = 0
                Q_in_m = np.zeros(3)
                Q_out_m = np.zeros(3)
                
                
                maxDt_temp = 250./(24.*3600.)
                maxRelShift = self.soilModel.MaxRelativeShift

                if verbose:
                    print('before solve:init [ topVal, topVal_mucil]',
                          np.array([ topVal, topVal_mucil]),'init QflowOut',QflowOut)
                    print('before solve: cyl.getSolutionHead()',np.array(cyl.getSolutionHead()))
                for ncomp in range(self.numSoluteComp):
                    try:
                        if verbose:
                            print('before solve: cyl.getSolution(ncomp + 1)',ncomp + 1,np.array(cyl.getSolution(ncomp + 1)).flatten())
                        assert (np.array(cyl.getSolution(ncomp + 1)).flatten() >= 0).all()
                    except:
                        print('(np.array(cyl.getSolution(ncomp + 1)).flatten() < 0).any()',np.array(cyl.getSolution(ncomp + 1)).flatten(),
                             cyl.getContentCyl(ncomp + 1, (ncomp < self.numDissolvedSoluteComp)))
                        raise Exception
                try:
                                        
                    redoSolve = True
                    n_iter_solve = 0
                    
                    while redoSolve:
                        verbose = True
                        try:
                            errorCOnly = False
                            cyl.ddt =min( 1.e-5,cyl.ddt)#or just reset to 1e-5?
                            
                            didReset = False
                            cyl.solve(dt)
                            # newton parameters are re-read at each 'solve()' calls
                            Q_in_m, Q_out_m = cyl.getSavedBC(self.radii[gId], 
                                                             np.array( self.outer_radii)[gId])     
                            
                            assert Q_out_m[0] == 0.
                            Q_out_m[0] = QflowOut_temp
                            assert (Q_out_m[1:] == 0.).all()
                            Q_out_m[1:(self.numFluidComp)] = QCflowOut_temp
                                
                            QflowIn_limited = Q_in_m[0]
                            QflowOut_limited = Q_out_m[0]
                            
                            
                            buWAfter_ =  cyl.getWaterVolumesCyl()
                            buWAfter = sum(buWAfter_ )
                            buCCAfter_ = cyl.getContentCyl(1, True)
                            buTotCAfter = self.getTotCContent(cyl)
                            diffWproposed = buWAfter - ( buWBefore + (QflowIn + QflowOut) * dt)
                            diffWtemp = buWAfter - ( buWBefore + (QflowIn_temp + QflowOut_temp) * dt)
                            diffWlimited = buWAfter - ( buWBefore + (QflowIn_limited + QflowOut_limited) * dt)

                            for ncomp in range(self.numSoluteComp):
                                try:
                                    assert (np.array(cyl.getSolution(ncomp + 1)).flatten() >= 0).all()
                                except:
                                    errorCOnly = True
                                    print(rank,
                                    '(np.array(cyl.getSolution(ncomp + 1)).flatten() < 0).any()', 
                                          ncomp,np.array(cyl.getSolution(ncomp + 1)).flatten() ,
                                         (np.array(cyl.getSolution(ncomp + 1)).flatten() >= 0).all())
                                    raise Exception
                            cyl.setParameter("Newton.MaxRelativeShift",
                                             str(self.soilModel.MaxRelativeShift))
                            redoSolve = False
                            cyl.setParameter("Newton.EnableResidualCriterion", "false") 
                            cyl.setParameter("Newton.EnableAbsoluteResidualCriterion", "false")
                            cyl.setParameter("Newton.SatisfyResidualAndShiftCriterion", "false")
                            
                            cyl.setParameter("Newton.MaxSteps", "18")
                            cyl.setParameter("Newton.MaxTimeStepDivisions", "10")
                        except Exception as e:     
                            print('solve Failed:', e,'rank',rank,'gId',gId,
                                  'n_iter_solve',n_iter_solve,
                                 (QflowIn_temp < 0 or QflowOut_temp < 0), 
                                  (min(QCflowOut_temp) <0. or min(QCflowIn_temp)<0.))
                            cyl.setParameter("Newton.EnableResidualCriterion", "false") # sometimes helps, sometimes makes things worse
                            cyl.setParameter("Newton.EnableAbsoluteResidualCriterion", "false")
                            cyl.setParameter("Newton.SatisfyResidualAndShiftCriterion", "false")
                            if n_iter_solve == 0:
                                cyl.setParameter("Newton.MaxSteps", "200")
                                cyl.setParameter("Newton.MaxTimeStepDivisions", "100")
                                cyl.reset();
                                cyl.resetNonlinearSolver()
                                didReset = True
                            elif n_iter_solve == 1:
                                print(rank,
                                      'soil.solve() failed. making the computation more precise')
                                cyl.setParameter("Newton.EnableResidualCriterion", "true") # sometimes helps, sometimes makes things worse
                                cyl.setParameter("Newton.EnableAbsoluteResidualCriterion", "true")
                                cyl.setParameter("Newton.SatisfyResidualAndShiftCriterion", "true")
                                cyl.setParameter("Newton.MaxRelativeShift", str(self.soilModel.MaxRelativeShift/10.))# 
                                cyl.reset();
                                cyl.resetNonlinearSolver()
                                didReset = True
                            elif (QflowIn_temp < 0 or QflowOut_temp < 0) or (min(QCflowOut_temp) <0. or min(QCflowIn_temp)<0.): #
                                cyl.reset();didReset = True # do reset before re.-distributing the sink
                                if (n_iter_solve < 20):  
                                    divVale = 0.9
                                    if verbose:
                                        print(rank,
                                          lId,gId,'errorCOnly',errorCOnly,
                                          'qflowIn',qIn,valueBotBC,
                                          'or qflowout',qOut,valueTopBC,
                                          'too low, and/or','maxDt',maxDt_temp,'too high.',
                                          'Decrease manually the lowests and maxDt_temp by',
                                              divVale,
                                             'QCflowOut_temp, QCflowIn_temp',
                                              QCflowOut_temp , QCflowIn_temp)
                                else:
                                    divVale = 0.
                                    if verbose:
                                        print('cannot solve',rank,
                                                          lId,gId,'errorCOnly',errorCOnly,
                                                          'qflowIn',qIn,valueBotBC,
                                                          'or qflowout',qOut,valueTopBC,
                                                          'too low, and/or','maxDt',maxDt_temp,'too high.',
                                                          'Decrease manually the lowests and maxDt_temp by',divVale,
                                                             'QCflowOut_temp, QCflowIn_temp',QCflowOut_temp , QCflowIn_temp)
                                
                                
                                # water
                                if (not errorCOnly) or (min(QCflowOut_temp) >0. and min(QCflowIn_temp)>0.):
                                    if verbose:
                                        print('adapt water', errorCOnly, min(QCflowOut_temp), min(QCflowIn_temp)>0.)
                                    if (QflowIn_temp < 0 or QflowOut_temp < 0) and (QflowIn_temp <= QflowOut_temp):
                                        QflowIn_temp *=divVale
                                        qIn = QflowIn_temp/ (2 * np.pi * self.radii[gId] * l) # [cm3/day] -> [cm /day]
                                        cyl.setInnerFluxCyl(qIn) 
                                        
                                    elif (QflowIn_temp < 0 or QflowOut_temp < 0):
                                        QflowOut_temp *= divVale
                                        qOut = cyl.distributeSource(QflowOut_temp, 0, l, self.numDissolvedSoluteComp)
                                            
                                # solutes
                                if (min(QCflowOut_temp) <0. or min(QCflowIn_temp)<0.) and (min(QCflowIn_temp) <= min(QCflowOut_temp)):
                                    if verbose:
                                        print('(min(QCflowOut_temp) <0. or min(QCflowIn_temp)<0.) and (min(QCflowIn_temp) <= min(QCflowOut_temp))')
                                    QCflowIn_temp[np.where(QCflowIn_temp == min(QCflowIn_temp))] *=divVale
                                    valueBotBC[:self.numDissolvedSoluteComp] = QCflowIn_temp/ (2 * np.pi * self.radii[gId] * l) # [cm3/day] -> [cm /day]
                                    cyl.setSoluteBotBC(typeBC, valueBotBC)
                                    if verbose:
                                        print('new SoluteBotBC', valueBotBC)
                                elif (min(QCflowOut_temp) <0. or min(QCflowIn_temp)<0.):
                                    if verbose:
                                        print('(min(QCflowOut_temp) <0. or min(QCflowIn_temp)<0.)')
                                    QCflowOut_temp[np.where(QCflowOut_temp == min(QCflowOut_temp))] *= divVale # or id where getSolution < 0
                                    valueTopBC = cyl.distributeSources(QCflowOut_temp, 
                                                              np.array([nc+1 for nc in range(self.numFluidComp)]), 
                                                                           self.numDissolvedSoluteComp)
                                    if verbose:
                                        print('new SoluteTopBC', valueTopBC)
                                        
                            elif (maxRelShift < 1e-5) :#random upper limit
                                if verbose:
                                    print(rank, lId,gId,' with qflowIn',qIn,
                                          ' or qflowout',qOut,'maxDt',maxDt_temp,'dt',dt,
                                      '. Increase NewtonMaxRelativeShift from',
                                          maxRelShift,'to',maxRelShift*10.)
                                maxRelShift *= 10.
                                # newton parameters are re-read at each 'solve()' calls
                                cyl.setParameter("Newton.MaxRelativeShift", str(maxRelShift))
                                cyl.reset();
                                cyl.resetNonlinearSolver()
                                didReset = True
                            else:
                                if True:#verbose:
                                    print(rank,
                                      'gId',gId,'lId',lId,'ERROR, GIVING UP. errorCOnly',errorCOnly,
                                      'qflowIn',qIn,valueBotBC,
                                      ' qflowout',qOut,valueTopBC)
                                    print("\ttheta_old",gId, WaterContentOld)
                                    print("\thead_old",gId, SolutionHeadOld)
                                    print("\ttheta_new",gId, cyl.getWaterContent())
                                    print("\thead_new", gId,cyl.getSolutionHead())
                                cyl.setParameter("Newton.MaxRelativeShift",
                                             str(self.soilModel.MaxRelativeShift))
                                redoSolve = False
                                self.solve_gave_up = True
                                cyl.reset();
                                cyl.resetNonlinearSolver()
                                didReset = True
                            
                            
                            assert didReset
                            for ncomp in range(self.numSoluteComp):
                                try:
                                    assert (np.array(cyl.getSolution(ncomp + 1)).flatten() >= 0).all()
                                except:
                                    raise Exception
                            n_iter_solve += 1
                            
                    verbose = False   
                    buWAfter_ =  cyl.getWaterVolumesCyl()
                    buWAfter = sum(buWAfter_ )
                    buCCAfter_ = cyl.getContentCyl(1, True)
                    buTotCAfter = self.getTotCContent(cyl)
                    diffW = buWAfter - ( buWBefore + (QflowIn_limited + QflowOut_limited) * dt)
                    diffC = (sum(buTotCAfter) - ( sum(buTotCBefore) \
                                    + (sum(QCflowOut_temp) + sum(QCflowIn_temp)) * dt))
                    self.diffW = max(diffW,abs(self.diffW))

                    self.seg_fluxes_limited[lId] = QflowIn_limited
                    self.seg_fluxes_limited_sol_In[lId] = Q_in_m[1]
                    self.seg_fluxes_limited_mucil_In[lId] = Q_in_m[2]
                    
                    self.seg_fluxes_limited_Out[lId] =  QflowOut_limited
                    
                    self.seg_fluxes_limited_sol_Out[lId] = Q_out_m[1] 
                    self.seg_fluxes_limited_mucil_Out[lId] = Q_out_m[2] 
                        
                    
                except:
                    print( "gId",gId,"l",l,"a_in", self.radii[gId] ,"a_out",self.outer_radii[gId],'maxRelShift',maxRelShift )
                    print("water conservation ",lId, "before", buWBefore , "all cells",buWBefore_)
                    print("QWflowIn", QflowIn ,"QWflowOut", QflowOut,"qWflowIn", qIn ,"qWflowOut", qOut, "dt", dt)
                    print("theta_old", WaterContentOld)
                    print("head_old", SolutionHeadOld)
                    print("theta_new", cyl.getWaterContent())
                    print("head_new", cyl.getSolutionHead())
                    print("QCflowIn", botVal, botVal_mucil ,"QCflowOut", topVal, topVal_mucil)
                    print("qCIn",valueBotBC,"qCOut",valueTopBC, 'diffC', 
                          ( sum(buTotCBefore) + (botVal + topVal + botVal_mucil+ topVal_mucil) * dt))
                    myStr = "RhizoMappedSegments.solve: dumux exception with boundaries in flow {:g} cm3/day, out flow {:g} cm3/day, segment radii [{:g}-{:g}] cm"
                    myStr = myStr.format(QflowIn,QflowOut, self.radii[gId], self.outer_radii[gId])
                    raise Exception(myStr)
            if verbose:
                print(rank, "FINISHED cyl no ",lId,"/",len(self.cyls),'gId', isinstance(cyl, AirSegment))
        
            
        
    def getTotCContent(self, cyl):
        if isinstance(cyl, AirSegment):
            return np.full(self.numSoluteComp,0.)
        
        return cyl.getTotCContent()
        
    
    
    def getDeff(self,Ds,phi,theta):
        return Ds* (theta**(10/3)) / (phi**2)
    
    def do1d1dFlows(self, soilVals,#wat. pot. (cm)
                        seg_valuesW,# theta (cm3/cm3)
                        seg_valuesSol,# mol/cm3 wat
                        seg_valuesmucil,# mol/cm3 wat
                        verbose=False):
        # manually set 1d-1d flow
        cellIds = self.getCellIds()
        self.flow1d1d_w = np.zeros(seg_valuesW.shape)
        self.flow1d1d_sol = np.zeros(seg_valuesW.shape)
        self.flow1d1d_mucil = np.zeros(seg_valuesW.shape)
        
            
        Ds = self.soilModel.Ds *(24*3600) *10000 # m^2/s to cm^2/d
        Dl = self.soilModel.Dl *(24*3600) *10000 # m^2/s to cm^2/d
        #Dynamic viscosity: 1e-3 [Pa s]
        #cm_per_hPa = 1/0.9806806
        #viscosity =  1e-3 /(24*3600)/100 * cm_per_hPa # cm*d
        for cyl in self.cylsSoil:
            gId = cyl.gId
            cellid = self.seg2cell[gId]
            pHeadMean = soilVals[cellid]
            thetaMean = vg.water_content( soilVals[cellid], self.vg_soil)
            # hydraulic conductivity [cm/day]
            Kf = vg.hydraulic_conductivity(pHeadMean, self.vg_soil)# 
            
            # get all seg ids
            segIds = self.cell2seg[cellid]
            #ots = self.organTypes[segIds]
            rootIds = np.array([sid for sid in segIds if ((self.organTypes[sid] == 2) and (sid != gId))])
            if verbose:
                print('\n\n1d1dflow', gId,cellid ,rootIds,segIds,'Kf',Kf,'cm/d')#,Kf*viscosity)
            #print('viscosity',viscosity/cm_per_hPa,'hPa*d')
            if len(rootIds) > 0:# at least 1 other segment
                
                #use position of node y
                pos0 =   np.array(self.nodesPos[gId+1])#cm
                posX =   np.array(self.nodesPos[rootIds+1])#cm
                if verbose:
                    print('posX_cm',pos0,posX)
                dist =( (pos0-posX)**2).sum(axis=1)**(1/2) #cm
                if verbose:
                    print('dist_cm',dist)
                pHead0 = vg.pressure_head(seg_valuesW[gId],self.soilModel.vg_soil)  #cm
                pHeadX = np.array([
                    vg.pressure_head(thetaCyl_,
                                    self.soilModel.vg_soil) for thetaCyl_ in seg_valuesW[rootIds]]) #cm
                if verbose:
                    print('pHead0_cm',pHead0,pHeadX,pHeadX - pHead0 )
                #exchange surface == 1/4 side surface of segs
                Surf = np.minimum(cyl.l*cyl.a_out,
                                  self.seg_length[rootIds]*self.outer_radii[rootIds] )*np.pi/2 
                # if flow1d1d_w > 0, gain for the cylinder
                flow1d1d_w = Kf * ( pHeadX - pHead0 )/dist *Surf #cm3/d
                self.flow1d1d_w[gId] = sum(flow1d1d_w)
                if verbose:
                    print('Surf_cm2',Surf/(np.pi/2 ),cyl.l*cyl.a_out,
                                  self.seg_length[rootIds]*self.outer_radii[rootIds] )
                if verbose:
                    print('flow1d1d_w',flow1d1d_w)
                
                CC = np.where(flow1d1d_w > 0,seg_valuesSol[rootIds],seg_valuesSol[gId])
                if verbose:
                    print('CC',CC,seg_valuesSol[gId],seg_valuesSol[rootIds])
                
                self.flow1d1d_sol[gId] =sum(flow1d1d_w * CC)
                phi = self.soilModel.vg_soil.theta_S
                
                diffS = self.getDeff(Ds,phi,thetaMean)                                                           
                # if flow1d1d_sol > 0, gain for the cylinder                               
                self.flow1d1d_sol[gId] += sum(diffS/dist * (seg_valuesSol[rootIds] - seg_valuesSol[gId])*Surf)
                if verbose:
                    print('self.flow1d1d_sol[gId]',self.flow1d1d_sol[gId],
                      self.flow1d1d_sol[gId],'diffS',diffS)
                diffL=  self.getDeff(Dl,phi,thetaMean)
                # cm^2/d /cm * [mol/cm^3] / cm = mol/cm^3/d
                # if flow1d1d_mucil > 0, gain for the cylinder 
                self.flow1d1d_mucil[gId] = sum(
                    diffL/dist * (seg_valuesmucil[rootIds] -seg_valuesmucil[gId] )*Surf)
                                                            
                if verbose:
                    print('self.flow1d1d_mucil[gId]',self.flow1d1d_mucil[gId],
                      self.flow1d1d_w[gId])
                    
                if verbose:
                    print('set1d1d',rank,'gId',gId, self.flow1d1d_w[gId] , 
                     'pHead0',pHead0,pHeadX,pHeadX - pHead0 ,'rootIds',rootIds, 
                      'watvol0',
                      seg_valuesW[gId]*np.pi*cyl.l*((cyl.a_out-cyl.a_in)**2), 
                      seg_valuesW[rootIds]*np.pi*self.seg_length[rootIds] *((self.outer_radii[rootIds]-np.array(self.radii)[rootIds])**2),
                      'cylvol',np.pi*cyl.l*((cyl.a_out-cyl.a_in)**2),
                         np.pi*self.seg_length[rootIds] *((self.outer_radii[rootIds]-np.array(self.radii)[rootIds])**2), 
                      flush=True)
        self.flow1d1d_wG = comm.gather(self.flow1d1d_w,root=0)  
        self.flow1d1d_solG = comm.gather(self.flow1d1d_sol,root=0)  
        self.flow1d1d_mucilG = comm.gather(self.flow1d1d_mucil,root=0) 
        if verbose:
            print(repr(self.flow1d1d_w))
            print(repr(self.flow1d1d_sol))
            print(repr(self.flow1d1d_mucil))
            
        gotError = 0
        #print('self.flow1d1d_wG',np.array(self.flow1d1d_wG))
        if rank == 0:
            self.flow1d1d_w = np.array(self.flow1d1d_wG).sum(axis=0)  
            self.flow1d1d_sol = np.array(self.flow1d1d_solG).sum(axis=0)  
            self.flow1d1d_mucil = np.array(self.flow1d1d_mucilG).sum(axis=0) 
            #print('self.flow1d1d_w',np.array(self.flow1d1d_w))
            assert self.flow1d1d_w.shape == seg_valuesW.shape
            assert self.flow1d1d_sol.shape == seg_valuesW.shape
            assert self.flow1d1d_mucil.shape == seg_valuesW.shape
        
            for cc in cellIds:
                segs = self.cell2seg[cc]
                if verbose:
                    print('\n\ncell',cc,segs)
                    print('W',repr(self.flow1d1d_w[segs]),sum(self.flow1d1d_w[segs]))
                    print('sol',repr(self.flow1d1d_sol[segs]),sum(self.flow1d1d_sol[segs]))
                    print('mucil',repr(self.flow1d1d_mucil[segs]),sum(self.flow1d1d_mucil[segs]))
                divw = 1 if np.max(abs(self.flow1d1d_w[segs])) == 0 else np.max(self.flow1d1d_w[segs])
                divsol = 1 if np.max(abs(self.flow1d1d_sol[segs])) == 0 else np.max(self.flow1d1d_sol[segs])
                divmucil = 1 if np.max(abs(self.flow1d1d_mucil[segs]) )== 0 else np.max(self.flow1d1d_mucil[segs])

                vals = np.array([sum(self.flow1d1d_w[segs]),
                             sum(self.flow1d1d_sol[segs]),
                             sum(self.flow1d1d_mucil[segs])])
                divs = np.array([divw,divsol,divmucil])
                try:
                    assert max(abs(vals/divs)) < 0.1
                except:
                    print(vals/divs,'vals',vals,'divs',divs)
                    print('maxs', np.max(abs(self.flow1d1d_w[segs])) , 
                          np.max(abs(self.flow1d1d_sol[segs])) , 
                          np.max(abs(self.flow1d1d_mucil[segs]))) 
                    print('FAILED')
                    gotError = 1
                    
        gotError = comm.bcast(gotError,root=0)
        if gotError:
            raise Exception
                    
                
                            
                
            
    
    def splitSoilVals(self, soilVals, #fluxes: can be <, =, or > 0
                      seg_values, # contents: are >= 0.
                      troubleShootId = -1, verbose = False, isWater = False):
        """ split soilFlux array according to the values in seg_values """
        #verbose = False
        try:
            if troubleShootId <0. :
                #cellIds = np.fromiter(self.cell2seg.keys(), dtype=int)
                cellIds = self.getCellIds()# np.array([x for x in cellIds if x >= 0])#take out air segments
            else:
                cellIds = np.array([troubleShootId])
                verbose = True            

            if verbose:
                print('enter splitsoilval',repr(soilVals), repr(seg_values), isWater)
                print('seg2cell',repr(self.cell2seg),'cellIds',repr(cellIds))
            
            assert min(seg_values) >= 0.
            assert seg_values.shape == (len(self.organTypes), )
            # sum seg_values per voxel
            organTypes = np.array(self.organTypes)

            if verbose:
                print('organTypes',repr(organTypes))
                
            if (organTypes != 2).any():
                try:
                    assert (seg_values[np.where(organTypes != 2)] == 0.).all()
                except:
                    print(seg_values, organTypes)
                    print(seg_values[np.where(organTypes != 2)])
                    raise Exception


            seg_values_voxel = np.full(len(self.organTypes), 0.)
            splitVals = np.full(len(self.organTypes), 0.)
            if (soilVals[cellIds] != 0.).any():
                for cellid in cellIds:
                    segIds = self.cell2seg[cellid]
                    if soilVals[cellid] != 0:#
                        assert (splitVals[segIds] == 0.).all()
                        ots = organTypes[segIds]
                        rootIds = np.array([sid for sid in segIds if (organTypes[sid] == 2)])
                        if len(rootIds) > 0:
                            seg_values_root = seg_values[rootIds]
                            if isWater:
                                if verbose:
                                    print('get preassuree head',seg_values_root,self.vg_soil.theta_S,self.vg_soil.theta_R)
                                seg_values_root = np.array([ vg.pressure_head(svr, self.vg_soil) for svr in seg_values_root])
                                seg_values_root = seg_values_root  - min(seg_values_root) - np.mean(seg_values_root)+1e-14
                            weightVals = np.full(len(segIds), 0.)

                            if (ots != 2).any():
                                try:
                                    assert (seg_values[segIds][np.where(ots != 2)] == 0.).all()
                                except:
                                    print(cellid, segIds, ots, seg_values[segIds])
                                    raise Exception
                            if verbose:
                                print("soilVals[cellid]", cellid, soilVals[cellid],
                                      seg_values_voxel[segIds],segIds,
                                     seg_values[segIds])
                            if soilVals[cellid] < 0:# goes away from the 1d models
                                seg_values_voxel[rootIds] = sum(seg_values_root)
                                weightVals[np.where(ots == 2)] = seg_values_root / seg_values_voxel[rootIds]

                                if verbose:
                                    print("sum(seg_values[segIds])", seg_values_voxel[segIds],
                                          repr(weightVals), sum(weightVals[np.where(ots == 2)]),
                                         (sum(weightVals) - 1.) )
                            else:# goes toward  the 1d models
                                if verbose: 
                                    print("soilVals[cellid] > 0A", seg_values_voxel[segIds],
                                          (seg_values[segIds] == 0.).all(),
                                          min(abs(seg_values_root)) == 0.,
                                          weightVals, seg_values_root)
                                if (seg_values[segIds] == 0.).all():# all seg_values == 0
                                    seg_values_voxel[rootIds] = float(len(seg_values_root))
                                    weightVals[np.where(ots == 2)] = 1./float(len(seg_values_root))
                                    if verbose: 
                                        print("(seg_values[segIds] == 0.).all()",
                                              float(len(seg_values_root)), 
                                              seg_values_voxel[rootIds] ,
                                              weightVals[np.where(ots == 2)] )
                                else:
                                    if min(abs(seg_values_root)) == 0.:# at least 1 seg_values = 0 
                                        if verbose: 
                                            print("if min(abs(seg_values_root))A",
                                                  seg_values_root,
                                                 np.where(seg_values_root== 0.), 
                                                  seg_values_root[np.where(seg_values_root == 0.)],
                                                 np.maximum(seg_values_root,1.e-14))
                                        seg_values_root = np.maximum(seg_values_root,1.e-14)
                                        #[np.where(seg_values_root == 0.)] = 1.e-14 # to avoid /0
                                        if verbose:# none == 0
                                            print("if min(abs(seg_values_root))B", seg_values_root)
                                        assert min(abs(seg_values_root)) != 0.
                                    seg_values_voxel[rootIds] = sum(1/seg_values_root)
                                    weightVals[np.where(ots == 2)] = (1 /seg_values_root) / seg_values_voxel[rootIds]
                                if verbose: 
                                    print("soilVals[cellid] < 0B", seg_values_voxel[segIds], 
                                          weightVals, seg_values[rootIds], seg_values_root)
                                    #soilVals[cellid] < 0B [1.e+14 1.e+14 0.e+00] [1.0000000e+00 8.9776452e-15] [0.20625049 0.20703416] [1.00000000e-14 1.11387784e+00

                            splitVals[segIds] = weightVals * soilVals[cellid]
                            if verbose:
                                print("splitVals[segIds]",splitVals[segIds])
                                print("sum(weightVals)", sum(weightVals), sum(splitVals[segIds]))
                            try:
                                assert ((sum(weightVals) - 1.) < 1e-13) or ( abs(sum(splitVals[segIds]) - soilVals[cellid]) < 1e-13) or (abs((sum(splitVals[segIds]) - soilVals[cellid])/(soilVals[cellid])) < 1e-13)
                            except:
                                print('ERROR SPLIT VALSA')
                                print('(sum(weightVals) - 1.) < 1e-13',rank,weightVals, sum(weightVals))
                                print(splitVals[segIds], soilVals[cellid])
                                print(sum(splitVals[segIds]), soilVals[cellid])
                                print(sum(splitVals[segIds]) - soilVals[cellid])
                                print(((sum(weightVals) - 1.) < 1e-13) , ( abs(sum(splitVals[segIds]) - soilVals[cellid]) < 1e-13) , (abs((sum(splitVals[segIds]) - soilVals[cellid])/(soilVals[cellid])) < 1e-13))
                                print(((sum(weightVals) - 1.) ) , ( abs(sum(splitVals[segIds]) - soilVals[cellid]) ) , (abs((sum(splitVals[segIds]) - soilVals[cellid])/(soilVals[cellid])) ))
                                raise Exception

                try:
                    if True:# anyway sum(splitVals)~0, sum(soilVals[cellIds])~0
                        if abs(sum(soilVals[cellIds])) > 1e-16:
                            assert abs((sum(splitVals) - sum(soilVals[cellIds]))/sum(soilVals[cellIds]))*100. < 0.1 
                        else:
                            assert abs((sum(splitVals) - sum(soilVals[cellIds]))) < 1e-16
                except:
                    print('ERROR SPLIT VALSB')
                    print('abs((sum(splitVals) - sum(soilVals[cellIds]))/sum(soilVals[cellIds]))*100. >= 0.1')
                    print(sum(splitVals), sum(soilVals),  sum(soilVals[cellIds]))
                    print('splitVals',splitVals,'soilVals',soilVals ,'seg_values',seg_values)
                    print(sum(splitVals) -  sum(soilVals[cellIds]))
                    splitVals_ = 0.
                    soilVals_ = 0.
                    for cellid in cellIds:
                        segIds = self.cell2seg[cellid]
                        splitVals_ += sum(splitVals[segIds])
                        soilVals_ +=  soilVals[cellid]
                        print(cellid,segIds,sum(splitVals[segIds]), soilVals[cellid],'organTypes[segIds]',organTypes[segIds],
                            "current",splitVals_,  soilVals_,splitVals_-  soilVals_,
                            "to",sum(splitVals),sum(soilVals[cellIds]),sum(splitVals)-sum(soilVals[cellIds]),
                            (sum(splitVals)-sum(soilVals[cellIds]))/sum(soilVals[cellIds]),
                            sum(soilVals))
                    raise Exception

                try:
                    assert (splitVals[np.where(organTypes != 2)] == 0.).all()
                except:
                    print('ERROR SPLIT VALSC')
                    print("np.where(organTypes != 2)", np.where(organTypes != 2))
                    print("splitVals",splitVals,organTypes )
                    print(splitVals[np.where(organTypes != 2)] )
                    print("seg_values",seg_values)
                    print(seg_values[np.where(organTypes != 2)] )
                    raise Exception
            return splitVals
        except:
            
            if verbose == False:
                print('\n\nrelaunch splitsoilVals with verbose = True')
                self.splitSoilVals(soilVals, #fluxes: can be <, =, or > 0
                      seg_values, # contents: are >= 0.
                      troubleShootId = -1, verbose = True, isWater = isWater)
            raise Exception 
    
    def _map(self, x):
        """Converts @param x to a numpy array and maps it to the right indices                 """
        
        indices = self.allgatherv(self.eidx, data2share_type_default = np.int64)  # gather segment indices from all threads

        
        if rank == 0:  # only for rank 0 it is not empty
            assert len(indices) == len(x), "RhizoMappedSegments._map: indices and values have different length"
            p = np.zeros((len(x),), dtype = np.float64)
            for i in range(0, len(indices)):  #
                p[indices[i]] = x[i]            
            return p
        else:
            return np.array([])

