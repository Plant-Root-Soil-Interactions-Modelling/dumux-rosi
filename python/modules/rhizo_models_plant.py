import sys


sys.path.append("../../../build-cmake/cpp/python_binding/");

import plantbox as pb
import functional.xylem_flux as xylem_flux
import sys
from functional.xylem_flux import XylemFluxPython
from richards_no_mpi_flat import RichardsNoMPIFlatWrapper as RichardsNoMPIWrapper  # Python part of cylindrcial model (a single cylindrical model is not allowed to run in parallel)
from fv.fv_grid import *
import fv.fv_richards as rich  # local pure Python cylindrical models
import functional.van_genuchten as vg
from helpful import write_file_array, write_file_float
# import scenario_setup# import setDefault
import numpy as np
import matplotlib.pyplot as plt
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank(); max_rank = comm.Get_size(); size = comm.Get_size()
import psutil
from air_models import AirSegment
from scipy import sparse
import scipy.sparse.linalg as LA
import helpful
import numbers 
from scipy.interpolate import PchipInterpolator,  CubicSpline
import pandas as pd
from plantbox import Perirhizal
# import FPItHelper


class RhizoMappedSegments(Perirhizal):#pb.MappedPlant):
    """
        Adds 1-dimensional rhizospere models to each root segment of a MappedSegments (or later MappedPlant)    
        
        modes:        
        "dumux_w"             
        "dumux_3c"        
        "dumux_10c"        
    """

    # TODO copy mapped segments constructors (!)...

    def __init__(self,  soilModel,
                 ms = None,hm = None,
                 #seedNum=None,
                 limErr1d3dAbs = 1e-11, RichardsNCCylFoam = None):
        """ @param file_name is either a pb.MappedRootSystem, pb.MappedSegments, or a string containing a rsml filename"""
        #if not seedNum is None:
        #    super().__init__(seednum = seedNum)
        if ms:
            super().__init__(ms)
        else:
            super().__init__()
        #self.ms.constantLoc = True # make that internal
        self.set_hydraulicModel(hm) # make that internal
        
        # simulation characteristics  
        self.RichardsNCCylFoam = RichardsNCCylFoam        
        self.results_dir = soilModel.results_dir # directory for output data
        self.mode = "dumux_10c" # mode of the simulation represented
        self.numFluidComp = soilModel.numFluidComp
        self.numDissolvedSoluteComp = soilModel.numDissolvedSoluteComp
        self.numSoluteComp = soilModel.numSoluteComp
        self.numComp = soilModel.numComp
        self.static_plant = False
        self.debugMode = False
        self.fixed_point_iteration = False
        self.RIFidx = np.array([n for n in range(self.numDissolvedSoluteComp)]) # solutes where we have a realised inner flux (from the 1D models)

        self.eps = 1e-14 # allow for small errors (limit of accuracy) 
        self.diff1d3dCurrant_rel_lim = 0.01
        self.maxdiff1d3dCurrant_rel_lim = 0.01
        self.MaxRelativeShift_1DS = 1e-10 # max error from within dumux

        
        
        # constants
        #self.molarDensityWat =  soilModel.molarDensityWat # [mol/cm3] 
        
        # changes with growing plant (to update)
        # plant shape        
        self.wilting_point = soilModel.wilting_point
        self.theta_wilting_point = vg.water_content( self.wilting_point, soilModel.vg_soils[0])
        self.IdCyllMPI = {}
        self.seg_length = np.array([]) 
        self.seg_length_old = np.array([])
        
        # 3d soil 
        self.soilModel = soilModel
        self.soils = self.soilModel.soils
        self.vg_soils = self.soilModel.vg_soils
        self.cellWithRoots = [] # id of the (3d) bulk soil with at least one root
        self.soilModel_cellVolumes = comm.bcast(self.soilModel.getCellVolumes(), root = 0) #cm3
        #self.solidDensity = 2.7 #self.soilModel.solidDensity #  [g/cm^3 solid]
        #self.solidMolarMass = 60.08 #self.soilModel.solidMolarMass# [g/mol]           
        #self.solidMolDensity = self.solidDensity/self.solidMolarMass # [mol / cm3 solid] =[g/cm^3 solid] / [g/mol] 
        #self.bulkDensity = self.solidMolDensity*(1.- self.soilModel.vg_soils[0].theta_S) #porosity == theta_s # [mol / cm3 scv] = [mol / cm3 solid] * [cm3 solid /cm3 space]
        
        
        
        # all 1d models
        self.eidx_all_ = [] # all 1d segs id, also the ones being created 
        self.eidx_all = [] # all created 1d segs id 
        self.eidx = np.array([], dtype = np.int64) # id of 1d segments on that thread
        self.cyls = [] # local cylinders
        self.outer_radii = None # outer radius of the 1d models (cm)
        self.rhizoVol = [] # volume of 1d models (cm3)
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
        self.cellIdleftover = np.array([]) # ids of cels with cyl that have decreased (for the leftover functions)
        self.finishedUpdate = True
        
        # error rates and fails
        self.limErr1d3dAbs = limErr1d3dAbs
        self.maxdiff1d3dCNW_absOld = np.full(self.numComp, 1e-6)
        self.maxdiff1d3dCNW_relOld = np.full(self.numComp, 0.)
        self.maxdiff1d3dCNW_abs = np.full(self.numComp, 1e-6)
        self.maxdiff1d3dCNW_rel = np.full(self.numComp, 0.)
        self.maxdiff1d3dCurrant_rel = np.inf
        self.maxdiff1d3dCurrant = np.inf
        self.sumdiff1d3dCNW_absOld = np.full(self.numComp, 0.)
        self.sumdiff1d3dCNW_relOld = np.full(self.numComp, 0.)
        self.diff1d3dCurrant_rel = 0.
        self.diff1d3dCurrant = np.inf
        self.solve_gave_up = False
        self.errWrsi = 1000.
        self.errWrsiRealInput = 1000.
        self.errW1ds = 1000.
        self.errCN1ds = 1000.
        
        self.soil_water3dAfter_old = 0 # 
        self.totCN3dAfter_eachVoxeleachComp_old = 0 # colute content [mol]
        self.rhizoWAfter_eachCyl_old = np.array([]) # water content in 1d models at the end of the time step
        self.rhizoWAfter_eachCyl = np.array([]) #cm3 water per 1d model
        self.rhizoTotCNAfter_eachCyl_old = np.array([]) # solute content in 1d models at the end of the time step
        
        # 1d1d flow
        self.do1d1dFlow = False
        self.flow1d1d_w = np.zeros(1)
        self.flow1d1d_sol = np.zeros(1)
        self.flow1d1d_mucil = np.zeros(1)
    

    def set_hydraulicModel(self, hydraulicModel):
        """ save the plant model"""
        self.hydraulicModel = hydraulicModel
        
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
            
    
    def _flat0(self, xx):
        """flattens the gathered list in rank 0, empty list for other ranks """
        return self.soilModel._flat0(xx)
            
    
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
        aboveGround = np.array([], dtype = int)
        if not (self.cell2seg.get(-1) is None):
            aboveGround = self.cell2seg.get(-1)                
        self.airSegs = np.array(list(set(np.concatenate((aboveGround,
                                                    np.where(np.array(self.organTypes) != 2)[0])) )))
        self.airSegs.sort()#to go from local segIdx to global segIdx easely
    
    def getNewCyldistribution(self):
        """ get new ids and shape of 1d models and define distribute between the threads """
        self.outer_radii = np.array(self.segOuterRadii(type = 0,vols = self.soilModel_cellVolumes )) 
                
        self.getAirSegsId()
        if len (self.airSegs) > 0:
            self.outer_radii[self.airSegs ] = np.array(self.radii)[self.airSegs ]*1.1 #dummy outer radii                
        
        ## check value outer radiis        
        if self.debugMode :
            if isinstance(self.outer_radii_old, type(self.outer_radii)):
                assert (self.outer_radii_old >= self.outer_radii[:len(self.outer_radii_old)]).all()
            
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
            vol_total = self.soilModel_cellVolumes[cellId]  # cm3.
            
            idCylsAll, idCyls =  self.getIdCyllMPI(cellId)
            
            if(len(idCylsAll) >0):
            
                # from dumux volume
                #localIdCyls =   self.getLocalIdCyls(idCyls)                     
                #V_rhizo = np.array([sum(self.cyls[i].getCellVolumes()) for i in localIdCyls ]) #cm3
                #vol_rhizo = self.getXcyl(data2share=V_rhizo,idCyll_ = idCyls, doSum=True, reOrder=False)
                # getVolumesCyl uses radii_in
                vol_rhizo =self.getVolumesCyl(idCyls=idCyls,idCylsAll=idCylsAll,doSum = True)
                if rank == 0:
                    assert (abs((vol_total - vol_rhizo)) <= 1e-10)
                
            if not finishedUpdate:
                idSegs = np.array(self.cell2seg[cellId])#all segments in the cell
                idCylsAll = np.array([ids for ids in idSegs if self.organTypes[ids]==2 ]) # 
            if(len(idCylsAll) >0):
                # from inner and outer radii
                vol_Rhizo = sum(self.rhizoVol[idCylsAll]) #self.getVolumesCyl(idCyls,doSum = True)
                if rank == 0:
                    try:
                        assert ((abs(((vol_total - vol_Rhizo)/vol_total)*100) < 1e-12) or #  error << size of soil
                                (abs(((vol_total - vol_Rhizo)/len(idCylsAll))) < 1e-12)) # error per root small
                    except:
                        print("checkVolumeAndRadii ", idCylsAll,cellId,
                                'vol_total',vol_total , vol_Rhizo,'diff',vol_total - vol_Rhizo,
                               abs(((vol_total - vol_Rhizo)/vol_total)*100) , 
                               abs(((vol_total - vol_Rhizo)/len(idCylsAll))) )
                        print('finishedUpdate',finishedUpdate)
                        print('self.rhizoVol',sum(self.rhizoVol[idCylsAll]))
                        print('self.getVolumesCyl(idCyls,doSum = False)',self.getVolumesCyl(idCyls=idCyls,idCylsAll=idCylsAll,doSum = False))
                        raise Exception
    
                    
    def check1d3dDiff(self, 
                        diff1d3dCNW_abs_lim = np.inf, # maximal axcepted absolute arror
                              verbose_ = False,
                              diff1d3dCNW_rel_lim = None, # maximal accepted relative error
                              ):
        """ check difference between soil content according to 1d and 3d soil """
        if diff1d3dCNW_rel_lim is None :
            diff1d3dCNW_rel_lim = np.full(self.numComp, np.inf)
        if diff1d3dCNW_abs_lim is None:
            diff1d3dCNW_abs_lim = self.limErr1d3dAbs
        cellIds = self.getCellIds()
        
        s =  self.soilModel
        self.alldiff1d3dCNW_rel = np.full((self.numComp,len(cellIds)), 0.)
        self.alldiff1d3dCNW_abs = np.full((self.numComp,len(cellIds)), 0.)
        self.alldiff1d3dCNW_ = np.full((self.numComp,len(cellIds)), 0.)
        self.alldiff1d3dCNW = np.full((self.numComp,len(s.getCellVolumes())), 0.)
        
        self.sumdiff1d3dCNW_abs = np.full(self.numComp, 0.)
        self.sumdiff1d3dCNW_rel = np.full(self.numComp, 0.)
        self.maxdiff1d3dCNW_abs = np.full(self.numComp, 0.)
        self.maxdiff1d3dCNW_rel = np.full(self.numComp, 0.)
        water_content = self.soilModel.getWaterContent() 
        self.soil_water3dAfter = np.multiply(water_content, s.getCellVolumes()) 
        self.totCN3dAfter_eachVoxeleachComp = np.array([s.getContent(i+1) for i in range(s.numSoluteComp)])
        self.rhizoWAfter_eachCyl = self.getWaterVolumesCyl(doSum = False, reOrder = True) #cm3 water per 1d model

        self.soil_W1d_perVoxelAfter = np.array([self.rhizoWAfter_eachCyl[
                        self.getIdCyllMPI(cellId, False)
                    ].sum()
                    for cellId in cellIds])
        self.soil_CN1d_eachCylAfter = self.getTotCNContentAll(doSum = False, reOrder = True) # total carbon content per 1d model
        self.soil_CN1d_perVoxelAfter = np.array([
            self.soil_CN1d_eachCylAfter[self.getIdCyllMPI(cellId, False)].sum(axis=0) for cellId in cellIds
        ]).T
        
        if rank == 0:
            wat_total_ =  self.soil_water3dAfter[cellIds]  # [cm3].

            self.alldiff1d3dCNW_[0] = self.soil_water3dAfter[cellIds] - self.soil_W1d_perVoxelAfter 
            self.alldiff1d3dCNW_abs[0] = abs(self.soil_water3dAfter[cellIds] - self.soil_W1d_perVoxelAfter )

            self.alldiff1d3dCNW_rel[0] = abs((
                    self.alldiff1d3dCNW_abs[0]/wat_total_
                ) *100)
            #self.sumdiff1d3dCNW_abs[0] = self.alldiff1d3dCNW_abs[0].sum()

            # only look at the relative error if the absolute error is high enough
            # self.sumdiff1d3dCNW_rel[0] = self.sumdiff1d3dCNW_abs[0]/sum(wat_total_)
            #self.maxdiff1d3dCNW_abs[0] = self.alldiff1d3dCNW_abs[0].max()
            #compErrorAboveLim = np.where(self.alldiff1d3dCNW_abs[0] > 1e-13)
            #if len(compErrorAboveLim) > 0:
            #    self.maxdiff1d3dCNW_rel[0] = np.concatenate(([0.],self.alldiff1d3dCNW_rel[0][compErrorAboveLim])).max()


            self.alldiff1d3dCNW_[1:] = np.array([
                self.totCN3dAfter_eachVoxeleachComp[idComp][cellIds] - self.soil_CN1d_perVoxelAfter[idComp]
             for idComp in range(0, self.numComp-1)])
            self.alldiff1d3dCNW_abs[1:] = np.array([abs(
                self.totCN3dAfter_eachVoxeleachComp[idComp][cellIds] - self.soil_CN1d_perVoxelAfter[idComp]
            ) for idComp in range(0, self.numComp-1)])
            self.alldiff1d3dCNW_rel[1:] = np.array([ (
                self.alldiff1d3dCNW_abs[idComp+1]/np.where(
                self.totCN3dAfter_eachVoxeleachComp[idComp][cellIds]>0.,
                self.totCN3dAfter_eachVoxeleachComp[idComp][cellIds],1.))*100. for idComp in range(0, self.numSoluteComp)])

            self.sumdiff1d3dCNW_abs = np.array([self.alldiff1d3dCNW_abs[idComp].sum() for idComp in range(self.numComp)])

            self.maxdiff1d3dCNW_abs = np.array([self.alldiff1d3dCNW_abs[idComp].max() for idComp in range(self.numComp)])

            divideTemp_ =  np.concatenate(([sum(wat_total_)], np.array([ sum(self.soil_CN1d_perVoxelAfter[idComp]) for idComp in range(self.numComp-1)])))
            divideTemp = np.array([ divideTemp_[idComp] if divideTemp_[idComp] != 0. else 1. for idComp in range(self.numComp)])
            
            self.sumdiff1d3dCNW_rel = np.array([self.sumdiff1d3dCNW_abs[idComp]/divideTemp[idComp] for idComp in range(self.numComp)])

            self.maxdiff1d3dCNW_rel = np.array([ np.concatenate(([0.],self.alldiff1d3dCNW_rel[idComp][self.alldiff1d3dCNW_abs[idComp] > 1e-13])).max() for idComp in range(self.numComp)])



            issueWater = ((self.maxdiff1d3dCNW_abs[0].max() > diff1d3dCNW_abs_lim) and 
                (self.maxdiff1d3dCNW_rel[0] > diff1d3dCNW_rel_lim[0]))
            issueSolute = [(self.maxdiff1d3dCNW_rel[idComp] > diff1d3dCNW_rel_lim[idComp]
                        ) and (self.maxdiff1d3dCNW_rel[idComp]  > 0.1) and (
                        self.maxdiff1d3dCNW_abs[idComp]  > 1e-13) for idComp in range(1, self.numComp)]

            if issueWater or any(issueSolute) :
                print("check1d3dDiff error", issueWater, issueSolute, 
                      'self.maxdiff1d3dCNW_rel',self.maxdiff1d3dCNW_rel,
                      'self.maxdiff1d3dCNW_abs',self.maxdiff1d3dCNW_abs)
                raise Exception
            
            for nc in range(self.numComp):
                self.alldiff1d3dCNW[nc][cellIds] = self.alldiff1d3dCNW_[nc]
    
    def broadcastPlantShape(self):
        """
            1) save plant shape at last time step
            2) get new relevant plant shape data on thread 0
            3) broadcst plant shape data across the threads
        """
        #backups
        self.outer_radii_old = self.outer_radii
        self.seg_length_old = self.seg_length
        self.rhizoVol_old = self.rhizoVol
        
        self.organTypes = self.ms.organTypes
        self.radii = self.ms.radii
        self.seg2cell = self.ms.seg2cell
        self.cell2seg = self.ms.cell2seg
        #node might have shifted: new length for pre-existing segments          
        self.seg_length = np.array(self.ms.segLength())
        
        if rank == 0:
            self.getNewCyldistribution()
          
        # bcast data
        self.nodesPos = comm.bcast(np.array(self.hydraulicModel.get_nodes()), root = 0) 
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

    def getCellIdleftover(self):
        cellIds = self.getCellIds()
        self.cellIdleftover = []
        for cid in cellIds:
            segIds = self.cell2seg[cid]
            cylSoil = np.intersect1d(segIds, self.rhizoSegsId)
            cylSoilOld = np.intersect1d(segIds, self.cylSoilidx_all)

            assert len(cylSoil) >= len(cylSoilOld)

            if len(cylSoil) > len(cylSoilOld):
                self.cellIdleftover.append(cid)
            else: # no new segments but they might have changed shape
                oldVol = self.rhizoVol_old[cylSoilOld]
                newVol = self.rhizoVol[cylSoil]
                Dvol = oldVol - newVol

                if ( (max(abs(Dvol)) > 1e-16) and (max(Dvol)> sum(oldVol - newVol) ) ) :
                    self.cellIdleftover.append(cid)
        self.cellIdleftover = np.array(self.cellIdleftover)

    def update(self):
        """ 
            creates new 1d models or update the sape of existing 1d models
        """
        if self.debugMode:
            self.printData('updateBefore')
            for cyl in self.cyls:                  
                gId = cyl.gId
                if self.changedLen(cyl):     
                    print('wrong segLength', gId, str(self.seg_length[gId]),
                          self.seg_length[gId], cyl.segLength,
                          (self.seg_length[gId] - cyl.segLength),
                          (self.seg_length[gId] - cyl.segLength) / cyl.segLength)
                    raise Exception

        try:
            maxlim1d3d = max(self.maxdiff1d3dCNW_abs[0]*10,self.limErr1d3dAbs)
        except:
            maxlim1d3d = self.limErr1d3dAbs
        
        self.broadcastPlantShape() # compute and share plant data from thread 0 
        self.finishedUpdate = False
        self.eidx = np.concatenate((self.eidx,np.array(self.newEidx, dtype = np.int64)), dtype = np.int64) # 1d model global id on this thread
        self.eidx_all_ = self._flat0(comm.gather(self.eidx, root=0))# all segs
        
        if rank==0:
            assert len(self.eidx_all_) == len(list(set(list(self.eidx_all_))))
        
        if self.debugMode :
            self.checkVolumeAndRadii(finishedUpdate=self.finishedUpdate) 
            if(len(self.cyls)>0):

                FPItHelper.storeNewMassData1d(self)
                FPItHelper.storeNewMassData3d(self.soilModel,self)
                self.check1d3dDiff( diff1d3dCNW_abs_lim = maxlim1d3d) # check mass balance before updating size
            
        cellIds = self.getCellIds()
        
        wat_total = self.soilModel.getWaterContent()  # m3/m3 #self.soilWatVol_old
        mol_total = np.array([self.soilModel.getContent(ncomp) for ncomp in range(1, self.numComp)]) 

        # update the shape only of air segments and rhizoseg with almost same shape

        self.getCellIdleftover()
        for i, cyl in enumerate(self.cyls):
            if not isinstance(cyl, AirSegment) :
                #update cylinder which have changed with abs(delta Vol) < epsilon
                self.updateOld(i, cyl,shapeOnly = True)
            else:
                cyl.segLength = self.seg_length[cyl.gId]

        for i, cyl in enumerate(self.cyls):
            if not isinstance(cyl, AirSegment) :
                #update cylinder which have shrunk
                self.updateOld(i, cyl,smaller = True)


        try:
            
            
            # compute the leftover volume after cylinder shrinkage
            volLeftOver = np.array([self.get_Vol_leftoverI(i) for i in range(len(self.soilModel_cellVolumes))])# cm3 scv 
            # compute the leftover water after cylinder shrinkage  
            watVol_leftover = np.array([self.get_watVol_leftoverI(i, wat_total) for i in range(len(self.soilModel_cellVolumes))])#cm3 wat

            # cm3 => cm and cm3/cm3 
            WatPsi_leftover, theta_leftover = self.get_vol2theta_leftover(watVol_leftover[cellIds], 
                                                               volLeftOver[cellIds],cellIds )
            
            c_content_leftover = dict([(i,np.array([self.getC_content_leftoverI(i, ncomp, mol_total[ncomp -1 ]) for ncomp in range(1, self.numComp )])) for i in cellIds])# mol    

            #molar content of the phase ()water or soil the element is in # m3 wat * (mol wat/m3 wat) or m3 scv * (mol scv/m3 scv)
            if rank == 0:
                phaseMol = dict([(i, np.array([theta_leftover[i] * volLeftOver[i] * self.phaseDensity(ncomp) if ncomp <= self.numDissolvedSoluteComp else volLeftOver[i] * self.phaseDensity(ncomp) for ncomp in range(1, self.numComp)])) for i in cellIds])
            else:
                phaseMol = None

            # in mol/mol wat or mol/mol bulk soil
            molFr_leftover, conc_leftover = self.getMolFrAndConcLeftover(c_content_leftover, # mol C
                                                              phaseMol,# mol wat or mol mineral soil
                                                              volLeftOver,# mol/molscv or mol/cm3 scv
                                                              cellIds)
                                                              
            
            if self.debugMode and (rank == 0):
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
            
            theta_leftover = comm.bcast(theta_leftover, root =0)
            conc_leftover = comm.bcast(conc_leftover, root =0)
            WatPsi_leftover = comm.bcast(WatPsi_leftover, root =0)
            molFr_leftover = comm.bcast(molFr_leftover, root =0)
            
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


            if self.debugMode  :
                for cyl in self.cyls:                                                       
                    gId = cyl.gId
                    if self.changedLen(cyl):       
                        print('wrong segLength',gId,str(self.seg_length[gId]),self.seg_length[gId], cyl.segLength,
                              (self.seg_length[gId] -  cyl.segLength),
                             (self.seg_length[gId] -  cyl.segLength) / cyl.segLength)
                        raise Exception

            for gId in self.newEidx:#only initialize the new eidx
                self.initialize_(gId,WatPsi_leftover,molFr_leftover)


            if self.debugMode :
                for cyl in self.cyls:           
                    gId = cyl.gId
                    if self.changedLen(cyl):
                        print('wrong segLength',gId,str(self.seg_length[gId]),self.seg_length[gId], cyl.segLength,
                              (self.seg_length[gId] -  cyl.segLength),
                             (self.seg_length[gId] -  cyl.segLength) / cyl.segLength)
                        raise Exception


        except:
            self.printData('updateAfter')
            raise Exception
            
        self.cylSoilidx = np.array([gId for lId, gId in enumerate(self.eidx) if (not isinstance(self.cyls[lId], AirSegment))])
        self.cylsSoil = np.array([cyl for cyl in self.cyls if (not isinstance(cyl, AirSegment))])

        # DO NOT SORT THEM, allgather or only gather?
        #self.eidx_all = self.allgatherv(np.array(self.eidx), data2share_type_default = np.int64)# all epreviously existsing segs
        self.cylSoilidx_all = self.allgatherv(np.array(self.cylSoilidx), data2share_type_default = np.int64)# all epreviously existsing segs
        self.eidx_all = self.eidx_all_.copy()#self.gather(np.array(self.eidx), dtype = int)# all epreviously existsing segs
        #self.cylSoilidx_all = self.gather(np.array(self.cylSoilidx))# all epreviously existsing segs
        
        self.storeIdCyllMPI()
        self.finishedUpdate = True
        
        if self.debugMode :
            self.printData('updateAfter')
            
            self.checkVolumeAndRadii(finishedUpdate=self.finishedUpdate)

            FPItHelper.storeNewMassData1d(self)
            FPItHelper.storeNewMassData3d(self.soilModel,self)
            
            self.check1d3dDiff(diff1d3dCNW_abs_lim = maxlim1d3d)
            
        
        
        self.nsSoilOld = sum(self.repartitionSoilOld )
        self.nsSoilOld = comm.bcast(self.nsSoilOld, root = 0)  
        self.nsAirOld = sum(self.repartitionAirOld )
        self.nsAirOld = comm.bcast(self.nsAirOld, root = 0)

        
    
    def printData(self, title):
        '''to see how the water and volume gets divided between the volumes'''
        cellIds = self.getCellIds()
        for cellId in cellIds:
            idCylsAll, idCyls = self.getIdCyllMPI(cellId, doSum = False)
            if len(idCylsAll)> 0:
                wat_rhizo_ = self.getWaterVolumesCyl(idCyls, doSum = False, reOrder=True) #cm3
                vol_rhizo_ = self.getVolumesCyl(idCyls=idCyls,idCylsAll=idCylsAll, doSum = False, reOrder=True) #cm3
                if rank == 0:
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
        if rank == 0:
            assert phaseMol[cellIds[0]].shape == (self.numSoluteComp,)#
            assert c_content_leftover[cellIds[0]].shape == (self.numSoluteComp,)
            CC_leftover = np.full((len(cellIds), self.numSoluteComp),np.nan)
            konz_leftover = np.zeros((len(cellIds), self.numSoluteComp))
            for i, cid in enumerate(cellIds):
                if cid in self.cellIdleftover:
                    pm = phaseMol[cid]
                    c_content = c_content_leftover[cid]
                    CC_leftover[i][np.where(pm != 0)] = c_content[np.where(pm != 0)]/pm[np.where(pm != 0)]

                    if volLeftOver[cid] != 0:
                        konz_leftover[i,:] = c_content/volLeftOver[cid]
                
            
            assert CC_leftover.shape == (len(cellIds), self.numSoluteComp)
            
            CC_leftover = dict([(cellIds[i],CC_leftover[i]) for i in range(len(cellIds)) ]) # molar fraction
            konz_leftover = dict([(cellIds[i],konz_leftover[i]) for i in range(len(cellIds)) ]) # 

            return CC_leftover, konz_leftover #mol/mol, mol/cm3 scv
        else:
            return None, None
        
    def get_vol2theta_leftover(self,watVol_leftover, #cm3, len(cellIds)
                        volLeftOver,     #cm3, len(cellIds)
                        cellIds):   
        if rank == 0:
            verbose =  False#(rank == 0) #(idCell== 916) and 
            theta_leftOver = np.array([watVol_leftover[idvol]/vol if ((cellIds[idvol] in self.cellIdleftover) and (vol>0.))
                                       else np.nan
                                       for idvol, vol in enumerate(volLeftOver) ])
            
            lowTheta = np.where(theta_leftOver <= self.vg_soils[0].theta_R )
            if len(volLeftOver[lowTheta]) > 0:
                try:
                    if verbose:
                        print('get_vol2theta_leftover, some theta too low', cellIds[lowTheta],theta_leftOver[lowTheta],
                              min((theta_leftOver - self.vg_soils[0].theta_R)/self.vg_soils[0].theta_R)*100)
                    if  self.fixed_point_iteration:
                        assert (
                            np.logical_or(
                            (((theta_leftOver[lowTheta] - self.vg_soils[0].theta_R)/self.vg_soils[0].theta_R)*100 > -1.),
                                ((-theta_leftOver[lowTheta] + self.vg_soils[0].theta_R)*volLeftOver[lowTheta]<=self.maxdiff1d3dCNW_abs[0] *10)
                        )
                        ).all()
                    theta_leftOver[lowTheta] = self.vg_soils[0].theta_R + self.eps
                except:
                    print('min((theta_leftOver - self.vg_soils[0].theta_R)/self.vg_soils[0].theta_R)*100 < -1.', volLeftOver,theta_leftOver )
                    print(volLeftOver[lowTheta],
                          theta_leftOver[lowTheta],self.vg_soils[0].theta_R, 
                          'lowTheta',lowTheta, 'cellIds',cellIds,self.maxdiff1d3dCNW_abs[0])
                    #raise Exception
            highTheta = np.where(theta_leftOver > self.vg_soils[0].theta_S )
            if len(volLeftOver[highTheta]) > 0:
            
                try:
                
                    if verbose:
                        print('get_vol2theta_leftover, some theta too high', volLeftOver[highTheta],
                          watVol_leftover[highTheta], theta_leftOver[highTheta],highTheta,
                          self.maxdiff1d3dCNW_abs[0],self.vg_soils[0].theta_S)
                    if  self.fixed_point_iteration:
                        assert (np.logical_or(
                            ((theta_leftOver[highTheta] - self.vg_soils[0].theta_S)/self.vg_soils[0].theta_S)*100 < 1.,
                            (theta_leftOver[highTheta] - self.vg_soils[0].theta_S)*volLeftOver[highTheta]<=self.maxdiff1d3dCNW_abs[0] * 10
                             )).all()
                    theta_leftOver[highTheta] = self.vg_soils[0].theta_S 
                except:
                    print('theta too high',volLeftOver[highTheta],
                          watVol_leftover[highTheta], theta_leftOver[highTheta],highTheta,
                          self.maxdiff1d3dCNW_abs[0],self.vg_soils[0].theta_S)
                    raise Exception
                    
            WatPsi_leftover = np.full(theta_leftOver.shape, np.nan)
            if len(theta_leftOver) > 0:
                try:
                    idtocompute = np.where(~np.isnan(theta_leftOver))
                    WatPsi_leftover[idtocompute] = np.array([vg.pressure_head( tlo, self.vg_soils[0]) for tlo in theta_leftOver[idtocompute]])
                except:
                    print(theta_leftOver[lowTheta],volLeftOver[lowTheta])
                    print(theta_leftOver[highTheta],volLeftOver[highTheta])
                    print(theta_leftOver )
                    raise Exception
            
            emptyCell = np.where(np.isnan(theta_leftOver) )
            theta_leftOver[emptyCell] = 0.
            
            return dict(zip(cellIds, WatPsi_leftover)), dict(zip(cellIds, theta_leftOver))
        else:
            return None, None
    
    

    def allgatherv(self,data2share,keepShape = False, data2share_type_default= float):
        return self.soilModel.allgatherv(data2share,keepShape = keepShape, data2share_type_default = data2share_type_default)
        
    
    def getIdCyllMPI(self,cellId, getidCylsAll = True, doSum =False):
        if cellId in self.IdCyllMPI.keys():
            idCylsAll, idCyls = self.IdCyllMPI[cellId]
        else:
            idCylsAll, idCyls = np.array([],dtype = int), np.array([],dtype = int)
        if getidCylsAll:
            if doSum:
                idCylsAll = len(idCylsAll)
            return idCylsAll, idCyls
        else:
            return idCyls
        
    def storeIdCyllMPI(self): 
        self.IdCyllMPI.clear()
        cellIds = self.getCellIds()
        self.IdCyllMPI = {cellId:self.getIdCyllMPI_(cellId, 
                                               getidCylsAll = True, 
                                               doSum =False
                                              ) for cellId in cellIds}
        
    def getIdCyllMPI_(self, cellId, getidCylsAll=True, doSum=False):    
        idSegs = set(self.cell2seg[cellId])  # Convert idSegs to a set

        # Perform intersection and convert to a sorted numpy array
        idCyls = np.sort(np.array(list(idSegs.intersection(self.cylSoilidx))))  

        if getidCylsAll:
            idCylsAll = np.sort(np.array(list(idSegs.intersection(self.cylSoilidx_all))))  # Get the newest segments at the end
            if doSum:
                return len(idCylsAll), idCyls
            return idCylsAll, idCyls
        else:
            return idCyls
            
    def getXcyl(self,data2share,idCyll_=None, doSum = True, reOrder = True):
    
        data2share = self._flat0(comm.gather(data2share,root=0))
        if not doSum and reOrder:
            if idCyll_ is None:
                idCyll =self.eidx_all
            else:
                idCyll = self._flat0(comm.gather(idCyll_))#, dtype = int) #  
            #try:
            #    assert data2share.shape == idCyll.shape
            #except:
            #    print(data2share.shape, idCyll.shape)
            #    raise Exception
            
            if len(idCyll) > 0 and (rank == 0):
                #idCyll = np.array(idCyll).reshape(-1)#, dtype = int) 
                try:
                    sorted_order = np.argsort(idCyll)
                    data2share = data2share[sorted_order]
                except:
                    print('issue data2share2[idCyll] = data2share','rank',rank,  'idCyll',idCyll_,idCyll, 
                    'data2share',data2share, 'self.eidx',self.eidx)
                    raise Exception

        if doSum and (rank == 0):
            data2share = sum(data2share)        
        return data2share
        
    def getLocalIdCyls(self,idCyls=None):# idCyls needs to be sorted
        """ get local from global cylinder id  for the thread """
        if idCyls is None:
            idCyls = self.eidx
        try:
            outout = np.array([ np.where(self.eidx == i)[0] for i in idCyls if  np.where(self.eidx == i)[0] < len(self.cyls)])
        except:
            print('error getLocalIdCyls', self.eidx, idCyls)
            raise Exception
        return outout.flatten()
        
        
    def getContentCyl(self,idCyls=None, idComp=1, doSum = True, reOrder = True, verbose=False):#mol
        localIdCyls =   self.getLocalIdCyls(idCyls)      
        mol_rhizo = np.array([self.cyls[i].getContent( idComp).sum() for i in localIdCyls ]) #cm3
        
        try:
            assert mol_rhizo.shape == (len(localIdCyls),)
        except:
            print(' mol_rhizo.shape , (len(localIdCyls),)',rank, mol_rhizo.shape , (len(localIdCyls),))
            raise Exception
        if verbose:
            print(' mol_rhizo.shape , (len(localIdCyls),)',rank, mol_rhizo.shape , (len(localIdCyls),),
                 'mol_rhizo',mol_rhizo)
        return self.getXcyl(data2share=mol_rhizo, idCyll_ = idCyls,doSum=doSum, reOrder=reOrder)
        
    def getTotCNContentAll(self,idCyls=None, doSum = True, reOrder = True):#mol
        localIdCyls =   np.array(self.getLocalIdCyls(idCyls))         
        # sum over all cells for each solute
        localCyls = np.array(self.cyls)[localIdCyls]
        mol_rhizo = np.array([ np.array([cyl.getContent(i+1) for i in range(self.numSoluteComp)]).sum(axis=1) for cyl in localCyls ]) #mol
        return self.getXcyl(data2share=mol_rhizo,idCyll_ = idCyls, doSum=doSum, reOrder=reOrder)
        
    def getThetaCyl(self,idCyls=None, doSum = True, reOrder = True, verbose = False):#cm3
        localIdCyls =   self.getLocalIdCyls(idCyls)           
        wat_rhizo = np.array([sum(self.cyls[i].getWaterVolumes()) for i in localIdCyls ]) #cm3
        V_rhizo = np.array([sum(self.cyls[i].getCellVolumes()) for i in localIdCyls ]) #cm3
        theta= wat_rhizo/V_rhizo#att for air segs
        if verbose and (len(wat_rhizo) > 0):
            print('getWaterVolumes',wat_rhizo,idCyls)#
        return self.getXcyl(data2share=theta, idCyll_ = idCyls,doSum=doSum, reOrder=reOrder)
    
    def getWaterVolumesCyl(self,idCyls=None, doSum = True, reOrder = True, verbose = False):#cm3
        localIdCyls =   self.getLocalIdCyls(idCyls)  
        wat_rhizo = np.array([self.cyls[i].getWaterVolumes().sum() for i in localIdCyls ]) #cm3
        return self.getXcyl(data2share=wat_rhizo, idCyll_ = idCyls,doSum=doSum, reOrder=reOrder)
        
    def getVolumesCyl(self,idCyls=None,idCylsAll=None,  doSum = True, reOrder = True):#cm3
        if self.finishedUpdate:
            if rank == 0:
                if idCylsAll is None:
                    # no need to use eidx_all_: rhizoVol is already in correct order
                    idCylsAll = np.array([i for i in range(len(self.eidx_all_))])# self.eidx_all_
                V_rhizo = self.rhizoVol[idCylsAll]         
                if doSum and (not isinstance(V_rhizo, numbers.Number)):
                    try:
                        V_rhizo = sum(V_rhizo)
                    except:
                        print('getVolumesCyl V_rhizo',V_rhizo,idCylsAll,
                              self.rhizoVol[idCylsAll],self.rhizoVol)
                        raise Exception
                return V_rhizo 
            else:
                return np.array([])   
        else:
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
    
    def getC_rhizo(self, idComp = 1, konz = True): # if konz:  mol/cm3 wat or mol/cm3 scv, else: mol
        """ return component concentration or content, only in voxel with 1D-domains"""
        
        isDissolved = (idComp <= self.numDissolvedSoluteComp)
        
        contentOrKonz = np.full(len(self.soilModel_cellVolumes), 0.)
        cellIds = self.getCellIds()
            
        for ncell, idCell in enumerate(cellIds):
            idSegs = self.cell2seg[idCell]#all segments in the cell
            idCylsAll, idCyls= self.getIdCyllMPI(idCell, doSum = True)
            if idCylsAll > 0:
                mol_rhizo = self.getContentCyl(idCyls, idComp)
                if konz:
                    if isDissolved:
                        V_rhizo = self.getWaterVolumesCyl(idCyls)
                    else:
                        V_rhizo = self.getVolumesCyl(idCyls=idCyls,idCylsAll=idCylsAll)
                    V_tot =  V_rhizo 
                else:
                    V_tot = 1
                    V_rhizo = np.nan
                if rank == 0:
                    res_CC = ( mol_rhizo)/(V_tot)#mol/cm3 or mol
                    if res_CC < 0:
                        print("res_CC < 0:",rank , idComp,  mol_rhizo,V_tot ,V_rhizo, flush=True)
                        print(idCyls, idSegs, idCell, flush=True)
                        raise Exception
                    contentOrKonz[idCell] =  res_CC   
        return contentOrKonz
        
        
    def update_concentration(self, totContent, changeRatio, gradient, volumes, isWater, verbose=False):
        """
        Update the concentration to get specific total content and gradient.
        @param totContent: new total content of the element in the cylinder (cm3 water or mol solute)
        @param changeRatio: ratio between toContent and old total content of the element (for trouble shooting)
        @param gradient: old concentration gradient
        @param volumes: volume of each cell of the cylinder (cm3)
        @param isWater : water element (True) or solute (False)
        @param verbose
        """

        def print_verbose(*args):
            if verbose:
                print(rank, *args)

        def solve_concentration_matrix(matrix, aB):
            SPmatrix = sparse.csc_matrix(sparse.coo_matrix(matrix))
            return LA.spsolve(SPmatrix, aB, use_umfpack=True)

        def validate_concentration(val_new, totContent, gradient, volumes):
            """ check that tot content and gradient is correct """
            try:
                assert abs(sum(val_new * volumes) - totContent) < 1e-13
                assert (abs((val_new[1:] - val_new[:-1]) - gradient) < 1e-13).all()
            except:
                print_verbose('update_concentration error', 'val_new', val_new, 'totContent', totContent, 'gradient', gradient, 'volumes', volumes, 'changeRatio', changeRatio)
                raise Exception

        # Verbose initial information
        print_verbose('isWater', isWater, 'update_concentration error', 'totContent', totContent, 'gradient', gradient,  'volumes', volumes, 'changeRatio', changeRatio)

        # Create and solve concentration matrix
        matrix_size = self.NC - 1
        matrix = np.diag(np.full(matrix_size - 1, -1.), k=1) + np.diag(np.full(matrix_size, 1.), k=0)
        matrix[-1, :] = volumes
        aB = np.append(-gradient, totContent)
        val_new = solve_concentration_matrix(matrix, aB)
        
        print_verbose('val_new', val_new)
        validate_concentration(val_new, totContent, gradient, volumes)

        # Set bounds based on water or solute
        if isWater:
            maxVal = self.vg_soils[0].theta_S
            minVal = self.theta_wilting_point
        else:
            maxVal = np.inf
            minVal = 0.
        # check that the mean concentration respects the boundary 
        assert ( (minVal - totContent/sum(volumes)) < self.eps) 
        assert ( (totContent/sum(volumes) - maxVal) < self.eps) 
            

        # Verbose before changes
        print_verbose('BEFORE possible changes.val_new',val_new, 'new content', val_new * volumes, 'totContent', totContent, 'sum new content', sum(val_new * volumes), 'sum(val_new * volumes) - totContent', sum(val_new * volumes) - totContent, 'maxVal',maxVal, 'minVal',minVal)

        # Adapt values if necessary
        try:
            val_new = np.array(self.adapt_values(val_new_ = val_new, minVal_ = minVal, maxVal_=maxVal, volumes_=volumes, 
                                                 divideEqually_=True, verbose_ = verbose))
        except:
            print("issue adapt_values val_new, minVal, maxVal, volumes",
                  val_new, minVal, maxVal, volumes)
            raise Exception
        #val_new = helpful.adapt_values(val_new.copy(), minVal, maxVal, volumes, divideEqually=True)
        
        #print('val_new_,val_new',val_new_,val_new)
        #assert max(abs(val_new_-val_new)) < 1e-13

        # Final verbose and checks
        #print_verbose('AFTER possible changes. new content', val_new * volumes, 'totContent', totContent, 'sum new content', sum(val_new * volumes), 'change ratio error', sum(val_new * volumes) - totContent)
        
        # check that concentration in each cell respects the boundary and totcontent is still correct
        assert (val_new >= minVal).all()
        assert (val_new <= maxVal).all()
        assert abs(sum(val_new * volumes) - totContent) < self.eps*10.

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

    def changedLen(self,cyl):
        gId = cyl.gId
        return (abs(self.seg_length[gId] -
             cyl.segLength) >= 1e-15 or abs((self.seg_length[gId] -
                                             cyl.segLength) / cyl.segLength) >= 1e-5)


    def updateOld(self, lId, cyl, smaller=False, shapeOnly=False,
                  thetaLeftOver=0., konzLeftOver=None,
                  verbose=False):
        """
        Update distribution of a cylinder if its volume has changed.
        
        Parameters:
            lId: Local cylinder id.
            cyl: Cylinder object.
            smaller: Apply to the cylinder which has shrunk (True) or grown (False).
            thetaLeftOver: Mean theta in the newly available space (if smaller is False).
            konzLeftOver: Mean solute concentration in the newly available space (if smaller is False).
        """
        konzLeftOver = konzLeftOver if konzLeftOver is not None else np.zeros(self.numSoluteComp)
        gId = self.eidx[lId]
        assert gId == cyl.gId
        cellId = self.seg2cell[gId]

        segsId = np.array([ids for ids in self.cell2seg[cellId] if self.organTypes[ids] == 2])
        oldPoints = np.array(cyl.getPoints()).flatten()
        lb = self.logbase
        a_in = self.radii[gId]
        a_out = self.outer_radii[gId]
        

        if not ((self.seg2cell[gId] >= 0) and (self.organTypes[gId] == 2)):
            points = np.array([a_in,a_out ])
            assert len(points) == len(oldPoints)
        else:
            points = np.logspace(np.log(a_in) / np.log(lb), np.log(a_out) / np.log(lb), self.NC, base = lb) # cm

            volOld = cyl.getCellVolumes()
            ## new shape: 
            volNew = self.getVolumes(points, self.seg_length[gId] ) # cm^3
            deltaVol = sum(volNew) - sum(volOld) 
            
            ## shape:
            centersNew = (points[1:] + points[:-1]) / 2  # cm

                                                         
                                                   
            changeShapeAndConc = ((((deltaVol > 0) and (not smaller)) or (
                    (deltaVol < 0) and  smaller)
                          )  and (cellId in self.cellIdleftover)) and (not shapeOnly)


            changeShape = (shapeOnly and ( # only update if shape changes but no total solutes or water concentration change needed
                        ((deltaVol != 0.) and (cellId not in self.cellIdleftover)) # volume changed but it does not receive/loose solutes/water (so no total concentration change)
                        or (self.changedLen(cyl) and (deltaVol == 0.)) # length changed but not the volume (so no total concentration change)
                              )) 

            if (changeShapeAndConc or changeShape) and (not self.static_plant):


                changeRatio = 1.
                if changeShapeAndConc:
                    assert (thetaLeftOver > -self.maxdiff1d3dCNW_abs[0]) or (deltaVol <= 0.)
                                                               
                    ## change ratio
                    changeRatio = min(sum(volNew)/sum(volOld), 1.)# we migth have volNew > volOld if the gain by L increase is higher than loss via r_out decrease


                ##  water:
                theta_old = cyl.getWaterContent() # cm3/cm3
                gradientOld = (theta_old[1:] - theta_old[:-1])
                gradientNew = self.interAndExtraPolation(points[1:-1],oldPoints[1:-1], gradientOld)
                
                wOld = sum(theta_old*volOld)
                maxVal = self.vg_soils[0].theta_S  # keep that or set lower higher bound?
                minVal = self.theta_wilting_point
                
                if smaller:
                    changeRatioW = max(min(changeRatio, sum(maxVal*volNew)/wOld),sum(minVal*volNew)/wOld)
                else:
                    changeRatioW = changeRatio
                    
                try:     
                    assert min(theta_old) + self.eps >=minVal 
                    assert max(theta_old) - self.eps <=maxVal
                    theta_old = np.clip(theta_old, minVal, maxVal) # mannually fix if error <= self.eps
                    assert ((changeRatioW <= 1.) and (changeRatioW > 0.))
                    assert wOld > 0.
                    #assert max(thetaLeftOver*deltaVol,0.) >= 0. # makes no sense? (was testing precision of python?)
                except:
                    print('changeRatioW error', changeRatioW, wOld, thetaLeftOver,deltaVol, maxVal, minVal)
                    print('sum(maxVal*volNew)/wOld)',sum(maxVal*volNew)/wOld, sum(minVal*volNew)/wOld)
                    print('theta_old',min(theta_old), max(theta_old), 
                          repr(theta_old),repr(changeRatio), repr(volNew),repr(volOld) )
                    print('points',repr(points),repr(oldPoints))
                    print("cause of error",min(theta_old) + self.eps >=minVal, max(theta_old) - self.eps <=maxVal,
                          ((changeRatioW <= 1.) and (changeRatioW > 0.)),wOld > 0.,max(thetaLeftOver*deltaVol,0.) >= 0.)
                    raise Exception
                try:
                    try:
                        # get theta for each cell
                        theta_new = self.update_concentration(totContent = wOld*changeRatioW + max(thetaLeftOver*deltaVol,0.),
                                                             changeRatio=changeRatioW, 
                                                         gradient =gradientNew, volumes = volNew,isWater = True,
                                                            verbose = False)
                    except:
                        print('1st faile of update_concentration for water, try again, try again with verbose = True. rank', rank)
                        theta_new = self.update_concentration(totContent = wOld*changeRatioW + max(thetaLeftOver*deltaVol,0.),
                                                             changeRatio=changeRatioW, 
                                                         gradient =gradientNew, volumes = volNew,isWater = True,
                                                            verbose = True)
                    newHead = np.array([vg.pressure_head(nt, self.vg_soils[0]) for nt in theta_new])# cm
                except:
                         
                    print('\t',gId,"x_old",theta_old,theta_old* volOld,volOld )
                    print('\t',gId,"points",oldPoints, points,centersNew)
                    print('\t',gId,"gradient",gradientOld, gradientNew)
                    print('\t',gId,"vg param", self.vg_soils[0].theta_R, self.vg_soils[0].theta_S,self.theta_wilting_point)   
                    print('\t',gId,"changeRatio",changeRatio, changeRatioW)    
                    #print('\t',gId,"newHead", newHead ,theta_new)
                    raise Exception
                    
                ## new contents:  
                molFrOld =np.array( [np.array(cyl.getSolution(nC+1)) for nC in range(self.numSoluteComp)])   #mol/mol , g/g
                volWatNew = theta_new *volNew
                molFrNew = []
                for nComp in range(1, self.numComp):
                    if (molFrOld[nComp -1] != 0.).any():
                        isDissolved = (nComp <= self.numDissolvedSoluteComp)
                        if isDissolved: # mol phase = [cm3 phase] * [mol phase /cm3 phase] 
                            molarPhaseOld = theta_old*volOld * self.phaseDensity(nComp ) 
                            molarPhaseNew = volWatNew * self.phaseDensity(nComp ) # 
                        else:
                            molarPhaseOld = volOld * self.phaseDensity(nComp )
                            molarPhaseNew = volNew * self.phaseDensity(nComp ) 
                        gradientOld = (molFrOld[nComp -1][1:] - molFrOld[nComp -1][:-1])   
                        gradientNew = self.interAndExtraPolation(points[1:-1],oldPoints[1:-1], gradientOld)
                        cOld = sum(molFrOld[nComp -1] *   molarPhaseOld  )    
                        assert abs(cOld - sum(cyl.getContent(nComp)))< 1e-13
                        
                        try:
                            try:
                                #get mole fraction for each cell
                                molFrNew.append(
                                    self.update_concentration(totContent = cOld*changeRatio+ max(konzLeftOver[nComp -1]*deltaVol,0.),
                                                                         changeRatio=changeRatio,gradient =gradientNew, 
                                        volumes = molarPhaseNew,isWater = False, verbose = False))                             
                            except:
                                print('1st faile of update_concentration for comp no',nComp,', try again, try again with verbose = True. rank', rank,
                                     'max(konzLeftOver[nComp -1]*deltaVol,0.)',max(konzLeftOver[nComp -1]*deltaVol,0.))
                                molFrNew.append(
                                    self.update_concentration(totContent = cOld*changeRatio + max(konzLeftOver[nComp -1]*deltaVol,0.),
                                                                         changeRatio=changeRatio,gradient =gradientNew, 
                                        volumes = molarPhaseNew,isWater = False, verbose = True)) 
                                        
                        except:
                            print('update_concentration failed','cOld',cOld,#'contentC[nComp-1][cellId] ',contentC[nComp-1][cellId] ,
                                  'changeRatio',changeRatio,'gradientNew',gradientNew,'molarPhaseNew',
                                 molarPhaseNew,'molFrOld[nComp -1]',molFrOld[nComp -1],'molFrOld',molFrOld[nComp -1] ,
                                  'molarPhaseOld',   molarPhaseOld,'sum(volNew)',sum(volNew),
                                  'sum(volOld)',sum(volOld) )
                            raise Exception
                        # final checks that the content is as expected after the updates
                        try:
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
                        
                # create new cylinder from the data
                self.cyls[lId] = self.initialize_dumux_nc_( gId, 
                                                            x = newHead,# cm
                                                            cAll = molFrNew, # mol/mol water or mol/mol scv
                                                            Cells = centersNew) # cm
            
                
                
            
    def get_Vol_leftoverI(self, idCell):# cm3
        """ difference between volume of cylinder and volume of 3d soil voxel  """
        idCylsAll = 0
        idCyls =[]
        idcs = self.getCellIds()


        vol_rhizo = 0
        if (idCell in self.cell2seg) and (idCell in self.cellIdleftover):
            idCylsAll, idCyls =   self.getIdCyllMPI(idCell, True, doSum = True)

            # we have old segments and leftover space (to avoid division by 0 later)
            if idCylsAll > 0:
                vol_rhizo= self.getVolumesCyl(idCyls=idCyls,idCylsAll=idCylsAll) #m3

        if rank == 0:
            newVol = self.soilModel_cellVolumes[idCell] - vol_rhizo
            try:
                assert newVol > -1e-10
            except:
                print('issue in get_Vol_leftoverI',rank, idCylsAll, 
                      len(np.array(self.outer_radii)), np.array(self.outer_radii))
                print( idCell,newVol,self.soilModel_cellVolumes[idCell], vol_rhizo)
                raise Exception
        else:
            newVol =0
            
        return max(newVol,0.)
        
    def get_watVol_leftoverI(self, idCell, wat_total_):# cm3
        """ difference between water volume of cylinder and volume of 3d soil voxel  
        """
        
        wat_rhizo = 0
        idSegs = []
        idCyls = []
        idCylsAll = 0
        if (idCell in self.cell2seg) and (idCell in self.cellIdleftover):
            idSegs = self.cell2seg[idCell]#all segments in the cell
            idSegs=np.array( [ids for ids in idSegs if (self.organTypes[ids] == 2)])#only root segments
            idCylsAll, idCyls =   self.getIdCyllMPI(idCell, True, doSum = True)   # id of segments which already have a 1d model
            if idCylsAll > 0: # we have old segments
                wat_rhizo= self.getWaterVolumesCyl(idCyls) #cm3
        if rank == 0:
            wat_total = wat_total_[idCell] * self.soilModel_cellVolumes[idCell] # [cm3] water contentw
            watVol_res = wat_total - wat_rhizo #cm3
            if watVol_res < 0. and self.fixed_point_iteration:            
                try:
                    # no new segments (watVol_res should be 0) and error >= max water error -1e-13
                    assert (len(idSegs) == idCylsAll) and (watVol_res > (-1e-13-self.maxdiff1d3dCNW_abs[0])) # rounding error probably
                    watVol_res = 0.
                except:
                    print("getWatPsi_leftoverI")
                    print(idCell, wat_total,wat_rhizo, wat_total - wat_rhizo,self.maxdiff1d3dCNW_abs[0])
                    print( len(idSegs), len(idCyls), idCylsAll)# how many old and new cylinders?
                    raise Exception
        else:
            watVol_res =0
                
        return watVol_res
        
    
    def getC_content_leftoverI(self, idCell, idComp, mol_total):# mol
        """ difference between C content of cylinder and of 3d soil voxel  
            @param: idCell: id of the cell
            @param: idComp: solute component id
            @param: mol_total: C content of component idComp according to 3d soil models
        
        """
        
        mol_rhizo = 0
                      
        idCyls = []
        idCylsAll = 0
        if (idCell in self.cell2seg) and (idCell in self.cellIdleftover):
            idSegs = self.cell2seg[idCell]#all segments in the cell
            idSegs=np.array( [ids for ids in idSegs if (self.organTypes[ids] == 2)])#only root segments
            idCylsAll, idCyls =  self.getIdCyllMPI(idCell, True, doSum = True)   
            
            if idCylsAll > 0:#we have old segments
                mol_rhizo = self.getContentCyl(idCyls, idComp=idComp, doSum = True)#adapted for MPI
        
        if rank == 0:
            mol_total = mol_total[idCell] # mol
            res_CC = mol_total - mol_rhizo
            
            if res_CC < 0. and self.fixed_point_iteration:
                if (res_CC >(-1e-13-self.maxdiff1d3dCNW_abs[idComp]*10)):# and (len(idSegs) == len(idCylsAll))): # rounding error probably 
                    res_CC = 0.
                else:
                    print("getC_content_leftoverI", idCell)
                    print("res_CC = ",res_CC," < ",(-self.maxdiff1d3dCNW_abs[idComp]*10),"or",len(idSegs),"=/=",idCylsAll,
                    ", idComp:",idComp,'mol_total',mol_total ,
                    'mol_rhizo', mol_rhizo,'self.maxdiff1d3dCNW_abs',self.maxdiff1d3dCNW_abs ) 
                    print('lens', len(idSegs),idCylsAll, len(idCyls))# how many old and new cylinders?
                    raise Exception
        else:
            res_CC =0
                
        return res_CC

    

                
    def initialize_(self,gId,x,cc ):
        """ create new segments for soil element content data
            @param gId: global segment index
            @param x: mean water potential per 3d soil cell (cm)
            @param cc: solute molar fraction per 3d soil cell (mol/mol)
        
        """
        if ((self.seg2cell[gId] >= 0) and (self.organTypes[gId] == 2)):
            try:
                self.cyls.append(self.initialize_dumux_nc_(gId, x = x[self.seg2cell[gId]], cAll = cc[self.seg2cell[gId]]))
            except:
                print('error at initialization',gId,self.seg2cell[gId],
                      self.organTypes[gId],
                    self.radii[gId], self.outer_radii[gId],self.seg_length[gId],
                    x[self.seg2cell[gId]],
                      cc[self.seg2cell[gId]] )
                raise Exception
        else:
            # shoot or root aboveground
            a_in = self.radii[gId]#get Perimeter instead? not used for now anyway
            a_out = self.outer_radii[gId]
            self.cyls.append(AirSegment(a_in, a_out,self.seg_length[gId],self.numSoluteComp)) #psi_air computed by the photosynthesis module.
            self.cyls[-1].gId = gId    
            self.cyls[-1].segLength = self.seg_length[gId]    
    

    def initialize_dumux_nc_(self, gId, x,                          # cm
                                    cAll = [0. for i in range(22)], # mol/mol wat or mol/mol scv
                                         Cells = []):               # cm
        verbose = False
        a_in = self.radii[gId]
        a_out = self.outer_radii[gId]
        lId =int( np.where(self.eidx == gId)[0])
                    
        if a_in < a_out:
            cyl = RichardsNoMPIWrapper(self.RichardsNCCylFoam())  # only works for RichardsCylFoam compiled without MPI
            cyl.results_dir = self.soilModel.results_dir            
            cyl.gId = gId    
            cyl.theta_wilting_point = self.theta_wilting_point    
            cyl.a_in = a_in    
            cyl.a_out = a_out
            cyl.MaxRelativeShift = self.MaxRelativeShift_1DS  
            helpful.setDefault(cyl) 
            cyl.initialize(verbose = False)
            lb = self.logbase
            
            points = np.logspace(np.log(a_in) / np.log(lb), np.log(a_out) / np.log(lb), 
                                 self.NC, base = lb)
            
            cyl.createGrid1d(points, self.seg_length[gId] )# cm
                
            if len(self.dx2) > lId:
                self.dx2[lId] = 0.5 * (points[1] - points[0]) #when updating
            else:
                self.dx2.append(0.5 * (points[1] - points[0]))
                
                    
            cyl.setOuterBC("fluxCyl", 0.)  #  [cm/day] Neumann boundary condition
            cyl.setInnerBC("fluxCyl", 0.) 

            cyl.setParameter("Problem.segLength", str(self.seg_length[gId]))  # cm      
            cyl.segLength = self.seg_length[gId]              
            cyl.setParameter("Soil.CriticalPressure", str(self.soilModel.wilting_point))            
            cyl.setICZ(x, list(Cells))#setParameter("Soil.IC.P", cyl.dumux_str(x))# cm
            for j in range( 1, self.numComp):  
                cyl.setICZ_solute(cAll[j-1], list(Cells))
            if len(Cells) > 0:
                if len(Cells)!= len(x):
                    print("Cells, x",Cells, x, len(Cells), len(x))
                    raise Exception              
                assert(len(cAll[j-1])==len(Cells))
                    
            self.setSoilParam(cyl)
            
            cyl.maxDt =  250./(3600.*24.)
            cyl.initializeProblem(cyl.maxDt)
            cyl.ddt = 1.e-3
            
            cyl.setCriticalPressure(cyl.wilting_point)  
            
            ThetaCyl = cyl.getWaterContent()
            cyl.createNewtonSolver() # make sure solver parameters are implemented.

            try:
                assert (ThetaCyl >= self.vg_soils[0].theta_R).all()
                assert (ThetaCyl <= self.vg_soils[0].theta_S).all()
            except:
                print('issue thetaCyl',rank,ThetaCyl, self.vg_soils[0].theta_R, self.vg_soils[0].theta_S )
                raise Exception
            pHeadcyl = cyl.getSolutionHead()
            
            if verbose:
                print('end initialize_',gId,self.seg2cell[gId],'wat vol?',sum(cyl.getWaterVolumes()),self.seg_length[gId],
                  pHeadcyl , x,pHeadcyl - x,'Cells',Cells, a_in, a_out, cyl.segLength)
            
                
            try:
                x_divide = np.where(np.array(x)!=0,x,1)
                thetainit = cyl.getWaterContent()
                thetainitth = np.array([vg.water_content( p_mean_, self.vg_soils[0]) for
                                        p_mean_ in pHeadcyl])
                theta_divide = np.where(thetainitth!=0,thetainitth,1)

                assert abs(self.seg_length[gId] -
                           cyl.segLength) < 1e-16 or abs((self.seg_length[gId] -
                           cyl.segLength)/cyl.segLength) < 1e-5
                # we could have some issue because of the water
                # retention curve regularisation in dumux
                assert len(pHeadcyl) == (self.NC - 1)
                assert (np.logical_or( (abs((pHeadcyl - x)/x_divide)*100 < 1e-5) , 
                                       (abs(pHeadcyl - x) < 1e-9) )).all()
                assert (np.logical_or( (abs((thetainit - thetainitth)/theta_divide)*100 < 1e-5) ,
                                       (abs(thetainit- thetainitth) < 1e-9) )).all()
            except:
                print('error: issue with cylinder creations', rank)
                print('len(pHeadcyl) == (self.NC - 1)?', len(pHeadcyl), (self.NC - 1))
                print('(abs((pHeadcyl - x)/x_divide)*100 > 1e-5).all()?',abs((pHeadcyl - x)/x_divide)*100, pHeadcyl,x,'x_divide',x_divide)
                print('theta',thetainit,thetainitth)
                raise Exception
            ##
            # PIFidx : index or element where we will have a proposed flux from the plant model.
            # ==> all the water elements except if we have a michaelis menten uptake (BC idx == 8)
            # TODO: update for the 10c and 22c model
            ##
            cyl.typeBC = int(cyl.getParameter(cyl.param_group + "BC.Bot.SType"))
            if cyl.typeBC == 8:
                self.PIFidx = []
            else: 
                self.PIFidx = self.RIFidx.copy()
                
            return cyl
            
        else:
            print("RhizoMappedSegments.initialize_dumux_: Warning, segment {:g} might not be in domain, radii [{:g}, {:g}] cm".format(gId, a_in, a_out))
            raise Exception
            return []



    def get_inner_heads(self,  weather:dict={}):#        
        """ matric potential at the root surface interface [cm]"""
        
        RH = weather['RH'] if 'RH' in weather.keys() else weather["ea"]/weather["es"]
        TairC = weather['TairC'] if 'TairC' in weather.keys() else weather["Tair"]
        rsx = np.array([
            cyl.getInnerHead() if not isinstance(cyl, AirSegment) else helpful.getPsiAir(RH,TairC) for cyl in self.cyls
        ])  # [cm] (in richards.py, then richards_cyl.hh)
        inner_heads = self.gather( rsx)
        return inner_heads

    def get_inner_solutes(self, compId = 1):
        """ matric potential at the root surface interface [mol/cm3]"""
        rsx = np.full(len(self.cyls),0.)
        for i, cyl in enumerate(self.cyls):  # run cylindrical models
            rsx[i] = cyl.getInnerSolutes( compId=compId)  # [cm]
        inner_solutes= self.gather( rsx)#_flat0(comm.) #self._map()   # gathers and maps correctly
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
        inner_concentrations = self.gather( rsx)#_flat0(comm.) #self._map()  # gathers and maps correctly
        return inner_concentrations


    def get_dx2(self):
        """ TODO doc me AND only for mode="dumux" yet (set in initialize)"""
        dx2 = self.gather( self.dx2)#_flat0(comm) #self._map() 
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
                    if (len(isAirSegs)>1) and (rank==0):
                        V_rhizo[isAirSegs] = 1 # avoid RuntimeWarning: invalid value encountered in true_divide                    
                else:
                    V_rhizo = self.getVolumesCyl(doSum = False, reOrder = True) #[cm3]                
                contentRhizo = self.getContentCyl(idComp=compId, doSum = False, reOrder = True)# mol
            if (len(isAirSegs)>1) and (rank==0):
                assert (contentRhizo[isAirSegs] == 0.).all() #no solutes in the space around the air segments

            if not concentration:
                return contentRhizo
            concentration = contentRhizo/V_rhizo
            if (len(isAirSegs)>1) and (rank==0):
                concentration[isAirSegs] = np.nan
            return concentration # mol/cm3
        else:
            raise Exception("RhizoMappedSegments.get_inner_concentrations: Warning, mode {:s} unknown".format(self.mode))
        #concentration = self.gather(rsx,root=0)#self._flat0() #self._map( )
        #return concentration

        
    def _calculate_and_set_initial_solute_flows(self,cyl, dt,inner_fluxes_CN_, outer_fluxes_CN_,verbose):
        """Retrieve solute values for the segment.
            inner_fluxes_CN_ : net plant small molecule releases,  [mol day-1] 
            outer_fluxes_CN_: 1d-3d or 1d-1d small molecule exchanges,  [mol day-1] 
            returns
            inner_fluxes_CN: net plant molecule releases,  [mol day-1] 
            outer_fluxes_CN: net soil molecule exchanges,  [mol day-1] 
        """
        gId = cyl.gId
        l = cyl.segLength
        
        # from list to float
        inner_fluxes_CN = np.zeros(self.numDissolvedSoluteComp) 
        inner_fluxes_CN[self.PIFidx] = np.array([ 
            inner_fluxes_CN_[ii][gId] if not isinstance(inner_fluxes_CN_[ii], numbers.Number) else inner_fluxes_CN_[ii]
            for ii in range(len(self.PIFidx))])
        outer_fluxes_CN = np.array([ 
            outer_fluxes_CN_[jj][gId] if not isinstance(outer_fluxes_CN_[jj], numbers.Number) else outer_fluxes_CN_[jj]
            for jj in range(self.numDissolvedSoluteComp)])
           
        if self.do1d1dFlow:
            outer_fluxes_CN += self.flow1d1d_sol[gId]
            #outer_fluxes_mucil += self.flow1d1d_mucil[gId]
        

        
                  
        qIn_solMucil = np.full(self.numSoluteComp, 0.)
        for jj in range(self.numDissolvedSoluteComp):
            qIn_solMucil[jj] = inner_fluxes_CN[jj] / (2 * np.pi * self.radii[gId] * l)
        
        cyl.setSoluteBotBC(cyl.typeBC, qIn_solMucil)
        
        QflowOutCellLim = cyl.distributeSources(dt, source = outer_fluxes_CN,
                               inner_fluxs=inner_fluxes_CN, 
                               eqIdx =  np.array([nc+1 for nc in range(self.numDissolvedSoluteComp)]),plantM=self)
        
        # reset after first limitation in _distribute_source
        outer_fluxes_CNbu = outer_fluxes_CN
        outer_fluxes_CN = np.array([sum(QflowOutCellLim[nc]) for nc in range(self.numDissolvedSoluteComp)])
        
        for nc in range(self.numDissolvedSoluteComp):
            divideQout = outer_fluxes_CNbu[nc] if outer_fluxes_CNbu[nc] != 0 else 1

            #print('here', outer_fluxes_CN[nc],outer_fluxes_CNbu[nc])
            test = abs((outer_fluxes_CN[nc] - (outer_fluxes_CNbu[nc]))/ divideQout)*100

            #if verbose or (test > 1.):
            if (test > 1.):
                print('rhizo_modelsPlant::_calculate_and_set_initial_solute_flows, gId',gId,
                      ' qIn',qIn_solMucil, 
                inner_fluxes_CN,'qout', QflowOutCellLim,'outer_fluxes_CN',
                outer_fluxes_CN,
                      'qout init',outer_fluxes_CNbu,'divideQout',divideQout,
                      'shape', self.radii[gId],  
                      np.array( self.outer_radii)[gId] , 'length', cyl.segLength,'nc',nc,self.plant_or_RS,test )
                raise Exception
        return  inner_fluxes_CN, outer_fluxes_CN                 
                     
        
        
    def _calculate_and_set_initial_water_flows(self,cyl,dt,inner_fluxes, outer_fluxes, verbose):
        """ set water bc and sinks for cyl and solves it 
            inner_fluxes: PWU,  [cm3 day-1] 
            outer_fluxes: 1d-3d or 1d-1d water exchanges,  [cm3 day-1]  
            
            returns
            inner_fluxes: initial PWU [cm3 day-1] 
            outer_fluxes: initial (limited) soil water exchange [cm day-1] 
        """
        l = cyl.segLength
        gId = cyl.gId
        
        # PWU
        inner_fluxes = inner_fluxes[gId]
        qIn = inner_fluxes / (2 * np.pi * self.radii[gId] * l)
        cyl.setInnerFluxCyl(qIn)
        
        # 1d-1d and 1d-3d soil water exchange
        outer_fluxesbu = outer_fluxes[gId]   
        outer_fluxes = outer_fluxes[gId]    
        if self.do1d1dFlow:
            outer_fluxes += self.flow1d1d_w[gId]
            outer_fluxesbu += self.flow1d1d_w[gId]
        # might have QflowOutLim != outer_fluxes if abs(outer_fluxes) is too high
        QflowOutCellLim = cyl.distributeSource(dt, outer_fluxes,inner_fluxes*dt, eqIdx=0,plantM=self)
        outer_fluxes = sum(QflowOutCellLim)  # get outer_fluxes limited  
        divideQout = outer_fluxesbu if outer_fluxesbu != 0 else 1
        if verbose or ((abs((outer_fluxes - (outer_fluxesbu))/ divideQout)*100 > 1.) and self.fixed_point_iteration):
            print('rhizo_modelsPlant::_calculate_and_set_initial_water_flows, gId',gId,' qIn',qIn, 
            inner_fluxes,'qout', outer_fluxes,'qout init',outer_fluxesbu,'shape', self.radii[gId] , 'length', l )
            raise Exception
            
        return inner_fluxes, outer_fluxes
           
            
    
    
    def solve(self, dt, *argv):
        """ set bc and sinks for cyl and solves it 
            proposed_inner_fluxes: PWU,  [cm3 day-1] 
            proposed_outer_fluxes: 1d-3d or 1d-1d water exchanges,  [mol C day-1] 
            proposed_inner_fluxes_sol : net plant small C molecule releases,  [mol C day-1] 
            proposed_outer_fluxes_sol: 1d-3d or 1d-1d small C molecule exchanges,  [mol C day-1] 
            proposed_inner_fluxes_mucil : net plant large C molecule releases,  [mol C day-1] 
            proposed_outer_fluxes_mucil: 1d-3d or 1d-1d large C molecul exchanges,  [mol C day-1] 
        """
        #inner = bot = plant
        #outer = top = soil
        if rank == 0:
            if len(argv) >= 4:
                proposed_inner_fluxes = argv[0]
                proposed_outer_fluxes = argv[1]
                proposed_inner_fluxes_sol = argv[2]
                proposed_outer_fluxes_sol = argv[3]     
            else:       
                proposed_inner_fluxes = argv[0]
                proposed_outer_fluxes = argv[1]
                proposed_outer_fluxes_sol = argv[2]                   
                proposed_inner_fluxes_sol = np.zeros(proposed_outer_fluxes_sol.shape)                
        else:
            proposed_inner_fluxes = None
            proposed_outer_fluxes = None
            proposed_inner_fluxes_sol = None
            proposed_outer_fluxes_sol = None
        
        proposed_inner_fluxes = comm.bcast(proposed_inner_fluxes, root=0)
        proposed_outer_fluxes = comm.bcast(proposed_outer_fluxes, root=0)
        proposed_inner_fluxes_sol = comm.bcast(proposed_inner_fluxes_sol, root=0)
        proposed_outer_fluxes_sol = comm.bcast(proposed_outer_fluxes_sol, root=0)
        
        if isinstance(proposed_inner_fluxes_sol[0], (int,float)):
            proposed_inner_fluxes_sol = [proposed_inner_fluxes_sol]
            proposed_outer_fluxes_sol = [proposed_outer_fluxes_sol]
        
        # initialize_flux_storage
        self.seg_fluxes_limited = np.full(len(self.cyls), 0.)
        self.seg_fluxes_limited_Out = np.full(len(self.cyls), 0.)
        # sources/sink caused by Q_3D soil
        self.seg_fluxes_limited_CN_Out = np.full((self.numDissolvedSoluteComp,len(self.cyls)), 0.)
        # plant uptake/release        
        self.seg_fluxes_limited_CN_In = np.full((len(self.RIFidx),len(self.cyls)), 0.)
        self.seg_sources_limited = np.full((self.numComp,len(self.cyls)), 0.)
        
        for lId, cyl in enumerate(self.cyls):  # run cylindrical models
            gId = cyl.gId # for one process gId == lId
            assert gId == self.eidx[lId]  
            verbose = False
            #if True:#verbose:
            #    print(rank, "cyl no ",lId+1,"/",len(self.cyls),'gId')
            #n = round(float(lId+1) / float(len(self.cyls)) * 100.) 
            #print("\t1D models: [" + ''.join(["*"]) * n + ''.join([" "]) * (100 - n) + "]")

                
            if isinstance(cyl, AirSegment):  
                self._handle_air_segment(cyl, lId, proposed_inner_fluxes,
                                    proposed_outer_fluxes,
                                    proposed_inner_fluxes_sol,
                                    proposed_outer_fluxes_sol)                
            else:
                self._handle_non_air_segment(cyl, lId, dt, proposed_inner_fluxes,
                                    proposed_outer_fluxes,
                                    proposed_inner_fluxes_sol,
                                    proposed_outer_fluxes_sol, verbose)
  

    
    def _handle_air_segment(self, cyl, lId, inner_fluxes, outer_fluxes, inner_fluxes_CN_, outer_fluxes_sol):
        """Handle the case where the segment is an AirSegment."""
        cyl.setInnerFluxCyl(inner_fluxes[cyl.gId])
        self.seg_fluxes_limited[lId] = inner_fluxes[cyl.gId]
        self.seg_fluxes_limited_Out[lId] = outer_fluxes[cyl.gId]
        for jj in range(self.numDissolvedSoluteComp):
            self.seg_fluxes_limited_CN_Out[jj][lId] = outer_fluxes_sol[jj][cyl.gId]
        for ii, jj in enumerate( self.PIFidx ):
            self.seg_fluxes_limited_CN_In[ii][lId] = inner_fluxes_CN_[ii][cyl.gId] if not isinstance(inner_fluxes_CN_[ii], numbers.Number) else inner_fluxes_CN_[ii]

        
            
    def _handle_non_air_segment(self, cyl,lId, dt, inner_fluxes, outer_fluxes, inner_fluxes_sol,
            outer_fluxes_sol, verbose = False):
        """ set bc and sinks for cyl and solves it 
            proposed_inner_fluxes: PWU,  [cm3 day-1] 
            proposed_outer_fluxes: 1d-3d or 1d-1d water exchanges,  [mol C day-1] 
            proposed_inner_fluxes_sol : net plant C molecule releases,  [mol C day-1] 
            proposed_outer_fluxes_sol: 1d-3d or 1d-1d C molecule exchanges,  [mol C day-1] 
        """

        inner_fluxes_water, outer_fluxes_water = self._calculate_and_set_initial_water_flows(
                                                    cyl,dt, inner_fluxes, outer_fluxes, 
                                                        verbose)         
        inner_fluxes_CN, outer_fluxes_CN = self._calculate_and_set_initial_solute_flows(cyl, dt, 
                                                       inner_fluxes_sol, outer_fluxes_sol, verbose)    
        
        self._run_solver(cyl, dt, lId,  inner_fluxes_water, 
                        outer_fluxes_water,  inner_fluxes_CN, outer_fluxes_CN, verbose)
    
    def _run_solver(self, cyl, dt, lId,  inner_fluxes_water,  outer_fluxes_water,   inner_fluxes_solMucil, outer_fluxes_solMucil, verbose):
        """Run the solver for the segment.
            cyl: 1ds object
            dt: time step
            lId: cylinder index on this thread (local index)
            inner_fluxes_water: PWU,  [cm3 day-1] 
            outer_fluxes_water: 1d-3d or 1d-1d water exchanges,  [cm3 day-1] 
            inner_fluxes_solMucil : net plant C molecule releases,  [mol C day-1] 
            outer_fluxes_solMucil: 1d-3d or 1d-1d C molecule exchanges,  [mol C day-1] 
        """

        #return inner_fluxes_real[0], outer_fluxes_water_temp, inner_fluxes_real[1:], outer_fluxes_solMucil_temp, cyl_sources       
        #return inner_fluxes_real, outer_fluxes_water_temp, outer_fluxes_solMucil_temp, cyl_sources 
        #inner_fluxes_water, outer_fluxes_water,   inner_fluxes_solMucil, outer_fluxes_solMucil

        inner_fluxes, outer_fluxes_water,outer_fluxes_solMucil, cyl_sources = self._attempt_solution(cyl, dt, 
                                                    lId, inner_fluxes_water, 
                                                    outer_fluxes_water,   inner_fluxes_solMucil, 
                                                    outer_fluxes_solMucil, verbose)
        # store boundary fluxes for mass balance check
        self.seg_fluxes_limited[lId] = inner_fluxes[0]
        self.seg_fluxes_limited_Out[lId] = outer_fluxes_water
        for jj in range(self.numDissolvedSoluteComp):
            self.seg_fluxes_limited_CN_Out[jj][lId] = outer_fluxes_solMucil[jj]
        for ii, jj in enumerate( self.RIFidx ):
            self.seg_fluxes_limited_CN_In[ii][lId] = inner_fluxes[jj+1] 
        for ii in range( self.soilModel.numComp ):
            self.seg_sources_limited[ii][lId] = cyl_sources[ii]
            
    def _check_Ccontent(self, cyl):
        for ncomp in range(self.numSoluteComp):
            try:
                assert (np.array(cyl.getSolution(ncomp + 1)).flatten() >= 0).all()
            except:
                print(rank, '(np.array(cyl.getSolution(ncomp + 1)).flatten() < 0).any()', 
                      ncomp, np.array(cyl.getSolution(ncomp + 1)).flatten(), 
                      (np.array(cyl.getSolution(ncomp + 1)).flatten() >= 0).all())
                return  False
        return True
        
        
    def _check_Wcontent(self, cyl):
        try:
            maxVal = self.vg_soils[0].theta_S
            minVal = self.theta_wilting_point
            theta = cyl.getWaterContent()
            phead = cyl.getSolutionHead()
            assert min(theta)>= minVal
            assert max(theta)<= maxVal
            assert max(phead)<= 0.
        except:
            print(rank, 'min(theta)< minVal or max(theta)>maxVal', 
                  cyl.gId, 
                  'or phead > 0',phead)
            return  False
        return True
        
    def _reset_newton_solver(self, cyl,MaxRelativeShift=None,
                                EnableResidualCriterion=None, 
                                EnableAbsoluteResidualCriterion=None,
                                SatisfyResidualAndShiftCriterion=None, 
                                max_steps=None, 
                                max_divisions=None):
        if MaxRelativeShift is not None:
            cyl.setParameter("Newton.MaxRelativeShift", str(MaxRelativeShift))
        if EnableResidualCriterion is not None:
            cyl.setParameter("Newton.EnableResidualCriterion", str(EnableResidualCriterion)) 
        if EnableAbsoluteResidualCriterion is not None:
            cyl.setParameter("Newton.EnableAbsoluteResidualCriterion", str(EnableAbsoluteResidualCriterion))
        if SatisfyResidualAndShiftCriterion is not None:
            cyl.setParameter("Newton.SatisfyResidualAndShiftCriterion", str(SatisfyResidualAndShiftCriterion))
        if max_steps is not None:
            cyl.setParameter("Newton.MaxSteps", str(max_steps))
        if max_divisions is not None:
            cyl.setParameter("Newton.MaxTimeStepDivisions", str(max_divisions))
        cyl.createNewtonSolver()
                
    def _attempt_solution(self, cyl, dt: float,  lId: int, inner_fluxes_water_temp, outer_fluxes_water_temp, 
                        inner_fluxes_solMucil_temp, outer_fluxes_solMucil_temp,   verbose):
        """
        Attempt the solution of the segment, retrying with different inputs if necessary.

            cyl: 1ds object
            dt: time step
            lId: cylinder index on this thread (local index)
            inner_fluxes_water: PWU,  [cm3 day-1] 
            outer_fluxes_water_temp: 1d-3d or 1d-1d water exchanges,  [cm3 day-1] 
            inner_fluxes_solMucil_temp : net plant C molecule releases,  [mol C day-1] 
            outer_fluxes_solMucil_temp: 1d-3d or 1d-1d C molecule exchanges,  [mol C day-1] 
            
            return
            inner_fluxes_real: inner BC as implemented in dumux [cm3 day-1]  or [mol C day-1] 
        """          
                
        maxRelShift = float(cyl.MaxRelativeShift)
        
        # use first ddt during last simulation
        #if len(cyl.base.BC_ddt)>0:
        #    cyl.ddt = cyl.base.BC_ddt[0]/24./60./60.
            
        initialDdt = float(cyl.ddt)
        gId = cyl.gId
        a_in = self.radii[gId]
        a_out = np.array( self.outer_radii)[gId]
        l = cyl.segLength
        
        if self.changedLen(cyl):
            print('wrong segLengthG',gId,str(self.seg_length[gId]),self.seg_length[gId], cyl.segLength,
                  (self.seg_length[gId] -  cyl.segLength),
                 (self.seg_length[gId] -  cyl.segLength) / cyl.segLength)
            raise Exception
        assert self._check_Ccontent(cyl)
        assert self._check_Wcontent(cyl)


        redoSolve = True
        n_iter_solve = 0
        verbose = False
        if verbose:
            oldWaterContent = cyl.getWaterContent()
            print('start solve, water in ',inner_fluxes_water_temp,'water out', outer_fluxes_water_temp, 
                        repr(oldWaterContent),
                        repr(cyl.getSolutionHead()),
                  'shape', repr(cyl.getPoints()),
                  repr( cyl.getCellVolumes()), a_in, a_out, l, 'dt', dt)
        inner_fluxes_real = np.full(self.numFluidComp, np.nan)
        outer_fluxes_real = np.full(self.numFluidComp, np.nan)
        errorWOnly = np.nan; errorCOnly = np.nan
        while redoSolve:
            
            try:
                #if True:
                # reseting cyl.ddt to 1.e-5 was actually creating problems for sand VG
                cyl.ddt = 1e-3#min( 1.e-5,cyl.ddt)#or just reset to 1e-5?
                #cyl.setMaxTimeStepSize(self.soilModel.maxDt_1DS/10.)
                
                didReset = False
                #valbefore = cyl.getTotCNContent_each().sum(axis = 1)
                
                cyl.solve(dt)
                #valafter = cyl.getTotCNContent_each().sum(axis = 1)
                #helpful.run_with_timeout(60., cyl.solve, dt) # after Xmn time out
                # error, time out error
                
                # cm3 water or  mol C /cm3
                inner_fluxes_real, outer_fluxes_real = cyl.getBCFlows()  
                inner_fluxes_real = -inner_fluxes_real.flatten()   
                inner_fluxes_real_ = np.zeros(self.numComp)
                inner_fluxes_real_[:self.numFluidComp] = inner_fluxes_real
                
                cyl_sources = np.array([cyl.getSource(nc) * cyl.getCellVolumes() * dt for nc in range(cyl.numComp)]).sum(axis=1)    
                
                assert (outer_fluxes_real == 0.).all() # outer exchange implemented via sinks not via outer BC
                     
                errorCOnly = not self._check_Ccontent(cyl)
                errorWOnly = not self._check_Wcontent(cyl)
                        
                assert not errorCOnly
                assert not errorWOnly
                redoSolve = False
                
                # newton parameters are re-read at each 'createNewtonSolver' calls
                self._reset_newton_solver( cyl,
                                MaxRelativeShift = cyl.MaxRelativeShift,
                                EnableResidualCriterion=cyl.EnableResidualCriterion,#"true", 
                                EnableAbsoluteResidualCriterion=cyl.EnableAbsoluteResidualCriterion,#"true",
                                SatisfyResidualAndShiftCriterion=cyl.SatisfyResidualAndShiftCriterion,#"true", 
                                max_steps=cyl.MaxSteps, #50
                                          max_divisions=cyl.MaxTimeStepDivisions#20
                                         )  
            except Exception as e:
                #if False:
                np.set_printoptions(precision=20)
                
                
                print('solve Failed:', e,'rank',rank,'gId',gId,'lId',lId,
                      'n_iter_solve',n_iter_solve)

                cyl.setParameter("Problem.verbose", "0")
                cyl.setParameter("Newton.Verbosity", "0")
                # newton parameters are re-read at each 'createNewtonSolver' calls
                self._reset_newton_solver( cyl,
                                EnableResidualCriterion=cyl.EnableResidualCriterion,#"true", 
                                EnableAbsoluteResidualCriterion=cyl.EnableAbsoluteResidualCriterion,#"true",
                                SatisfyResidualAndShiftCriterion=cyl.SatisfyResidualAndShiftCriterion,#"true", 
                                          max_divisions=cyl.MaxTimeStepDivisions#20
                                         )
                
                if n_iter_solve == 0:
                    # just reset and see if the updated cyl.ddt is enough to solve the issue
                    # get info for debug and reset
                    print(rank,'dt',dt,'initialDdt',initialDdt,
                          'errorWOnly',errorWOnly,'errorCOnly',errorCOnly,
                      'QflowIn',inner_fluxes_water_temp, 
                      ' Qflowout',outer_fluxes_water_temp,'inner_fluxes_real',inner_fluxes_real)
                    print('inner_fluxes_solMucil_temp',
                      inner_fluxes_solMucil_temp, 'outer_fluxes_solMucil_temp',outer_fluxes_solMucil_temp)
                    
                    print('plant point',repr(cyl.getPoints()), 'length',cyl.segLength)
                    cyl.reset();
                    print('pheadold',cyl.getSolutionHead())
                    print('solute old')
                    for ncomp in range(self.numSoluteComp):
                        print(repr(np.array(cyl.getSolution(ncomp + 1))), ',')
                    #cyl.base.printParams()
                    didReset = True
                    
                if n_iter_solve == 1: # try with small init ddt
                    cyl.ddt = min( 1.e-3,cyl.ddt)
                    cyl.reset();
                    didReset = True
                    
                if n_iter_solve == 2:
                    cyl.setParameter("Newton.MaxSteps", "200")
                    cyl.setParameter("Newton.MaxTimeStepDivisions", "25")
                    if  errorWOnly or errorCOnly:
                        maxRelShift /= 10
                        cyl.setParameter("Newton.MaxRelativeShift", str(maxRelShift))# 
                    cyl.reset();
                    cyl.createNewtonSolver()
                    didReset = True
                elif n_iter_solve == 3:
                    print(rank,
                          'soil.solve() failed. making the computation more precise')
                    cyl.setParameter("Newton.EnableResidualCriterion", "true") # sometimes helps, sometimes makes things worse
                    cyl.setParameter("Newton.EnableAbsoluteResidualCriterion", "true")
                    cyl.setParameter("Newton.SatisfyResidualAndShiftCriterion", "true")
                    if errorWOnly or errorCOnly:
                        maxRelShift /= 10
                        cyl.setParameter("Newton.MaxRelativeShift", str(maxRelShift))# 
                    cyl.reset();
                    cyl.createNewtonSolver()
                    didReset = True
                elif (n_iter_solve < 8) and (( inner_fluxes_water_temp < 0 or outer_fluxes_water_temp < 0) or (min(outer_fluxes_solMucil_temp) <0. or min(inner_fluxes_solMucil_temp)<0.)): #
                    cyl.reset();didReset = True # do reset before re.-distributing the sink
                    inner_fluxes_water_temp, outer_fluxes_water_temp, inner_fluxes_solMucil_temp, outer_fluxes_solMucil_temp = self._adjust_fluxes(cyl,dt,  inner_fluxes_water_temp, outer_fluxes_water_temp, inner_fluxes_solMucil_temp, outer_fluxes_solMucil_temp,n_iter_solve, errorCOnly, errorWOnly)
                            
                        
                elif (maxRelShift < 1e-5) :#random upper limit
                    maxRelShift *= 10.
                    cyl.setParameter("Newton.MaxRelativeShift", str(maxRelShift))
                    cyl.reset();
                    cyl.createNewtonSolver()
                    didReset = True
                else:
                    print(rank,
                      'gId',gId,'ERROR, GIVING UP. errorCOnly',errorCOnly,errorWOnly,
                      'QflowIn',inner_fluxes_water_temp, 
                      ' Qflowout',outer_fluxes_water_temp)
                        
                    cyl.setParameter("Newton.MaxRelativeShift",
                                 str(cyl.MaxRelativeShift))
                    redoSolve = False
                    self.solve_gave_up = True
                    cyl.reset();
                    cyl.createNewtonSolver()
                    didReset = True
                    
                
                
                assert didReset
                assert self._check_Ccontent(cyl)
                assert self._check_Wcontent(cyl)
                    
                n_iter_solve += 1
            
            assert self._check_Ccontent(cyl)
            assert self._check_Wcontent(cyl)
            
        if verbose:
            print('finished solve 1ds gId', cyl.gId, repr(cyl.getWaterContent()), 
                    repr(cyl.getSolutionHead()),
                  'inner_fluxes_real',inner_fluxes_real, 
                  'inner_fluxes_water_temp', inner_fluxes_water_temp, 
                  repr(cyl.getCellVolumes()),
                  'changeW', sum((cyl.getWaterContent() - oldWaterContent )*
                   cyl.getCellVolumes()))
        return inner_fluxes_real, outer_fluxes_water_temp, outer_fluxes_solMucil_temp, cyl_sources        
     
    def _adjust_fluxes(self, cyl, dt, inner_fluxes_water_temp, outer_fluxes_water_temp, 
                       inner_fluxes_solMucil_temp, outer_fluxes_solMucil,  n_iter_solve, 
                       errorCOnly, errorWOnly):
        l = cyl.segLength
        gId = cyl.gId
        if (n_iter_solve < 8):                   
            divVale = 0.9
        else:
            print(rank,'gId',gId,'lId',lId,'ERROR, setting fluxes to 0. errorCOnly',errorCOnly,errorWOnly)
            divVale = 0.
        # water
        if (not errorCOnly) or (errorWOnly) or (min(outer_fluxes_solMucil) >0. and min(inner_fluxes_solMucil_temp)>0.):
            if (inner_fluxes_water_temp < 0 or outer_fluxes_water_temp < 0) and (inner_fluxes_water_temp <= outer_fluxes_water_temp):
                inner_fluxes_water_temp *=divVale
                qIn = inner_fluxes_water_temp/ (2 * np.pi * self.radii[gId] * l) # [cm3/day] -> [cm /day]
                cyl.setInnerFluxCyl(qIn) 
                
            elif (inner_fluxes_water_temp < 0 or outer_fluxes_water_temp < 0):
                outer_fluxes_water_temp *= divVale
                qOut = cyl.distributeSource(dt, outer_fluxes_water_temp,inner_fluxes_water_temp * dt,eqIdx = 0,plantM=self)
                outer_fluxes_water_temp = sum(qOut) # maybe outer_fluxes_water_temp was re-adjusted inside distributeSource
                    
        # solutes
        if (min(outer_fluxes_solMucil) <0. or min(inner_fluxes_solMucil_temp)<0.) and (min(inner_fluxes_solMucil_temp) <= min(outer_fluxes_solMucil)):
            inner_fluxes_solMucil_temp[np.where(inner_fluxes_solMucil_temp == min(inner_fluxes_solMucil_temp))] *=divVale
            valueBotBC = np.array(cyl.base.getbcBotValue())
            valueBotBC[len(inner_fluxes_solMucil_temp)] = inner_fluxes_solMucil_temp/ (2 * np.pi * self.radii[gId] * l) # [cm3/day] -> [cm /day]
            cyl.setSoluteBotBC(cyl.typeBC, valueBotBC)
            
        elif (min(outer_fluxes_solMucil) <0. or min(inner_fluxes_solMucil_temp)<0.):
            outer_fluxes_solMucil[np.where(outer_fluxes_solMucil == min(outer_fluxes_solMucil))] *= divVale # or id where getSolution < 0
            QflowOutCellLim = cyl.distributeSources(dt, outer_fluxes_solMucil, inner_fluxes_solMucil_temp * dt,
                                      np.array([nc+1 for nc in range(self.numDissolvedSoluteComp)]), plantM=self)
            for jj in range(self.numDissolvedSoluteComp):
                outer_fluxes_solMucil[jj] = sum(QflowOutCellLim[jj]) # maybe outer_fluxes_solMucil was re-adjusted inside distribSource
                                                   
        return inner_fluxes_water_temp, outer_fluxes_water_temp, inner_fluxes_solMucil_temp, outer_fluxes_solMucil
        


        
    def getTotCNContent(self, cyl):
        if isinstance(cyl, AirSegment):
            return np.full(self.numSoluteComp,0.)
        
        return cyl.getTotCNContent()
        
    
    
    def getDeff(self,Ds,phi,theta):
        return Ds* (theta**(10/3)) / (phi**2)
    
    def do1d1dFlow_(self, soilVals,#wat. pot. (cm)
                        seg_valuesW,# theta (cm3/cm3)
                        seg_valuesSol,# mol/cm3 wat
                        seg_valuesmucil,# mol/cm3 wat
                        verbose=False):
        raise Exception
        # manually set 1d-1d flow
        cellIds = self.getCellIds()
        self.flow1d1d_w = np.zeros(seg_valuesW.shape)
        self.flow1d1d_sol = np.zeros(seg_valuesW.shape)
        self.flow1d1d_mucil = np.zeros(seg_valuesW.shape)
        
            
        Ds = self.soilModel.Ds *(24*3600) *10000 # m^2/s to cm^2/d
        Dl = self.soilModel.Dl *(24*3600) *10000 # m^2/s to cm^2/d
        
        for cyl in self.cylsSoil:
            gId = cyl.gId
            cellid = self.seg2cell[gId]
            pHeadMean = soilVals[cellid]
            thetaMean = vg.water_content( soilVals[cellid], self.vg_soils[0])
            # hydraulic conductivity [cm/day]
            Kf = vg.hydraulic_conductivity(pHeadMean, self.vg_soils[0])# 
            
            # get all seg ids
            segIds = self.cell2seg[cellid]
            rootIds = np.array([sid for sid in segIds if ((self.organTypes[sid] == 2) and (sid != gId))])
            
            if len(rootIds) > 0:# at least 1 other segment
                
                #use position of node y
                pos0 =   np.array(self.nodesPos[gId+1])#cm
                posX =   np.array(self.nodesPos[rootIds+1])#cm
                
                dist =( (pos0-posX)**2).sum(axis=1)**(1/2) #cm
                
                pHead0 = vg.pressure_head(seg_valuesW[gId],self.soilModel.vg_soils[0])  #cm
                pHeadX = np.array([
                    vg.pressure_head(thetaCyl_,
                                    self.soilModel.vg_soils[0]) for thetaCyl_ in seg_valuesW[rootIds]]) #cm
                
                #exchange surface == 1/4 side surface of segs
                Surf = np.minimum(cyl.l*cyl.a_out,
                                  self.seg_length[rootIds]*self.outer_radii[rootIds] )*np.pi/2 
                # if flow1d1d_w > 0, gain for the cylinder
                flow1d1d_w = Kf * ( pHeadX - pHead0 )/dist *Surf #cm3/d
                self.flow1d1d_w[gId] = sum(flow1d1d_w)
                
                CC = np.where(flow1d1d_w > 0,seg_valuesSol[rootIds],seg_valuesSol[gId])
                
                self.flow1d1d_sol[gId] =sum(flow1d1d_w * CC)
                phi = self.soilModel.vg_soils[0].theta_S
                
                diffS = self.getDeff(Ds,phi,thetaMean)                                                           
                # if flow1d1d_sol > 0, gain for the cylinder                               
                self.flow1d1d_sol[gId] += sum(diffS/dist * (seg_valuesSol[rootIds] - seg_valuesSol[gId])*Surf)
                
                diffL=  self.getDeff(Dl,phi,thetaMean)
                # cm^2/d /cm * [mol/cm^3] / cm = mol/cm^3/d
                # if flow1d1d_mucil > 0, gain for the cylinder 
                self.flow1d1d_mucil[gId] = sum(
                    diffL/dist * (seg_valuesmucil[rootIds] -seg_valuesmucil[gId] )*Surf)
                                   
        self.flow1d1d_wG = comm.gather(self.flow1d1d_w,root=0)  
        self.flow1d1d_solG = comm.gather(self.flow1d1d_sol,root=0)  
        self.flow1d1d_mucilG = comm.gather(self.flow1d1d_mucil,root=0) 
        
        gotError = 0
        
        if rank == 0:
            self.flow1d1d_w = np.array(self.flow1d1d_wG).sum(axis=0)  
            self.flow1d1d_sol = np.array(self.flow1d1d_solG).sum(axis=0)  
            self.flow1d1d_mucil = np.array(self.flow1d1d_mucilG).sum(axis=0) 
            
        
            for cc in cellIds:
                segs = self.cell2seg[cc]
                divw = 1 if np.max(abs(self.flow1d1d_w[segs])) == 0 else np.max(self.flow1d1d_w[segs])
                divsol = 1 if np.max(abs(self.flow1d1d_sol[segs])) == 0 else np.max(self.flow1d1d_sol[segs])
                divmucil = 1 if np.max(abs(self.flow1d1d_mucil[segs]) )== 0 else np.max(self.flow1d1d_mucil[segs])

                vals = np.array([sum(self.flow1d1d_w[segs]),
                             sum(self.flow1d1d_sol[segs]),
                             sum(self.flow1d1d_mucil[segs])])
                divs = np.array([divw,divsol,divmucil])
                
                gotError =  max(abs(vals/divs)) > 0.1
                
        gotError = comm.bcast(gotError,root=0)            
        assert not gotError
    
    def getRealisedInnerFluxes(self, compId)    :
        if compId == 0:
            return self.getXcyl(data2share=self.seg_fluxes_limited, doSum = False, reOrder = True) 
        else:
            return self.getXcyl(data2share=self.seg_fluxes_limited_CN_Out[compId - 1], doSum = False, reOrder = True) 
                
    
    #def splitSoilVals(self, soilVals,seg_values, seg_volume, dt, verbose=False, isWater=False):
    def splitSoilVals(self, soilVals,compId, dt): 
        """ 
        Split soilFlux array according to the values in seg_values.
        
        Parameters:
        - soilVals: array-like, fluxes which can be <, =, or > 0
        - seg_values: array-like, concentration
        - verbose: bool, whether to print detailed debug information
        - isWater: bool, whether the operation is related to water
        """
        isWater = compId == 0
        isDissolved =  (compId <= self.numDissolvedSoluteComp)
        
        if isWater:
            seg_content = self.getWaterVolumesCyl(doSum = False, reOrder = True) #cm3 water per 1d model
        else:
            seg_content = self.getContentCyl(idComp= 1 , doSum = False, reOrder = True) 
        
        if isDissolved:
            seg_volume = self.getWaterVolumesCyl(doSum = False, reOrder = True) #cm3 water per 1d model
            seg_volume[self.airSegs] = 1. # to avoid division by 0.
        else:
            seg_volume = self.getVolumesCyl(doSum = False)
            
        seg_values = seg_content/seg_volume
        if rank ==0:
            
            cellIds = self.getCellIds()    
            organTypes = np.array(self.organTypes)

            if self.debugMode and (organTypes != 2).any():
                # assert min(seg_values) >= 0. # as we have initial concentration + flow, could have value < 0
                assert seg_values.shape == (len(organTypes), )
                assert (seg_values[np.where(organTypes != 2)] == 0.).all()
                
            if (soilVals[cellIds] != 0.).any():
                splitVals = np.array(self.splitSoilVals_(
                    soilVals, cellIds, isWater, seg_values, seg_volume, 
                    dt, self.vg_soils[0].theta_S, 
                      self.theta_wilting_point))
            else:
                splitVals = np.zeros(len(self.organTypes))

            return splitVals
        
    
    def _map(self, x, dtype = np.float64, doMap = True):
        """Converts @param x to a numpy array and maps it to the right indices                 """
        
        indices = self.eidx_all_ # gather segment indices from all threads
        
        if rank == 0:  # only for rank 0 it is not empty
            if doMap:
                if (len(indices) != len(x)):
                    print("len(indices) != len(x)",'indices',indices,'x',x)
                    print(rank, f"RhizoMappedSegments._map: indices and values have different length:{len(indices)} vs {len(x)}")
                    raise Exception
                
                p = np.zeros((len(x),), dtype = dtype)
                for i in range(0, len(indices)):  #
                    p[indices[i]] = x[i]  
            else:
                p = x
            return p
        else:
            return np.array([], dtype = dtype)
    
    def gather(self,x, dtype = np.float64, doMap = True):
        return self._map(self._flat0(comm.gather(x, root=0)),dtype, doMap)

