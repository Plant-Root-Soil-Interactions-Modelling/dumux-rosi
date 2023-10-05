# import sys; sys.path.append("../modules/"); sys.path.append("../modules/fv/"); sys.path.append("../../../CPlantBox/"); sys.path.append("../../../CPlantBox/src");
# sys.path.append("../../build-cmake/cpp/python_binding/")

import plantbox as pb
import functional.xylem_flux as xylem_flux
import sys
from functional.xylem_flux import XylemFluxPython
from rosi_richards_cyl import RichardsCylFoam  # C++ part (Dumux binding) of cylindrcial model
from rosi_richards2c_cyl import Richards2CCylFoam  # C++ part (Dumux binding)
from rosi_richardsnc_cyl import RichardsNCCylFoam  # C++ part (Dumux binding)
from rosi_richards10c_cyl import Richards10CCylFoam # C++ part (Dumux binding)
from richards_no_mpi import RichardsNoMPIWrapper  # Python part of cylindrcial model (a single cylindrical model is not allowed to run in parallel)
from fv.fv_grid import *
import fv.fv_richards as rich  # local pure Python cylindrical models
import functional.van_genuchten as vg
from scenario_setup import write_file_array

import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import norm
from scipy.interpolate import griddata as gd
from mpl_toolkits.mplot3d import Axes3D
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank(); max_rank = comm.Get_size()
import multiprocessing
from multiprocessing import Process, active_children
import psutil
from threading import Thread 
from air_modelsPlant import AirSegment
from scipy.interpolate import PchipInterpolator,  CubicSpline

class RhizoMappedSegments(pb.MappedPlant):#XylemFluxPython):#
    """
        Adds 1-dimensional rhizospere models to each root segment of a MappedSegments (or later MappedPlant)    
        
        modes:        
        "dumux"              inner boundary is FluxCyl,                    <- currently best
        "dumux_exact"        inner boundary is rootSystemExact (same as "dumux", but slower, for experiments...) 
        "python"             inner boundary is flux_in,
        "python_exact"       inner boundary is rootsystem_exact (slower, for experiments...)
    """

    # TODO copy mapped segments constructors (!)...

    def __init__(self, wilting_point, NC, logbase, mode, soil, recreateComsol_, 
                usemoles, seedNum=None):
        """ @param file_name is either a pb.MappedRootSystem, pb.MappedSegments, or a string containing a rsml filename"""
        if not seedNum is None:
            super().__init__(seednum = seedNum)
        else:
            super().__init__()
        self.NC = NC #dof +1
        self.logbase = logbase
        self.mode = mode  # more precise RichardsCylFoam, mode="dumux"
        #print(self.mode)

        # changes with growing plant (to update)
        #self.max_rank = max_rank
        self.nsOld = 0
        self.toAdd = np.array([])
        self.repartitionOld = np.array([0 for i in range( max_rank)])
        
        
        self.cyls = []
        self.cylsSoil = []
        self.outer_radii = None
        self.rhizoVol = None
        self.seg_length = np.array([]) 
        self.eidx = np.array([]) #root segments id
        self.cylidx = [] # perirhizal cylinders id
        self.cylSoilidx = [] # perirhizal cylinders id without airSegments
        self.dx2 = []
        self.soilModel = soil
        self.recreateComsol = recreateComsol_
        self.useMoles = usemoles
        
    
        self.wilting_point = self.soilModel.wilting_point
        self.molarMassWat = self.soilModel.molarMassWat # [g/mol]
        self.densityWat_m3 = self.soilModel.densityWat_m3 #[g/m3]
        # [mol/m3] = [g/m3] /  [g/mol] 
        self.molarDensityWat_m3 =  self.densityWat_m3 /self.molarMassWat # [mol/m3]    
        self.numFluidComp = self.soilModel.numFluidComp
        self.numComp = self.soilModel.numComp
        self.soil = self.soilModel.soil
        self.vg_soil = self.soilModel.vg_soil
        
        self.solidDensity = self.soilModel.solidDensity#2700 # [kg/m^3 solid]
        self.solidMolarMass = self.soilModel.solidMolarMass#60.08e-3 # [kg/mol] 
            
        # [mol / m3 solid] =[kg/m^3 solid] / [kg/mol] 
        self.solidMolDensity = self.solidDensity/self.solidMolarMass
        # [mol / m3 scv] = [mol / m3 solid] * [m3 solid /m3 space]
        self.bulkDensity_m3 = self.solidMolDensity*(1.- self.soilModel.vg_soil.theta_S) #porosity == theta_s
        
        # not sorption atm
        self.canBeNull = np.where(self.soilModel.ICcc == 0.)[0] + 1
        
        self.setSoilData()# initial soil condition
            
        self.dx2 = []
        self.eidx = np.array([], dtype=np.int64)
        # additional variables
        self.last_dt = 0.

        
        
    def initializeRhizo(self):#, cc = None):
        """ calls the specific initializer according to @param mode (no parallelisation)  
        @param soil     van genuchten parameters as list        
        @param x        is the solution (or initial condition) of the soil model
        @þaram eidx     in case of mpi the cylinder indices per process
        """
        pass
        
        
        # self.update(x, eidx, cc)
        
    def broadcastPlantShape(self):
        #backups
        self.outer_radii_old = self.outer_radii
        self.seg_length_old = self.seg_length
        self.rhizoVol_old = self.rhizoVol
        
        print('passing shape', rank)
        if rank == 0:
    
            self.outer_radii = np.array(self.segOuterRadii(type = 0))  # in the future, segment radius might vary with time. TODO: how to handle that?
            airSegs = np.array(list(set(np.concatenate((self.cell2seg[-1],np.where(np.array(self.organTypes) != 2)[0])) )))#aboveground
            self.outer_radii[airSegs] = np.array(self.radii)[airSegs]*1.1 #dummy outer radii for segments in the aire, also do it for shoot seg in soil
            if isinstance(self.outer_radii_old, type(self.outer_radii)):
                try:
                    assert (self.outer_radii_old >= self.outer_radii[:len(self.outer_radii_old)]).all()
                except:
                    print("diffradii", self.outer_radii_old - self.outer_radii[:len(self.outer_radii_old)])
                    print('old',self.outer_radii_old ,'radii', self.outer_radii, 'lenold',len(self.outer_radii_old))
                    segsWithIssues = np.where(self.outer_radii_old >= self.outer_radii[len(self.outer_radii_old)])[0]
                    print(segsWithIssues)
                    print([self.seg2cell[ss] for ss in segsWithIssues])
                    raise Exception
                    
            self.seg_length = np.array(self.segLength())#node might have shifted: new length for pre-existing segments          

            ns = len(self.seg_length)
            dcyl = int(np.floor(ns /  max_rank))
            repartition = np.array([dcyl for i in range( max_rank)])
            repartition[ max_rank -1] += ns - max_rank * dcyl #residual
            self.toAdd = repartition - self.repartitionOld
            self.repartitionOld = repartition
            self.organTypes =np.array( self.organTypes)
            assert sum(repartition) == ns
            assert sum(self.toAdd) == (ns - self.nsOld )
            print()
            #for next time step
        
        
        self.toAdd = comm.bcast(self.toAdd, root = 0) 
        self.seg_length = comm.bcast(self.seg_length, root = 0)  
        self.outer_radii = comm.bcast(self.outer_radii, root = 0)
        self.radii =np.array( comm.bcast(self.radii, root = 0)    )
        self.seg2cell = comm.bcast(self.seg2cell, root = 0)    
        self.cell2seg = comm.bcast(self.cell2seg, root = 0)         
        self.organTypes =np.array( comm.bcast(self.organTypes, root = 0) ) 
        self.organTypesArray =np.array( comm.bcast(self.organTypes, root = 0)   )
        self.organTypesArray =np.array( self.organTypesArray)
        self.rhizoVol = (self.outer_radii**2 - np.array(self.radii)**2)*np.pi*self.seg_length
        print('finished plant shape', rank)
        if rank > 0:
            print('passed values', self.toAdd, self.nsOld, self.cell2seg)
        
        
    def update(self):
        """ calls the specific initializer according to @param mode (no parallelisation)  
        @param soil     van genuchten parameters as list        
        @param x        is the solution (or initial condition) of the soil model
        @þaram newEidx  in case of mpi the cylinder indices per process
        """
        
        self.broadcastPlantShape() 
        newEidx = np.array([i for i in range(self.toAdd[rank])]) + sum(self.toAdd[:rank]) +self.nsOld
        if newEidx is None:
            newEidx = np.array(range(len(self.eidx), len(self.radii)), np.int64)  # segment indices for process
        self.eidx = np.concatenate((self.eidx,np.array(newEidx, dtype = np.int64)), dtype = np.int64)
        
        self.seg_length_ = self.seg_length[self.eidx]
        self.radii_ = np.array(self.radii)[self.eidx]
        self.outer_radii_ = self.outer_radii[self.eidx]
        self.rhizoVol_ = self.rhizoVol[self.eidx]
        
        cellIds = np.fromiter(self.cell2seg.keys(), dtype=int)
        cellIds =  np.array([x for x in cellIds if x >= 0])#take out air segments
        
        sizeSoilCell = comm.bcast(self.soilModel.getCellVolumes()/1e6, root = 0) #m3
        try:
            maxlim1d3d = max(self.maxDiff1d3dCW_rel[0]*10, 1e-13)
        except:
            maxlim1d3d = 1e-13
        self.checkMassOMoleBalance2( sourceWat = np.full(len(sizeSoilCell),0.), # cm3/day 
                                     sourceSol = np.full((self.numComp, len(sizeSoilCell)),0.), # mol/day
                                     dt = 0.,        # day    
                                     seg_fluxes = 0.,# [cm3/day]
                                     doWater = True, doSolute = True, doSolid = True,
                                     useSoilData = True,
                                     diff1d3dCW_rel_lim = maxlim1d3d)
        print("do vol balance", rank)
        self.checkVolumeBalance(finishedUpdate = False)
        print("checkRadii, rank", rank)
        self.checkRadii()
        
        print("updateOld", rank)
        for i, cyl in enumerate(self.cyls):
            self.updateOld(i, cyl)
        
        volLeftOver = np.array([self.get_Vol_leftoverI(i) for i in range(len(sizeSoilCell))])# m3
        watVol_leftover = np.array([self.get_watVol_leftoverI(i) for i in range(len(sizeSoilCell))])#m3
        try:
            assert watVol_leftover.shape == (len(sizeSoilCell),)
            assert isinstance(watVol_leftover[cellIds[0]], float)
        except:
            print(watVol_leftover)
            raise Exception
        
        c_content_leftover = dict([(i,np.array([self.getC_content_leftoverI(i, ncomp) for ncomp in range(1, self.numComp + 1)])) for i in cellIds])# mol    
            
        #m3 wat * (mol wat/m3 wat) or m3 scv * (mol scv/m3 scv)
        phaseMol = dict([(i, np.array([watVol_leftover[i]* self.molarDensityWat_m3 if ncomp <= self.numFluidComp else volLeftOver[i]*self.bulkDensity_m3 for ncomp in range(1, self.numComp + 1)])) for i in cellIds])
        
        # in mol/mol wat or mol/mol bulk soil
        CC_leftover = self.get_CC_leftover(c_content_leftover, phaseMol, cellIds)
        
        try:
            assert CC_leftover[cellIds[0]].shape == (self.numComp,)
        except:
            print(c_content_leftover,  CC_leftover)
            raise Exception
        
        
        
        XX_leftover = self.get_XX_leftover(watVol_leftover[cellIds], volLeftOver[cellIds], cellIds )# dict(zip(cellIds, XX_leftover))
        
        print('start init',rank, newEidx)
        for gId in newEidx:#only initialize the new eidx
            print('do init',rank, gId)
            self.initialize_(gId,XX_leftover,CC_leftover)
            print('done init',rank, gId)
        print('finished init',rank)
        self.cylidx = self.eidx
        self.cylSoilidx = np.array([gId for lId, gId in enumerate(self.cylidx) if (not isinstance(self.cyls[lId], AirSegment))])
        self.cylsSoil = np.array([cyl for cyl in self.cyls if (not isinstance(cyl, AirSegment))])
        
        self.checkVolumeBalance(finishedUpdate = True)
        print("afterUpdate")
        self.checkMassOMoleBalance2( sourceWat = np.full(len(sizeSoilCell),0.), # cm3/day 
                                     sourceSol = np.full((self.numComp, len(sizeSoilCell)),0.), # mol/day
                                     dt = 0.,        # day    
                                     seg_fluxes = 0.,# [cm3/day]
                                     doWater = True, doSolute = True, doSolid = True,
                                     useSoilData = True,
                                     diff1d3dCW_rel_lim = max(self.maxDiff1d3dCW_rel[0]*2, 1e-12) )# error after update should be = error before update
        # maybe check here that rhizo concentration + outerBC_flow == 3D concentration
        # self.checkMassOMoleBalance()#doSolute = False)
        self.nsOld = sum(self.repartitionOld)
        self.nsOld = comm.bcast(self.nsOld, root = 0)  
    
    def get_CC_leftover(self,c_content_leftover, #mol
                        phaseMol,     #mol water or mol scv
                        cellIds):   
        
        assert phaseMol[cellIds[0]].shape == (self.numComp,)
        assert c_content_leftover[cellIds[0]].shape == (self.numComp,)
        CC_leftover = np.full((len(cellIds), self.numComp),np.nan)
        for i, cid in enumerate(cellIds):
            pm = phaseMol[cid]
            c_content = c_content_leftover[cid]
            CC_leftover[i][np.where(pm != 0)] = c_content[np.where(pm != 0)]/pm[np.where(pm != 0)]
            # print(CC_leftover[i][np.where(pm != 0)] , c_content[np.where(pm != 0)],pm[np.where(pm != 0)])
            
        
        assert CC_leftover.shape == (len(cellIds), self.numComp)
        
        CC_leftover = dict([(cellIds[i],CC_leftover[i]) for i in range(len(cellIds)) ])

        return CC_leftover #mol/mol
        
    def get_XX_leftover(self,watVol_leftover, #m3
                        volLeftOver,     #m3
                        cellIds):   
        theta_leftOver = np.array([watVol_leftover[idvol]/vol if vol > 0 else np.nan for idvol, vol in enumerate(volLeftOver) ])
        lowTheta = np.where(theta_leftOver <self.vg_soil.theta_R )
        if len(volLeftOver[lowTheta]) > 0:
            assert max(volLeftOver[lowTheta]) < 1e-13 # low theta, no new segments in the soil voxel
        highTheta = np.where(theta_leftOver >self.vg_soil.theta_S )
        if len(volLeftOver[highTheta]) > 0:
            assert max(volLeftOver[highTheta]) < 1e-13 # low theta, no new segments in the soil voxel
        correctTheta = np.where(((theta_leftOver <= self.vg_soil.theta_S) & (theta_leftOver >= self.vg_soil.theta_R)))
        XX_leftover = np.full(theta_leftOver.shape, np.nan)
        if len(theta_leftOver[correctTheta]) > 0:
            try:
                XX_leftover[correctTheta] = np.array([vg.pressure_head( tlo, self.vg_soil) for tlo in theta_leftOver[correctTheta]])
            except:
                print(theta_leftOver[lowTheta],volLeftOver[lowTheta])
                print(theta_leftOver[highTheta],volLeftOver[highTheta])
                print(theta_leftOver[correctTheta])
                raise Exception
        
        # using dict() and zip() to convert arrays to dictionary
        return dict(zip(cellIds, XX_leftover))
        
    
    def phaseDensity(self, compId):# mol / m3
        return self.soilModel.phaseDensity( isDissolved = (compId <= self.numFluidComp))
    
    def checkRadii(self):
        """ vol soil == vol perirhizal (without the root volume)"""
        comm.Barrier()
        cellIds = np.fromiter(self.cell2seg.keys(), dtype=int)
        cellIds =  np.array([x for x in cellIds if x >= 0])#take out air segments
        for cellId in cellIds:
            #print("cellid",cellId, rank)
            vol_total = self.soilModel.getCellVolumes()
            if rank == 0:
                vol_total = vol_total[cellId] # cm3
                
            vol_total = comm.bcast(vol_total, root = 0) 
            
            idSegs = np.array(self.cell2seg[cellId])#all segments in the cell
            idCylsAll = np.array([ids for ids in idSegs if self.organTypes[ids]==2 ]) #all the segments which already have cylinder
            #print('rank',rank,idSegs ,idCylsAll)
            lengths_I = np.array(self.seg_length)[idCylsAll]
            radii_in = np.array(self.radii)[idCylsAll]
            radii_out = np.array(self.outer_radii)[idCylsAll]
            vol_rootRhizo = radii_out * radii_out * np.pi * lengths_I
            vol_root = radii_in * radii_in * np.pi * lengths_I 
            try:
                assert abs(((vol_total - sum(vol_rootRhizo - vol_root))/vol_total)*100) < 1e-12
            except:
                for rr in range(max_rank):
                    if rank == rr:
                        comm.Barrier()
                        print("\n",rank,cellId,"checkRadii error",  vol_total, 'vol_rootRhizo',vol_rootRhizo,sum(vol_rootRhizo))
                        print("\n",rank,cellId,'voltot-volrootrhizo',vol_total - sum(vol_rootRhizo))
                        print("\n",rank,cellId,'vol_root',vol_root, sum(vol_root))
                        print("\n",rank,cellId,'%error', ((vol_total - sum(vol_rootRhizo - vol_root))/vol_total)*100)
                        print("\n",rank,cellId,'radi in and out, lengths_I', radii_in,radii_out, lengths_I)
                        print("\n",rank,cellId, 'ots',self.organTypesArray,idSegs,type(self.organTypesArray), type(idSegs))
                        print("\n",rank,cellId, 'ots',self.organTypesArray[idSegs])
                        comm.Barrier()
                raise Exception
        comm.Barrier()
    
    def checkVolumeBalance(self, finishedUpdate):
        comm.Barrier()
        cellIds = np.fromiter(self.cell2seg.keys(), dtype=int)
        cellIds =  np.array([x for x in cellIds if x >= 0])#take out air segments
        for cellId in cellIds:
            vol_total = self.soilModel.getCellVolumes()
            if rank == 0:
                vol_total = vol_total[cellId] # cm3.
            vol_total =     comm.bcast(vol_total, root=0)
            idCylsAll, idCyls =  self.getIdCyllMPI(cellId)
				
            if(len(idCylsAll) >0):
                
                # localIdCyls = np.array([ np.where(self.eidx == i)[0] for i in idCyls4Vol]).flatten()
                # vol_rhizo_ = np.array([0])
                # try:
                    # if len(localIdCyls)>0:
                        # vol_rhizo_ = np.array([sum(self.cyls[i].getCellSurfacesCyl()) for i in localIdCyls ])*np.array(self.seg_length_)[localIdCyls] #cm3
                # except:
                    # print('error vol_rhizo_', rank, localIdCyls, idCyls4Vol, idCylsAll)
                    # raise Exception
                # vol_rhizo_ = comm.bcast(self._flat0(comm.gather(vol_rhizo_, root=0)), root=0)
                vol_rhizo = self.getVolumesCyl(idCyls)#cm3 sum(vol_rhizo_)
                    
                
                if (abs((vol_total - vol_rhizo)) > 1e-10):# or ((vol_total - vol_rhizo) < 0.):
                    print("checkVolumeBalance--error")
                    print(rank, cellId, vol_total, vol_rhizo, (vol_total - vol_rhizo),vol_rhizo_  )
                    localIdCyls = np.array([ np.where(self.eidx == i)[0] for i in idCyls4Vol]).flatten()
                    idSegs = np.array(self.cell2seg[cellId])#all segments in the cell
                    idCylsAllrr = np.array([ids for ids in idSegs if self.organTypes[ids]==2 ]) #all the segments which already have cylinder
                    lengths_I = np.array(self.seg_length)[idCylsAllrr]
                    radii_in = np.array(self.radii)[idCylsAllrr]
                    radii_out = np.array(self.outer_radii)[idCylsAllrr]
                    vol_rootRhizo = radii_out * radii_out * np.pi * lengths_I
                    vol_root = radii_in * radii_in * np.pi * lengths_I 
                    if rank == 0:
                        print(self.soilModel.getCellVolumes()[cellId])
                        print('radii vs volA',idCylsAllrr ,idCylsAll,idCyls4Vol,self.eidx[localIdCyls],localIdCyls,len(idCylsAllrr) ,len(idCylsAll))
                        print('radii vs volB', vol_rootRhizo-vol_root,  sum(vol_rootRhizo-vol_root))
                        print('radii vs volC',np.array(self.radii)[idCyls4Vol], np.array(self.outer_radii)[idCyls4Vol], np.array(self.seg_length)[idCyls4Vol])
                        print('radii vs volD',np.array(self.seg_length_)[localIdCyls] , localIdCyls)
                        print(np.array([cc.getCellSurfacesCyl() for cc in np.array(self.cyls)[localIdCyls] ]))
                        print(np.array([sum(cc.getCellSurfacesCyl()) for cc in np.array(self.cyls)[localIdCyls] ]))
                        # print(np.array([np.pi*(cc.base.rOut*cc.base.rOut-cc.base.rIn*cc.base.rIn)*10000 for cc in self.cyls ])*self.seg_length)
                        # print(sum(np.array([np.pi*(cc.base.rOut*cc.base.rOut-cc.base.rIn*cc.base.rIn)*10000 for cc in self.cyls ])*self.seg_length))
                        print("vol cyls",np.array([sum(self.cyls[i].getCellSurfacesCyl()) for i in localIdCyls])*np.array(self.seg_length_)[localIdCyls])
                        print("points", np.array([self.cyls[i].getPoints() for i in localIdCyls]))
                    raise Exception
        comm.Barrier()
    
                    
    def checkMassOMoleBalance2(self, sourceWat, # cm3/day 
                                     sourceSol, # mol/day
                                     dt,        # day    
                                     seg_fluxes,# [cm3/day]
                                     doWater = True, doSolute = True, doSolid = True,
                                     useSoilData = False,
                                     doKonz = True,
                                     diff1d3dCW_rel_lim = 1e-11):#would need to do it for each cell, not overall
        print("checkMassOMoleBalance2")
        
        cellIds = np.fromiter(self.cell2seg.keys(), dtype=int)
        cellIds =  np.array([x for x in cellIds if x >= 0])#take out air segments
        print("cellIds",rank, cellIds)
        
        self.sumDiff1d3dCW_abs = np.full(self.numComp+1, 0.)
        self.sumDiff1d3dCW_rel = np.full(self.numComp+1, 0.)
        self.maxDiff1d3dCW_abs = np.full(self.numComp+1, 0.)
        self.maxDiff1d3dCW_rel = np.full(self.numComp+1, 0.)
        if doSolid:
            NC = self.numComp
        elif doSolute:
            NC = self.numFluidComp
        else:
            NC = 0
        
        comm.Barrier()
        for cellId in cellIds:
            comm.Barrier()
            idCylsAll, idCyls = self.getIdCyllMPI(cellId)
            print("rank",rank,"cell",cellId,  "\tidCylsAll",idCylsAll)
            if len(idCylsAll)> 0:
                if useSoilData:
                    wat_total = self.soilWatVol_old * 1e6# [cm3].
                else:
                    wat_total = self.soilWatVol_old * 1e6# [cm3].
                if rank == 0:
                    wat_total = wat_total[cellId] # [cm3].
                wat_total = comm.bcast(wat_total, root=0)
                localIdCyls = np.array([ np.where(self.eidx == i)[0] for i in idCyls]).flatten()#np.where(idCyls in self.eidx)[0]
                #wat_rhizo_ = np.array([sum(self.cyls[i].getWaterVolumesCyl(self.seg_length_[i])) for i in localIdCyls ]) #cm3
                wat_rhizo = self.getWaterVolumesCyl(idCyls)# sum(comm.bcast(self._flat0(comm.gather(wat_rhizo_, root=0)), root=0))
                    
                print("in water",rank, wat_total,wat_rhizo, idCyls , self.eidx,idCyls in self.eidx)
                
                if rank == 0:
                    try:
                        testwat = sourceWat[cellId]
                    except:
                        print(sourceWat,  cellId)
                        raise Exception
                        
                    #print(cellId, 0,wat_total , sourceWat[cellId] * dt , wat_rhizo, wat_rhizo_)
                    diff1d3dCW_abs = abs(wat_total+ sourceWat[cellId] * dt - wat_rhizo)
                    diff1d3dCW_rel = abs(diff1d3dCW_abs/(wat_total+ sourceWat[cellId] * dt) *100)
                    self.sumDiff1d3dCW_abs[0] += diff1d3dCW_abs
                    self.sumDiff1d3dCW_rel[0] += diff1d3dCW_rel
                    self.maxDiff1d3dCW_abs[0] = max(self.maxDiff1d3dCW_abs[0], diff1d3dCW_abs)
                    self.maxDiff1d3dCW_rel[0] = max(self.maxDiff1d3dCW_rel[0], diff1d3dCW_rel)
                    if diff1d3dCW_rel > diff1d3dCW_rel_lim:#1e-10 :
                        print("checkMassOMoleBalance2error")
                        print(cellId, 0, wat_total, wat_rhizo, sourceWat[cellId], dt)
                        print(0, wat_total, wat_rhizo, useSoilData)
                        #print(self.soilModel.getWaterVolumes()[cellId])
                        print(self.soilWatVol_old[cellId] * 1e6)
                        print(np.array([sum(self.cyls[i].getWaterVolumesCyl(self.seg_length[i])) for i in idCyls ]))
                        print("Diff1d3dW",diff1d3dCW_abs, diff1d3dCW_rel, self.sumDiff1d3dCW_abs[0],
                        self.sumDiff1d3dCW_rel[0], self.maxDiff1d3dCW_abs[0], self.maxDiff1d3dCW_rel[0],
                        diff1d3dCW_rel_lim)

                        raise Exception          
                    
                for idComp in range(1, NC +1): 
                    isDissolved = (idComp <= self.numFluidComp)
                    print("check comp", rank, idComp, len(idCylsAll))
                    if useSoilData:
                        mol_total = self.soilContent_old[idComp - 1] # solute content [mol]
                    else:
                        mol_total = self.soilModel.getContent(idComp, isDissolved) # solute content [mol].
                    if rank == 0:
                        mol_total = mol_total[cellId] # solute content [mol]
                    
                    mol_total = comm.bcast(mol_total, root=0)
                    
                    if idComp not in self.canBeNull:
                        assert mol_total > 0.
                    else:
                        assert mol_total >= 0.
                                        
                    #mol_rhizo_ = np.array([sum(self.cyls[i].getContentCyl( idComp, isDissolved, self.seg_length_[i] )) for i in localIdCyls ]) #cm3
                    mol_rhizo = self.getContentCyl(idCyls, idComp)#sum(comm.bcast(self._flat0(comm.gather(mol_rhizo_, root=0)), root=0))
                        
                            
                    if rank == 0:
                        try:
                            sourceSol_ = sourceSol[idComp-1][cellId]             
                        except:
                            print(sourceSol, idComp, cellId)
                            raise Exception
                        if (mol_total  + sourceSol_ * dt  == 0) :
                            diff1d3dCW_abs = abs(mol_total  + sourceSol_ * dt)
                            diff1d3dCW_rel = np.nan
                            #print(cellId, idComp,mol_total , sourceSol[idComp-1][cellId] * dt , mol_rhizo, mol_rhizo_)
                            if ((mol_rhizo > 0)):
                                print(idComp, mol_total, mol_rhizo)
                                raise Exception
                        else:
                            #print(cellId, idComp,mol_total , sourceSol[idComp-1][cellId] * dt , mol_rhizo, mol_rhizo_)
                            #print(cellId, idComp,mol_total , sourceSol[idComp-1][cellId] * dt)
                            diff1d3dCW_abs = abs(mol_total + sourceSol_ * dt - mol_rhizo)
                            diff1d3dCW_rel = abs(self.sumDiff1d3dCW_abs[idComp]/(mol_total +sourceSol_ * dt) *100)
                            #print("content", cellId, idComp, mol_total, mol_rhizo,"sourceSol_ * dt",sourceSol_ * dt, 
                            #    abs((mol_total + sourceSol_ * dt - mol_rhizo)/(mol_total + sourceSol_ * dt)*100))
                            if self.sumDiff1d3dCW_rel[idComp] > 0.1:
                                print(cellId, idComp,mol_total , sourceSol_ , dt)
                                print(mol_rhizo, abs((mol_total + sourceSol_ * dt - mol_rhizo)/(mol_total + sourceSol_* dt)*100))
                                raise Exception
                                
                        
                        self.sumDiff1d3dCW_abs[idComp] += diff1d3dCW_abs
                        self.sumDiff1d3dCW_rel[idComp] += diff1d3dCW_rel
                        self.maxDiff1d3dCW_abs[idComp] = max(self.maxDiff1d3dCW_abs[idComp], diff1d3dCW_abs)
                        self.maxDiff1d3dCW_rel[idComp] = max(self.maxDiff1d3dCW_rel[idComp], diff1d3dCW_rel)
                        print("check compBis", rank, mol_rhizo,mol_total, sourceSol_ * dt,"error",diff1d3dCW_abs, diff1d3dCW_rel)
            comm.Barrier()
        print("finished checkmass2", rank)
        comm.Barrier()        
    
    def getIdCyllMPI(self,cellId, getidCylsAll = True):    
        idSegs = self.cell2seg[cellId]#all segments in the cell
        idCyls =np.array( list(set(idSegs).intersection(self.cylSoilidx))) #all the segments which already have cylinder
        if getidCylsAll:
            return np.array( comm.bcast(self._flat0(comm.gather(idCyls, root=0)), root=0)), idCyls
        else:
            return idCyls
            
    def getXcyl(self,X_rhizo, doSum = True, reOrder = False):
        X_rhizo = np.array( comm.bcast(self._flat0(comm.gather(X_rhizo, root=0)), root=0))
        if reOrder:
            idCyll = np.array( comm.bcast(self._flat0(comm.gather(self.eidx, root=0)), root=0))
            X_rhizo[idCyll] = X_rhizo
        if doSum:
            X_rhizo = sum(X_rhizo)
        return X_rhizo
    def getLocalIdCyls(self,idCyls=None):
        if idCyls is None:
            idCyls = self.eidx
        return np.array([ np.where(self.eidx == i)[0] for i in idCyls]).flatten()#np.where(idCyls in self.eidx)[0] 
        
    def getContentCyl(self,idCyls=None, idComp=1, doSum = True, reOrder = False):#mol
        isDissolved = (idComp <= self.numFluidComp)
        localIdCyls =   self.getLocalIdCyls(idCyls)             
        mol_rhizo = np.array([sum(self.cyls[i].getContentCyl( idComp, isDissolved, self.seg_length_[i] )) for i in localIdCyls ]) #cm3
        return self.getXcyl(mol_rhizo, doSum, reOrder)
        
    def getTotCContentAll(self,idCyls=None, doSum = True, reOrder = False):#mol
        localIdCyls =   self.getLocalIdCyls(idCyls)                      
        mol_rhizo = np.array([sum(self.getTotCContent(self.cyls[i], self.seg_length_[i] )) for i in localIdCyls ]) #mol
        return self.getXcyl(mol_rhizo, doSum, reOrder)
        
    def getWaterVolumesCyl(self,idCyls=None, doSum = True, reOrder = False):#cm3
        localIdCyls =   self.getLocalIdCyls(idCyls)                      
        wat_rhizo = np.array([sum(self.cyls[i].getWaterVolumesCyl(self.seg_length_[i])) for i in localIdCyls ]) #cm3
        return self.getXcyl(wat_rhizo, doSum, reOrder)
        
    def getVolumesCyl(self,idCyls=None, doSum = True, reOrder = False):#cm3
        localIdCyls =   self.getLocalIdCyls(idCyls)                     
        V_rhizo = np.array([sum(self.cyls[i].getCellSurfacesCyl()*self.seg_length_[i]) for i in localIdCyls ]) #cm3
        return self.getXcyl(V_rhizo, doSum, reOrder)
        
    def getC_rhizo(self,soilShape, idComp = 1, konz = True): # if konz:  mol/m3 wat or mol/m3 scv, else: mol
        """ return component concentration or content, only in voxel with 1D-domains"""
        comm.Barrier()
        isDissolved = (idComp <= self.numFluidComp)
        
        contentOrKonz = np.full(soilShape, 0.)
        
            
        cellIds = np.fromiter(self.cell2seg.keys(), dtype=int)
        cellIds =  np.array([x for x in cellIds if x >= 0])#take out air segments
        for idCell in cellIds:
            idSegs = self.cell2seg[idCell]#all segments in the cell
            idCylsAll, idCyls= self.getIdCyllMPI(idCell)
            if len(idCylsAll) > 0:
                mol_rhizo = self.getContentCyl(idCyls, idComp)
                if konz:
                    if isDissolved:
                        V_rhizo = self.getWaterVolumesCyl(idCyls)/1e6
                    else:
                        V_rhizo = self.getVolumesCyl(idCyls)/1e6#sum(np.array([sum(cc.getCellSurfacesCyl()*self.seg_length[idCyls[i]] )/1e6 for i, cc in enumerate(cyls_i) ])) # m3
                    V_tot =  V_rhizo 
                else:
                    V_tot = 1
                    V_rhizo = np.nan
                res_CC = ( mol_rhizo)/(V_tot)#mol/m3 or mol
                if res_CC < 0:
                    print("res_CC < 0:",idComp,  mol_rhizo,V_tot ,V_rhizo , V_root)
                    print(idCyls, idSegs, idCell)
                    raise Exception
                
                contentOrKonz[idCell] =  res_CC
        comm.Barrier()
        return contentOrKonz
        
    def checkSoilBC(self, outer_R_bc_wat, #cm3 wat
                    outer_R_bc_sol):#mol
        comm.Barrier()
        sizeSoilCell = self.soilModel.getCellVolumes() /1e6#m3   
        self.diffSoilData_abs = np.zeros((self.numComp+1, len(sizeSoilCell)))
        self.diffSoilData_rel = np.zeros((self.numComp+1, len(sizeSoilCell)))
        self.maxdiffSoilData_abs = np.zeros(self.numComp+1)
        self.maxdiffSoilData_rel = np.zeros(self.numComp+1)
        self.sumdiffSoilData_abs = np.zeros(self.numComp+1)
        self.sumdiffSoilData_rel = np.zeros(self.numComp+1)
        
        soilTheta_new = self.soilModel.getWaterContent() #cm3/cm3 or m3/m3
        soilWatVolNew = soilTheta_new * sizeSoilCell
        if rank == 0:
            self.diffSoilData_abs[0] = abs((self.soilWatVol_old + outer_R_bc_wat/1e6) - soilWatVolNew)
            self.diffSoilData_rel[0] = abs(self.diffSoilData_abs[0]/soilWatVolNew)*100     
            self.maxdiffSoilData_abs[0] = max(self.diffSoilData_abs[0])
            self.maxdiffSoilData_rel[0] = max(self.diffSoilData_rel[0])
            self.sumdiffSoilData_abs[0] = sum(self.diffSoilData_abs[0])
            self.sumdiffSoilData_rel[0] = sum(self.diffSoilData_rel[0])
            
        for nc in range(1,self.numComp+1):
            isDissolved = nc <= self.numFluidComp
            soilContentNew = self.soilModel.getContent(nc, isDissolved)
            if rank == 0:
                self.diffSoilData_abs[nc] = abs((self.soilContent_old[nc - 1] + outer_R_bc_sol[nc -1]) - soilContentNew)
                soilContentNew[np.where(soilContentNew == 0.)] = 1.
                self.diffSoilData_rel[nc] = abs(self.diffSoilData_abs[nc]/soilContentNew)*100
                
                self.maxdiffSoilData_abs[nc] = max(self.diffSoilData_abs[nc])
                self.maxdiffSoilData_rel[nc] = max(self.diffSoilData_rel[nc])
                self.sumdiffSoilData_abs[nc] = sum(self.diffSoilData_abs[nc])
                self.sumdiffSoilData_rel[nc] = sum(self.diffSoilData_rel[nc])
        comm.Barrier()
        
    def setEmptySoilVoxel(self, emptySoilVoxels):
        comm.Barrier()
        sizeSoilCell = self.soilModel.getCellVolumes()/1e6 #m3 
        self.soilTheta_old[emptySoilVoxels] = self.soilModel.getWaterContent()[emptySoilVoxels]
        self.soilWatVol_old[emptySoilVoxels] = self.soilTheta_old[emptySoilVoxels] * sizeSoilCell[emptySoilVoxels] #m3 water
        soilCC_new = np.array([self.soilModel.getSolution(i+1)*self.phaseDensity(i+1) for i in range(self.numComp)])# mol/m3, shape : [cell][numcomp]?
        soilvolumes_new = np.array([self.soilWatVol_old if i <self.numFluidComp else sizeSoilCell for i in range(self.numComp)])#m3 water or m3 scv
        for nc in range(self.numComp):
            self.soilvolumes_old[nc][emptySoilVoxels] = soilvolumes_new[nc][emptySoilVoxels] #m3 water or m3 scv
            self.soilContent_old[nc][emptySoilVoxels] = soilCC_new[nc][emptySoilVoxels] * self.soilvolumes_old[nc][emptySoilVoxels] # mol
        comm.Barrier()        
        
    def setSoilData(self, 
                    sourceWat=[],   # cm3/day
                    sourceSol=[],   # mol/day
                    dt=0):          # day
        """set soil values before bulk flow """
        sizeSoilCell = self.soilModel.getCellVolumes()/1e6 #m3        
        
        #self.soilXX_old = self.soilModel.getSolutionHead()
        self.soilTheta_old = self.soilModel.getWaterContent() #cm3/cm3 or m3/m3
        self.soilWatVol_old = self.soilTheta_old * sizeSoilCell #m3 water
        try:
            assert self.soilWatVol_old.shape == (len(sizeSoilCell),)
        except:
            print(rank, sizeSoilCell)
            print(self.soilTheta_old)
            print(self.soilWatVol_old)
            raise Exception
        #mol/m3 wat or mol/m3 soil
        #maybe use content wather than concentration
        soilCC_old = np.array([self.soilModel.getSolution(i+1)*self.phaseDensity(i+1) for i in range(self.numComp)])# mol/m3, shape : [cell][numcomp]?
        self.soilvolumes_old = np.array([self.soilWatVol_old if i <self.numFluidComp else sizeSoilCell for i in range(self.numComp)])#m3 water or m3 scv
        self.soilContent_old = soilCC_old * self.soilvolumes_old # mol
        
        
        try:
            assert self.soilContent_old.shape == (self.numComp, len(sizeSoilCell))
        except:
            print(self.soilModel.getSolution(1).shape)
            print(self.soilModel.molarDensityWat_m3)
            print(self.soilWatVol_old, self.soilTheta_old*sizeSoilCell )
            print(self.soilContent_old, self.soilContent_old.shape)
            raise Exception
            
        if (len(sourceWat) > 0) and (rank == 0):
            cellIds = np.fromiter(self.cell2seg.keys(), dtype=int)
            cellIds =  np.array([x for x in cellIds if x >= 0])#take out air segments

            for cellId in cellIds:
                self.soilWatVol_old[cellId] += sourceWat[cellId]/1e6 * dt #m3
                for i in range(self.numComp):
                    self.soilContent_old[i][cellId] += sourceSol[i][cellId] * dt # mol += mol/day * day
                    try:
                        assert self.soilContent_old[i].shape == (len(sizeSoilCell),)
                    except:
                        print(self.soilContent_old[i].shape)
                        print( sourceSol[i][cellId] , dt)
                        raise Exception
        comm.Barrier()
    
            
        
    def updateOld(self, lId,  cyl):
        """ update distribution of cylinder if it s volume has changed
        """
        gId = self.eidx[lId]
        cellId = self.seg2cell[gId]
        segsId =  np.array(self.cell2seg[cellId])
        hasNewSegs = max(segsId) >= len(self.cyls)
        oldPoints = np.array(cyl.getPoints()).flatten() # cm
        lb = self.logbase            
        a_in = self.radii[gId] # cm
        a_out = self.outer_radii[gId]  # cm
        # print(a_in, a_out)
        if ((self.seg2cell[gId] > 0) and (self.organTypes[gId] == 2)):
            points = np.logspace(np.log(a_in) / np.log(lb), np.log(a_out) / np.log(lb), self.NC, base = lb) # cm
        else:
            points = np.array([a_in,a_out ])
        if (len(points) != len(oldPoints)) or (max(abs(points - oldPoints)) > 1e-13):#new vertex distribution
            
            ## old shape:
            centersOld = np.array(cyl.getCellCenters()).flatten();assert len(centersOld)==9      # cm 
            volOld = self.getVolumes(oldPoints, self.seg_length_old[gId] )# cm3
            surfaceOld =   np.pi * (oldPoints[1:] ** 2 - oldPoints[:-1]**2) # cm2
            ## new shape: 
            volNew = self.getVolumes(points, self.seg_length[gId] ) # cm^3
            centersNew = (points[1:] + points[:-1])/2  # cm
            surfaceNew = np.pi * (points[1:] ** 2 - points[:-1]**2) # cm2      
            ## change ratio
            if hasNewSegs:# currently, only new segments can gain content from old segmetns.
                changeRatio = min(sum(volNew)/sum(volOld), 1.)# we migth have volNew > volOld if the gain by L increase is higher than loss via r_out decrease
            else:
                changeRatio = 1.
           
            try:
                assert ((changeRatio <= 1.) and (changeRatio > 0.))
            except:
                print('volNew', volNew,'volOld',volOld)
                print("points",oldPoints, points)
                print('lengths', self.seg_length_old[gId],self.seg_length[gId])
                print('radii', self.radii[gId], self.outer_radii[gId],self.outer_radii_old[i])
                raise Exception
            ##  old concentration:
            theta_old = cyl.getWaterContent() # cm3/cm3
            cOld =np.array( [np.array(cyl.getConcentration(nC+1, (nC < self.numFluidComp))) for nC in range(self.numComp)])   #mol/cm3 scv    
            #m_totOld = [sum(c) for c in cOld] # cm3 or mol
            ## new contents:
            theta_new = np.array([theta_old[cellnum]* changeRatio * volOld[cellnum] /volNew[cellnum] for cellnum in range(len(volNew))])
            
            if min(theta_new) < self.vg_soil.theta_R:
                print("updateOld: min(theta_new) < self.vg_soil.theta_R. min theta set to  theta_R")
                theta_new_ = theta_new
                theta_new[np.where(theta_new < self.vg_soil.theta_R)] = self.vg_soil.theta_R
                try:
                    assert sum(theta_new*volNew) <= sum(theta_old*volOld)
                except:
                    print("to much water in voxel after setting min theta to theta_R")
                    print(gId,lId, cellId, volOld,theta_old, theta_new, theta_new_ , volNew) # see what to do in this circumstances
                    raise Exception
            
            try:
                assert min(theta_new) >= self.vg_soil.theta_R
                assert max(theta_new) <= self.vg_soil.theta_S
                xNew = np.array([vg.pressure_head(nt, self.vg_soil) for nt in theta_new])
            except:
                print("x_old",x_old,volOld )
                print("xnew", cxNew[0],self.getVolumes(points, self.seg_length[i] ) )
                print("points",oldPoints, points,centersNew)
                print("theta",theta_new, theta_old)
                print(self.vg_soil.theta_R, self.vg_soil.theta_S)
                raise Exception
            try:
                assert abs(sum(theta_new*volNew)/sum(theta_old*volOld) - changeRatio) < 1e-13 
            except:
                print("volOld", volOld, sum(volOld), "volnew", volNew, sum(volNew), "ratio", volOld/volNew,"totratio",  changeRatio)
                print("old content", theta_old*volOld, sum(theta_old*volOld) )
                print("new content", theta_new*volNew,sum(theta_new*volNew))
                print("ratio content",  sum(theta_new*volNew)/sum(theta_old*volOld),np.array(theta_new*volNew)/np.array(theta_old*volOld) )
                print("theta new", theta_new)
                print("theta old",theta_old)
                print(abs(sum(theta_new*volNew)/sum(theta_old*volOld) - changeRatio))
                raise Exception    
                
            cwNew = np.array([ np.array([cOld[nc][cellnum]* changeRatio * (theta_old[cellnum]* volOld[cellnum]) /(theta_new[cellnum]*volNew[cellnum]) 
                    for cellnum in range(len(volNew))]) for nc in range(self.numFluidComp)])
            csNew = np.array([ np.array([cOld[nc][cellnum]* changeRatio * volOld[cellnum] /volNew[cellnum] 
                    for cellnum in range(len(volNew))]) for nc in range(self.numFluidComp, self.numComp)])
            cNew = np.concatenate((cwNew, csNew))
            for ii in range(self.numFluidComp):
                ccnew = cNew[ii]
                ccold = cOld[ii]
                try:
                    assert abs(sum(ccnew*theta_new*volNew)/sum(ccold*theta_old*volOld) - changeRatio) < 1e-13
                except:
                    print("volOld", volOld, sum(volOld), "volnew", volNew, sum(volNew), "ratio", volOld/volNew,"totratio",  changeRatio)
                    print("old content",ccold*theta_old*volOld, sum(ccold*theta_old*volOld) )
                    print("new content",  ccnew*theta_new*volNew,sum(ccnew*theta_new*volNew))
                    print("ratio content",  sum(ccnew*theta_new*volNew)/sum(ccold*theta_old*volOld),np.array(ccnew*theta_new*volNew)/np.array(ccold*theta_old*volOld) )
                    print("ccnew", ccnew)
                    print("ccold",ccold)
                    print(abs(sum(ccnew*theta_new*volNew)/sum(ccold*theta_old*volOld) - changeRatio))
                    raise Exception
            for ii in range(self.numFluidComp, self.numComp):
                ccnew = cNew[ii]
                ccold = cOld[ii]
                try:
                    assert abs(sum(ccnew*volNew)/sum(ccold*volOld) - changeRatio) < 1e-13
                except:
                    print("volOld", volOld, sum(volOld), "volnew", volNew, sum(volNew), "ratio", volOld/volNew,"totratio",  changeRatio)
                    print("old content",ccold*volOld, sum(ccold*volOld) )
                    print("new content",  ccnew*volNew,sum(ccnew*volNew))
                    print("ratio content",  sum(ccnew*volNew)/sum(ccold*volOld),np.array(ccnew*volNew)/np.array(ccold*volOld) )
                    print("ccnew", ccnew)
                    print("ccold",ccold)
                    print(abs(sum(ccnew*volNew)/sum(ccold*volOld) - changeRatio))
                    raise Exception
                    
            m_totOld = np.concatenate((np.array([sum(cx * theta_old*volOld) for cx in cOld[:2]]),np.array([sum(cx * volOld) for cx in cOld[2:]])))
            m_totNew = np.concatenate((np.array([sum(cx * theta_new*volNew) for cx in cwNew]),np.array([sum(cx * volNew) for cx in csNew])))
            molFrNew = np.concatenate((np.array([cx / (self.molarDensityWat_m3/1e6) for cx in cwNew]),np.array([cx / (self.bulkDensity_m3/1e6)  for cx in csNew]) ))
            
            
            self.cyls[lId] = self.initialize_dumux_nc_( gId, 
                                                        x = xNew,# cm
                                                        cAll = molFrNew, # mol/mol water or mol/mol scv
                                                        Cells = centersNew) # cm
            x_new = self.cyls[lId].getWaterVolumesCyl( self.seg_length[gId]) # cm3
            cxNew_ = [np.array(self.cyls[lId].getContentCyl(nC+1, (nC < self.numFluidComp), self.seg_length[gId])) for nC in range(self.numComp)]    
            #cxNew_.insert(0, x_new)        
            m_totNew_ = [sum(cx) for cx in cxNew_]
            
            # check
            
            try:
                # assert (abs(m_totNew_ - m_totOld) < 1e-13 ).all() #3rd assert better, if we ever want to allow changeRatio > 1.
                # assert (sum(x_new) <= sum(theta_old*volOld) ).all() #3rd assert better
                assert ((abs(m_totNew_/m_totOld - changeRatio) < 1e-13).all() and (abs(sum(x_new)/sum(theta_old*volOld) - changeRatio) < 1e-13))
            except:                
                print("len(points) != len(oldPoints)",len(points) != len(oldPoints),len(points), len(oldPoints))
                print("max(abs(points - oldPoints)) > 0.",max(abs(points - oldPoints)) > 0.,max(abs(points - oldPoints)) ,abs(points - oldPoints))
                print("check results", abs(m_totNew_/m_totOld - changeRatio), abs(sum(x_new)/sum(theta_old*volOld) - changeRatio))
                print("updateOld",self.seg2cell[gId] ,gId,lId, self.organTypes[gId])
                print("m_totOld",m_totOld,"prescibed",m_totNew, "actual", m_totNew_, 'diff',  np.array(m_totNew_) - np.array(m_totOld),  np.array(m_totNew_) - np.array(m_totNew))
                print('ratios, prescribed',sum(theta_new*volNew)/ sum(theta_old*volOld), 'actual', sum(x_new)/sum(theta_old*volOld))
                print("x_old", theta_old*volOld,sum(theta_old*volOld),"x_new" ,theta_new*volNew, sum(theta_new*volNew) ,"actual",x_new, sum(x_new))
                print("points, old",oldPoints, 'new', points,'actual',np.array(self.cyls[lId].getPoints()).flatten())#,centersNew)
                print("volumes", self.getVolumes(oldPoints, self.seg_length_old[gId] ) ,sum(self.getVolumes(oldPoints, self.seg_length_old[gId] ) ), 
                        self.getVolumes(points, self.seg_length[gId] ),sum(self.getVolumes(points, self.seg_length[gId] ))  )
                print('lengths', self.seg_length_old[gId],self.seg_length[gId])
                print('radii', self.outer_radii[gId],self.outer_radii_old[gId])
                print('moleFr, prescibed:', molFrNew,'actual',[self.cyls[lId].getSolution(nc+1) for nc in range(self.numComp)])
                print('newHead', xNew)
                print('actalHead', self.cyls[lId].getSolutionHead())
                raise Exception
                
            
                                                 
    def getVolumes(self,vertices, length):
        return   np.pi * (vertices[1:] ** 2 - vertices[:-1]**2) * length
    
    def get_Vol_leftoverI(self, idCell):# m3
        comm.Barrier()
        
        V_total = self.soilModel.getCellVolumes()# [m3] cell volume
        if rank == 0:
            V_total = V_total[idCell] /1e6
        V_total = comm.bcast(V_total,root = 0)
        V_rhizo = 0
        if idCell in self.cell2seg:
            idSegs = self.cell2seg[idCell]#all segments in the cell
            idCyls =   self.getIdCyllMPI(idCell, False)            
            if len(idCyls) > 0: # we have old segments
                V_rhizo = sum(np.array([sum(self.cyls[np.where(self.eidx == gId)].getCellSurfacesCyl()*self.seg_length[gId])/1e6 for gId in idCyls ])) #[m3]
        V_rhizo = sum(comm.bcast(comm.gather(V_rhizo,root =0),root = 0)  )      
        Vol_res = V_total - V_rhizo #m3
        if Vol_res <= 0.:
            if ((Vol_res > -1e-13) and (len(idSegs) == len(idCyls))): # rounding error probably
                Vol_res = 0.
            else:
                print("get_Vol_leftoverI")
                print(idCell, V_total,V_rhizo, V_total - V_rhizo)
                print( len(idSegs), len(idCyls))# how many old and new cylinders?
                raise Exception
              
        if False:#idCell  == 486:
            print("get_watVol_leftoverI",idCell)
            print("V_total",V_total, V_rhizo)  
        comm.Barrier()
        return Vol_res
        
    def get_watVol_leftoverI(self, idCell):# m3
        comm.Barrier()
        wat_total = self.soilWatVol_old
        if rank == 0:
            wat_total = wat_total[idCell] # [m3] water content before bulk flow
        wat_total = comm.bcast(wat_total,root = 0)
        wat_rhizo = 0
        idCyls = []
        if idCell in self.cell2seg:
            idSegs = self.cell2seg[idCell]#all segments in the cell
            idCyls =   self.getIdCyllMPI(idCell, False)   
            if len(idCyls) > 0: # we have old segments
                wat_rhizo = sum(np.array([sum(self.cyls[np.where(self.eidx == gId)].getWaterVolumesCyl(self.seg_length[gId]))/1e6 for lId, gId in enumerate(idCyls) ])) #[m3]
        wat_rhizo = sum(comm.bcast(comm.gather(wat_rhizo,root=0),root=0))
        watVol_res = wat_total - wat_rhizo #m3
        if watVol_res < 0.:
            if watVol_res > -1e-13: # rounding error probably
                watVol_res = 0.
            else:
                print("getXX_leftoverI")
                print(idCell, wat_total,wat_rhizo, wat_total - wat_rhizo)
                print( len(idSegs), len(idCyls))# how many old and new cylinders?
                raise Exception
        if False:#idCell  == 486:
            print("get_watVol_leftoverI",idCell)
            print("wat_total",wat_total, wat_rhizo, idCyls, watVol_res)
        comm.Barrier()        
        return watVol_res
        
    
    def getC_content_leftoverI(self, idCell, idComp):# mol
        comm.Barrier()
        isDissolved = (idComp <= self.numFluidComp)
        mol_total = self.soilContent_old[idComp - 1]
        if rank == 0:
            mol_total = mol_total[idCell] # [m3] water content before bulk flow
        mol_total = comm.bcast(mol_total,root = 0)
        mol_rhizo = 0
        
        if idCell in self.cell2seg:
            idSegs = self.cell2seg[idCell]#all segments in the cell
            idCyls =  self.getIdCyllMPI(idCell, False)   
            
            if len(idCyls) > 0:#we have old segments
                mol_rhizo = sum(np.array([sum(self.cyls[np.where(self.eidx == gId)].getContentCyl(idComp, isDissolved, self.seg_length[gId])) for lId, gId in enumerate(idCyls) ])) #[m3]
        
        res_CC = mol_total - mol_rhizo
        if res_CC < 0.:
            if res_CC > -1e-13: # rounding error probably 
                res_CC = 0.
            else:
                print("getC_content_leftoverI", idCell)
                print("res_CC = ",res_CC," < 0, idComp:",idComp,'mol_total',mol_total ,'mol_rhizo', mol_rhizo,"res_CC",res_CC ) 
                print( len(idSegs), len(idCyls))# how many old and new cylinders?
                raise Exception
        comm.Barrier()
        return res_CC

                
    def initialize_(self,gId,x,cc ):
        if ((self.seg2cell[gId] > 0) and (self.organTypes[gId] == 2)):
            if self.mode == "dumux":
                self.cyls.append(self.initialize_dumux_(gId, x[self.seg2cell[gId]], False, False))
            elif self.mode == "dumux_exact":
                self.cyls.append(self.initialize_dumux_(gId, x[self.seg2cell[gId]], True, False))
            elif self.mode == "dumux_dirichlet":
                self.cyls.append(self.initialize_dumux_(gId, x[self.seg2cell[gId]], False, True))
            elif self.mode == "dumux_dirichlet_2c":
                self.cyls.append(self.initialize_dumux_nc_(gId, x[self.seg2cell[gId]], cc[self.seg2cell[gId]]))
            elif (self.mode == "dumux_dirichlet_nc") or(self.mode == "dumux_dirichlet_10c") or (self.mode == "dumux_10c") :
                self.cyls.append(self.initialize_dumux_nc_(gId, x[self.seg2cell[gId]], cc[self.seg2cell[gId]]))
            elif self.mode == "dumux_nc":
                self.cyls.append(self.initialize_dumux_(gId, x[self.seg2cell[gId]], False, True))
            elif self.mode == "python" or self.mode == "python_exact":
                x0 = x[self.seg2cell[gId]]
                self.cyls.append(self.initialize_python_(gId, x0))
            else:
                print("let s see", self.mode,"dumux_dirichlet_2c",self.mode == "dumux_dirichlet_2c" )
                raise Exception("RhizoMappedSegments.initialize: unknown solver {}".format(self.mode))
            if False:#self.seg2cell[gId]  == 486:
                print("initialize_",gId,x,cc, sum(self.cyls[-1].getWaterVolumesCyl(self.seg_length[gId])))
        else:
            a_in = self.radii[gId]#get Perimeter instead? not used for now anyway
            a_out = self.outer_radii[gId]
            self.cyls.append(AirSegment(a_in, a_out)) #psi_air computed by the photosynthesis module.
            #print("seg",gId,"is out of the soil domain and has no rhizosphere")        
    

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

        a_in = self.radii[gId]
        a_out = self.outer_radii[gId]
        lId =int( np.where(self.eidx == gId)[0])
        print('to wrap', rank,gId, a_in,a_out )
        if a_in < a_out:
            if self.mode == "dumux_dirichlet_2c":
                cyl = RichardsNoMPIWrapper(Richards2CCylFoam())  # only works for RichardsCylFoam compiled without MPI
            if self.mode == "dumux_dirichlet_nc":
                cyl = RichardsNoMPIWrapper(RichardsNCCylFoam())  # only works for RichardsCylFoam compiled without MPI
            if self.mode == "dumux_dirichlet_10c":
                cyl = RichardsNoMPIWrapper(Richards10CCylFoam(), self.useMoles)  # only works for RichardsCylFoam compiled without MPI
            if self.mode == "dumux_10c":
                cyl = RichardsNoMPIWrapper(Richards10CCylFoam(), self.useMoles)  # only works for RichardsCylFoam compiled without MPI
            cyl.initialize(verbose = True)
            cyl.setVGParameters([self.soil])
            lb = self.logbase
            
            if self.recreateComsol:
                nCells = self.NC
                cyl.createGrid([0.02], [0.6], [nCells])# cm
            else:
                points = np.logspace(np.log(a_in) / np.log(lb), np.log(a_out) / np.log(lb), self.NC, base = lb)
                cyl.createGrid1d(points)# cm
                if len(self.dx2) > lId:
                    self.dx2[lId] = 0.5 * (points[1] - points[0]) #when updating
                else:
                    self.dx2.append(0.5 * (points[1] - points[0]))#what is this?
            # print(a_in,a_out,lb)
            # print( a_in)
            if False:#i==1:
                cyl.setParameter("Problem.verbose", "0")
            else:
                cyl.setParameter("Problem.verbose", "-1")
            cyl.setParameter("Problem.reactionExclusive", "0")    
            cyl.setParameter( "Soil.Grid.Cells",str( self.NC-1)) # -1 to go from vertices to cell (dof)
            if self.recreateComsol:
                cyl.setHomogeneousIC(-100.)  # cm pressure head
            else:
                # cyl.setHomogeneousIC(x)  # cm pressure head
                #print("Soil.IC.P", cyl.dumux_str(x), Cells)
                cyl.setParameter("Soil.IC.P", cyl.dumux_str(x))# cm
            # cyl.setICZ_solute(c)  # [kg/m2] 
            
            #default: no flux
            cyl.setInnerBC("fluxCyl", 0.)  # [cm/day] #Y 0 pressure?
            #cyl.setInnerBC_solute("fluxCyl", 0.)  # [kg/m2], that s fair
            cyl.setOuterBC("fluxCyl", 0.)
            #cyl.setOuterBC_solute("fluxCyl", 0.)
            
            
            cyl.setParameter( "Soil.MolarMass", str(self.soilModel.solidMolarMass))
            cyl.setParameter( "Soil.solidDensity", str(self.soilModel.solidDensity))
            
            cyl.setParameter("Soil.betaC", str(self.soilModel.betaC ))
            cyl.setParameter("Soil.betaO", str(self.soilModel.betaO))
            cyl.setParameter("Soil.C_S_W_thresC", str(self.soilModel.C_S_W_thresC )) #mol/cm3
            cyl.setParameter("Soil.C_S_W_thresO", str(self.soilModel.C_S_W_thresO )) #mol/cm3
            cyl.setParameter("Soil.k_decay", str(self.soilModel.k_decay))
            cyl.setParameter("Soil.k_decay2", str(self.soilModel.k_decay2 ))
            #cyl.setParameter("Soil.k_decay3", str(1 ))
            cyl.setParameter("Soil.k_DC", str(self.soilModel.k_DC  )) # 1/d
            cyl.setParameter("Soil.k_DO", str(self.soilModel.k_DO  )) # 1/d
            #cyl.setParameter("Soil.k_growthBis", str(1 )) #bool
            cyl.setParameter("Soil.k_growthC", str(self.soilModel.k_growthC))
            cyl.setParameter("Soil.k_growthO", str(self.soilModel.k_growthO))
            cyl.setParameter("Soil.K_L", str(self.soilModel.K_L))#[mol/cm3]
            cyl.setParameter("Soil.k_phi", str(self.soilModel.k_phi ))
            cyl.setParameter("Soil.k_RC", str(self.soilModel.k_RC))
            cyl.setParameter("Soil.k_RO", str(self.soilModel.k_RO ))

            cyl.setParameter("Soil.k_SC", str(self.soilModel.k_SC )) #cm^3/mol/d
            #cyl.setParameter("Soil.k_SCBis", str(k_SC )) #cm^3/mol/d
            cyl.setParameter("Soil.k_SO", str(self.soilModel.k_SO )) #cm^3/mol/d
            #cyl.setParameter("Soil.k_SOBis", str(k_SO )) #cm^3/mol/d
 
            cyl.setParameter("Soil.m_maxC", str(self.soilModel.m_maxC  ))# 1/d
            cyl.setParameter("Soil.m_maxO", str(self.soilModel.m_maxO  ))# 1/d
            cyl.setParameter("Soil.micro_maxC", str(self.soilModel.micro_maxC ))# 1/d
            cyl.setParameter("Soil.micro_maxO", str(self.soilModel.micro_maxO ))# 1/d
            cyl.setParameter("Soil.v_maxL", str(self.soilModel.v_maxL))#[d-1]

            cyl.setParameter("Soil.k_sorp", str(self.soilModel.k_sorp)) # mol / cm3
            cyl.setParameter("Soil.f_sorp", str(self.soilModel.f_sorp)) #[-]
            cyl.setParameter("Soil.CSSmax", str(self.soilModel.CSSmax)) #[mol/cm3 scv]
            cyl.setParameter("Soil.alpha", str(self.soilModel.alpha)) #[1/d]


            cyl.setParameter( "Soil.BC.Bot.C1Type", str(3)) #put free flow later
            cyl.setParameter( "Soil.BC.Top.C1Type", str(3))
            cyl.setParameter( "Soil.BC.Bot.C1Value", str(0)) 
            cyl.setParameter( "Soil.BC.Top.C1Value", str(0 )) 

            cyl.setParameter("1.Component.LiquidDiffusionCoefficient", str(self.soilModel.Ds)) #m^2/s

            cyl.setParameter( "Soil.BC.Bot.C2Type", str(3))
            cyl.setParameter( "Soil.BC.Top.C2Type", str(3))
            cyl.setParameter( "Soil.BC.Bot.C2Value", str(0)) 
            cyl.setParameter( "Soil.BC.Top.C2Value", str(0 )) 
            cyl.setParameter("2.Component.LiquidDiffusionCoefficient", str(self.soilModel.Dl)) #m^2/s

            # cyl.decay = 0. #1.e-5
            
            for j in range(self.numFluidComp + 1, self.numComp+1):
                cyl.setParameter( "Soil.BC.Bot.C"+str(j)+"Type", str(3))
                cyl.setParameter( "Soil.BC.Top.C"+str(j)+"Type", str(3))
                cyl.setParameter( "Soil.BC.Bot.C"+str(j)+"Value", str(0)) 
                cyl.setParameter( "Soil.BC.Top.C"+str(j)+"Value", str(0 )) 
                
            for j in range( 1, self.numComp+1):       
                #print("cAll[j-1]",j-1, cAll[j-1],cyl.dumux_str(cAll[j-1])  )
                cyl.setParameter("Soil.IC.C"+str(j), cyl.dumux_str(cAll[j-1]) ) 

            if len(Cells) > 0:#in case we update a cylinder, to allow dumux to do interpolation
                assert(len(cAll[j-1])==len(Cells))
                CellsStr = cyl.dumux_str(Cells/100)#cm -> m
                cyl.setParameter("Soil.IC.Z",CellsStr)# m
                if len(Cells)!= len(x):
                    print("Cells, x",Cells, x, len(Cells), len(x))
                    raise Exception
                for j in range( 1, self.numComp+1):
                    cyl.setParameter("Soil.IC.C"+str(j)+"Z",CellsStr) # m
                    if len(Cells)!= len( cAll[j-1]):
                        print("Cells,  cAll[j-1]",Cells,  cAll[j-1], 
                                len(Cells), len(cAll[j-1]), j)
                        raise Exception
                        
            cyl.setParameter("Newton.EnableAbsoluteResidualCriterion", "True")
            cyl.setParameter("Newton.MaxAbsoluteResidual", "1.e-10")
                        
            cyl.setParameter("Problem.EnableGravity", "false")
            cyl.setParameter("Flux.UpwindWeight", "0.5")
            cyl.setParameter("Newton.EnableAbsoluteResidualCriterion", "true")
            cyl.setParameter("Newton.MaxAbsoluteResidual", "1.e-10")
            cyl.setParameter("Newton.EnableChop", "true")
            cyl.setParameter("Newton.EnableResidualCriterion", "true")
            cyl.setParameter("Newton.EnableShiftCriterion", "true")
            cyl.setParameter("Newton.MaxAbsoluteResidual", "1e-10")

            cyl.setParameter("Newton.MaxRelativeShift", "1e-10")

            cyl.setParameter("Newton.MaxSteps", "30")
            cyl.setParameter("Newton.ResidualReduction", "1e-10")
            cyl.setParameter("Newton.SatisfyResidualAndShiftCriterion", "true")
            cyl.setParameter("Newton.TargetSteps", "10")
            cyl.setParameter("Newton.UseLineSearch", "false")
            cyl.setParameter("Newton.EnablePartialReassembly", "true")
            cyl.setParameter("Grid.Overlap", "0")  #no effec5
            cyl.initializeProblem()
            cyl.setCriticalPressure(self.soilModel.wilting_point)  # cm pressure head
            cyl.bulkDensity_m3 = self.soilModel.bulkDensity_m3
            cyl.solidDensity =self.soilModel.solidDensity 
            cyl.solidMolarMass =self.soilModel.solidMolarMass
            cyl.solidMolDensity =self.soilModel.solidMolDensity         
            
            cyl.gId = gId    
            
            assert len(cyl.getSolutionHead()) == (self.NC - 1)
            return cyl
        else:
            print("RhizoMappedSegments.initialize_dumux_: Warning, segment {:g} might not be in domain, radii [{:g}, {:g}] cm".format(i, a_in, a_out))
            return []

    def initialize_dumux_(self, i, x, exact, dirichlet = False):
        """ Dumux RichardsCylFoam solver"""
        raise Exception
        a_in = self.radii[i]
        a_out = self.outer_radii[i]
        if a_in < a_out:
            cyl = RichardsNoMPIWrapper(RichardsCylFoam())  # only works for RichardsCylFoam compiled without MPI
            cyl.initialize()
            lb = self.logbase
            points = np.logspace(np.log(a_in) / np.log(lb), np.log(a_out) / np.log(lb), self.NC, base = lb)
            cyl.createGrid1d(points)
            self.dx2.append(0.5 * (points[1] - points[0]))
            cyl.setHomogeneousIC(x)  # cm pressure head
            cyl.setVGParameters([self.soil])
            cyl.setOuterBC("fluxCyl", 0.)  # [cm/day]
            if exact:
                cyl.setInnerBC("rootSystemExact")  # parameters are passed with cyl.setRootSystemBC (defined in richards_cyl.hh)
            else:
                if dirichlet:
                    cyl.setInnerBC("pressure", 0.)  # [cm/day]
                else:
                    cyl.setInnerBC("fluxCyl", 0.)  # [cm/day]
            cyl.setParameter("Newton.EnableAbsoluteResidualCriterion", "True")
            cyl.setParameter("Newton.MaxAbsoluteResidual", "1.e-10")
            cyl.initializeProblem()
            cyl.setCriticalPressure(self.wilting_point)  # cm pressure head
            return cyl
        else:
            raise Exception("RhizoMappedSegments.initialize_dumux_: Warning, segment {:g} might not be in domain, radii [{:g}, {:g}] cm".format(i, a_in, a_out))
            return []

    def initialize_python_(self, i, x):
        """ Python home grown richards fv solver"""
        a_in = self.radii[i]
        a_out = self.outer_radii[i]
        if a_in < a_out:
            lb = self.logbase
            points = np.logspace(np.log(a_in) / np.log(lb), np.log(a_out) / np.log(lb), self.NC, base = lb)
            grid = FVGrid1Dcyl(points)
            cyl = rich.FVRichards1D(grid, self.soil)
            cyl.x0 = np.ones((self.NC - 1,)) * x
            return cyl
        else:
            raise Exception("RhizoMappedSegments.initialize_python_: Warning, segment {:g} might not be in domain, radii [{:g}, {:g}] cm".format(i, a_in, a_out))
            return []

    def set_phloem_flux(self, rs):
        """ we need XylemFlux for the kr and kx call back functions (for python)"""
        self.rs = rs

    def get_inner_heads(self, shift = 0, weather:dict={}):#        
        """ matric potential at the root surface interface [cm]"""
        #print(self.mode, weather)
        rsx = np.zeros((len(self.cyls),))
        for i, cyl in enumerate(self.cyls):  # run cylindrical models
            if isinstance(cyl, AirSegment):
                rsx[i] = cyl.get_inner_head(val = self.rs.getPsiAir(weather["ea"]/weather["es"], weather["TairC"]))  # [cm]
            elif self.mode.startswith("dumux"):
                rsx[i] = cyl.getInnerHead(shift)  # [cm] (in richards.py, then richards_cyl.hh)
            elif self.mode.startswith("python"):
                rsx[i] = cyl.get_inner_head()  # [cm]
            else:
                raise Exception("RhizoMappedSegments.get_inner_heads: Warning, mode {:s} unknown".format(self.mode))
        return self._map(self._flat0(comm.gather( rsx, root = 0))) # gathers and maps correctly

    def get_inner_solutes(self, shift = 0, compId = 1):
        """ matric potential at the root surface interface [mol/cm3]"""
        rsx = np.full(len(self.cyls),0.)
        isDissolved = (compId <= self.numFluidComp)
        for i, cyl in enumerate(self.cyls):  # run cylindrical models
            rsx[i] = cyl.getInnerSolutes(shift, compId, isDissolved)  # [cm]
            #print('rsx', rsx[i])
        return self._map(self._flat0(comm.gather(rsx, root = 0))).flatten()  # gathers and maps correctly

    def get_soil_k(self, rx):
        """ TODO """
        soil_k = np.zeros((len(self.cyls),))
        if self.mode.startswith("dumux"):
            for i, idx in enumerate(self.eidx):  # run cylindrical models
                rsx = self.cyls[i].getInnerHead()
                nidx = idx + 1  # segment index+1 = node index
                try:
                    soil_k[i] = ((vg.fast_mfp[self.vg_soil](rsx) - vg.fast_mfp[self.vg_soil](rx[nidx])) / (rsx - rx[nidx])) / self.dx2[i]
                except:
                    print(rsx, rx[nidx])
                    raise Exception
        else:
            raise Exception("RhizoMappedSegments.get_soil_k: Warning, mode {:s} not implemented or unknown".format(self.mode))
        return self._map(self._flat0(comm.gather(soil_k, root = 0)))

    def get_dx2(self):
        """ TODO doc me AND only for mode="dumux" yet (set in initialize)"""
        return self._map(self._flat0(comm.gather(self.dx2, root = 0)))

    def get_inner_fluxes(self):
        """ fluxes [cm3/day] at the root surface, i.e. inner boundary  """
        fluxes = np.zeros((len(self.cyls),))#
        for i, cyl in enumerate(self.cyls):  # run cylindrical models
            if isinstance(cyl, AirSegment):   
                j = self.eidx[i]  # for one process j == i
                fluxes[i] = float(cyl.getInnerFlux()) 
            elif self.mode.startswith("dumux"):            
                j = self.eidx[i]  # for one process j == i
                fluxes[i] = -float(cyl.getInnerFlux()) * (2 * np.pi * self.radii[j] * self.seg_length[j]) / self.radii[j]  # [cm/day] -> [cm3/day], ('/inner_radii' comes from cylindrical implementation)
            elif self.mode.startswith("python"):
                j = self.eidx[i]  # for one process j == i
                fluxes[i] = cyl.get_inner_flux() * (2 * np.pi * self.radii[j] * self.seg_length[j]) / self.last_dt  # divide by dt is correct here! getInnerFlux only gives the source in cm3/cm2
            else:
                raise Exception("RhizoMappedSegments.get_inner_fluxes: Warning, mode {:s} unknown".format(self.mode))
        return self._map(self._flat0(comm.gather(fluxes, root = 0)))  # gathers and maps correctly

    def get_concentration(self, compId, konz = True):  
        """ solute or water concentration """
        isAirSegs = np.array(list(set(np.concatenate((self.cell2seg[-1],np.where(np.array(self.organTypes) != 2)[0])) )))#np.array([isinstance(cc , AirSegment) for cc in self.cyls])
        if self.mode.startswith("dumux"):
            if compId == 0:
                V_rhizo = self.getVolumesCyl(doSum = False, reOrder = True)#np.array([sum(cc.getCellSurfacesCyl()*self.seg_length[i]) for i, cc in enumerate(self.cyls) ]) #[cm3]
                contentRhizo = self.getWaterVolumesCyl(doSum = False, reOrder = True)#np.array([sum(cc.getWaterVolumesCyl(self.seg_length[i])) for i, cc in enumerate(self.cyls) ]) # cm3
            else:
                isDissolved = (compId <= self.numFluidComp)
                if isDissolved: 
                    V_rhizo = self.getWaterVolumesCyl(doSum = False, reOrder = True)#np.array([sum(cc.getWaterVolumesCyl(self.seg_length[i])) for i, cc in enumerate(self.cyls) ]) # cm3     
                    V_rhizo[isAirSegs] = 1 # avoid RuntimeWarning: invalid value encountered in true_divide                    
                else:
                    V_rhizo = self.getVolumesCyl(doSum = False, reOrder = True)#np.array([sum(cc.getCellSurfacesCyl()*self.seg_length[i]) for i, cc in enumerate(self.cyls) ]) #[cm3]                
                contentRhizo = self.getContentCyl(idComp=compId, doSum = False, reOrder = True)#np.array([sum(cc.getContentCyl(compId, isDissolved, self.seg_length[i])) for i, cc in enumerate(self.cyls) ]) # mol
            assert (contentRhizo[isAirSegs] == 0.).all()
            if not konz:
                return contentRhizo
            concentration = contentRhizo/V_rhizo
            concentration[isAirSegs] = np.nan
            #print("organTypes",self.organTypes)
            #print("concentration",concentration)
            return concentration # mol/cm3
        else:
            raise Exception("RhizoMappedSegments.get_inner_concentrations: Warning, mode {:s} unknown".format(self.mode))
        return self._map(self._flat0(comm.gather(rsx, root = 0)))  # gathers and maps correctly

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
        return self._map(self._flat0(comm.gather(rsx, root = 0)))  # gathers and maps correctly

    def get_inner_mass_fluxes(self):  # TODO check
        """ fluxes [cm3/day] at the root surface, i.e. inner boundary  """
        fluxes = np.zeros((len(self.cyls),))
        if self.mode.startswith("dumux"):
            for i, cyl in enumerate(self.cyls):  # run cylindrical models
                j = self.eidx[i]  # for one process j == i
                fluxes[i] = -float(cyl.getInnerFlux(1)) * (2 * np.pi * self.radii[j] * self.seg_length[j]) / self.radii[j]  # [cm/day] -> [cm3/day], ('/inner_radii' comes from cylindrical implementation)
        else:
            raise Exception("RhizoMappedSegments.get_inner_fluxes: Warning, mode {:s} unknown".format(self.mode))
        return self._map(self._flat0(comm.gather(fluxes, root = 0)))  # gathers and maps correctly

    def get_water_volumes(self):
        """ fluxes [cm3/day] at the root surface, i.e. inner boundary  """
        volumes = np.zeros((len(self.cyls),))
        if self.mode.startswith("dumux"):
            for i, cyl in enumerate(self.cyls):  # run cylindrical models
                j = self.eidx[i]  # for one process j == i
                volumes[i] = cyl.getWaterVolume()
        else:
            raise Exception("RhizoMappedSegments.get_water_volumes: Warning, mode {:s} unknown".format(self.mode))
        return self._map(self._flat0(comm.gather(volumes, root = 0)))  # gathers and maps correctly

    def solve(self, dt, *argv):#dt, seg_rx, proposed_outer_fluxes,kex,proposed_outer_sol_fluxes
        """ set bc for cyl and solves it """
        self.last_dt = dt

        if ((self.mode == "dumux_nc") or (self.mode == "dumux_10c")):
            #inner = bot = plant
            #outer = top = soil
            proposed_inner_fluxes = argv[0]
            proposed_outer_fluxes = argv[1]


            
            for lId, cyl in enumerate(self.cyls):  # run cylindrical models
                gId = self.eidx[lId]  # for one process gId == lId
                print(rank, "cyl no ",lId,"/",len(self.cyls),'gId')
                if isinstance(cyl, AirSegment):  
                    cyl.setInnerFluxCyl(proposed_inner_fluxes[gId])
                else:
                    l = self.seg_length[gId]
                    if self.recreateComsol:
                        raise Exception
                        cyl.setInnerFluxCyl(-0.26)#proposed_inner_fluxes[j] / (2 * np.pi * self.radii[j] * l))  # [cm3/day] -> [cm /day]
                        cyl.setOuterFluxCyl(0)#proposed_outer_fluxes[j] / (2 * np.pi * self.outer_radii[j] * l))  # [cm3/day] -> [cm /day]
                    else:   
                        QflowIn = proposed_inner_fluxes[gId] 
                        qIn = QflowIn/ (2 * np.pi * self.radii[gId] * l) # [cm3/day] -> [cm /day]
                        QflowOut = proposed_outer_fluxes[gId] 
                        qOut = QflowOut/(2 * np.pi * self.outer_radii[gId] * l)# [cm3/day] -> [cm /day]
                        cyl.setInnerFluxCyl(qIn)  
                        cyl.setOuterFluxCyl(qOut)
                         
                    if not isinstance(argv[2],float): 
                        # cyl.setInnerBC_solute("constantFluxCyl", argv[2][i])
                        botVal = argv[2][gId]
                        topVal = argv[3][gId]
                        botVal_mucil = argv[4][gId]
                        topVal_mucil = argv[5][gId]
                    else:
                        # cyl.setInnerBC_solute("constantFluxCyl", argv[2])
                        botVal = argv[2]
                        topVal = argv[3]
                        botVal_mucil = argv[4]
                        topVal_mucil = argv[5]
                        
                    typeBC = np.full(self.numComp,3)
                    
                    valueTopBC = np.full(self.numComp,0.)
                    valueBotBC = np.full(self.numComp,0.)
                    if self.recreateComsol:
                        valueBotBC[0] = 1.; valueBotBC[1] = 1.
                    else :
                        #print("solutes_bot",  botVal / (2 * np.pi * self.radii[j] * l), 
                        #    botVal / (2 * np.pi * self.radii[j] * l) )
                        valueBotBC[0] = botVal / (2 * np.pi * self.radii[gId] * l) # [mol/day] -> [mol/cm2 /day]
                        valueBotBC[1] = botVal_mucil / (2 * np.pi * self.radii[gId] * l) # [mol/day] -> [mol/cm2 /day]
                        
                        valueTopBC[0] = topVal / (2 * np.pi * self.outer_radii[gId] * l) # [mol/day] -> [mol/cm2 /day]
                        valueTopBC[1] = topVal_mucil / (2 * np.pi * self.outer_radii[gId] * l) # [mol/day] -> [mol/cm2 /day]
                        
                    #print("botSBC",typeBC, valueBotBC) 
                    cyl.setSoluteBotBC(typeBC, valueBotBC)
                    cyl.setSoluteTopBC(typeBC, valueTopBC)
                    # print("inner flux", proposed_inner_fluxes[j] / (2 * np.pi * self.radii[j] * l), valueBotBC)
                    # print("outer flux",proposed_outer_fluxes[j] / (2 * np.pi * self.outer_radii[j] * l), valueTopBC)
                    
                    # write_file_array("pressureHead_no"+str(i),np.array(cyl.getSolutionHead()).flatten())
                    # write_file_array("coord_no"+str(i), cyl.getDofCoordinates().flatten())
                    # for jj in range(self.numFluidComp):
                        # write_file_array("solute_conc"+str(jj+1)+"_no"+str(i), 
                                        # np.array(cyl.getSolution(jj+1)).flatten()* self.molarDensityWat ) 
                    # for jj in range(self.numFluidComp, self.numComp):
                        # write_file_array("solute_conc"+str(jj+1)+"_no"+str(i), 
                                        # np.array(cyl.getSolution(jj+1)).flatten()* self.bulkDensity_m3 /1e6 ) 
                    
                
                    # print("summary", i, j, botVal, topVal, proposed_inner_fluxes[j],proposed_outer_fluxes[j],
                    #         l, self.outer_radii[j] , self.radii[j])
                    buWBefore_ = cyl.getWaterVolumesCyl(l)
                    buWBefore = sum( buWBefore_ )
                    WaterContentOld = cyl.getWaterContent()
                    SolutionHeadOld = cyl.getSolutionHead()
                    buCCBefore_ = cyl.getContentCyl(1, True, l)
                    buTotCBefore = self.getTotCContent(cyl,l)
                    solsBefore = [  np.array(cyl.getSolution(ooo + 1)).flatten() for ooo in range(self.numComp)]
                    css1_before = np.array(cyl.base.getCSS1_out()).flatten()/1e6 

                    try:
                        print(rank, lId, 'start solve')
                        cyl.solve(dt, maxDt = 2500/(24*3600))
                        print(rank, lId, 'end solve')
                        buWAfter_ =  cyl.getWaterVolumesCyl(l)
                        buWAfter = sum(buWAfter_ )
                        buCCAfter_ = cyl.getContentCyl(1, True, l)
                        buTotCAfter = self.getTotCContent(cyl,l)
                        
                        try:
                            assert abs(buWAfter - ( buWBefore + (QflowIn + QflowOut) * dt)) < 1e-5
                        except:#check water mass balance
                            print("before",buWBefore_ )
                            print("after",buWAfter_ )
                            print("water conservation ",lId,gId, buWAfter ,"before", buWBefore )
                            print("QflowIn", QflowIn ,"QflowOut", QflowOut,"dt", dt)
                            print( "l",l,"a_in", self.radii[gId] ,"a_out",self.outer_radii[gId] )
                            print("diff",( buWBefore + (QflowIn + QflowOut) * dt))
                            #print("issue water conservation ",buWAfter , buWBefore)
                            print("qflowout",qOut , "qflowin",qIn)
                            print("theta", cyl.getWaterContent())
                            raise Exception
                        try:#check solute mass balance
                            assert abs((sum(buTotCAfter) - ( sum(buTotCBefore) + (botVal + topVal + botVal_mucil+ topVal_mucil) * dt))/sum(buTotCAfter)*100) < 1.
                        except:
                            print("error mass C")
                            print("before",buTotCBefore )
                            print("after",buTotCAfter )
                            print("CC conservation ",lId, sum(buTotCAfter) ,"before", sum(buTotCBefore) )
                            print("QCflowIn", botVal, botVal_mucil ,"QCflowOut", topVal, topVal_mucil,"dt", dt)
                            print("qCIn",valueBotBC,"qCOut",valueTopBC)
                            print( "l",l,"a_in", self.radii[gId] ,"a_out",self.outer_radii[gId] )
                            print("diff",( sum(buTotCAfter) -( sum(buTotCBefore) + (botVal + topVal + botVal_mucil+ topVal_mucil) * dt))/sum(buTotCAfter)*100)
                            #print("issue water conservation ",buWAfter , buWBefore)
                            #print("qflowout",qOut , "qflowin",qIn)
                            #print("theta", cyl.getWaterContent())
                            for ncomp in range(self.numComp):
                                print("ncomp_before", ncomp +1, solsBefore[ncomp])
                                print("ncomp_after",ncomp +1 , np.array(cyl.getSolution(ncomp + 1)).flatten() )
                            print("css1 before", css1_before)
                            print("css1 after",  np.array(cyl.base.getCSS1_out()).flatten()/1e6 )
                            raise Exception
                            
                        for ncomp in range(self.numComp):
                            try:
                                assert (np.array(cyl.getSolution(ncomp + 1)).flatten() >= 0).all()
                            except:
                                print("CC <0:",lId,gId,ncomp, np.array(cyl.getSolution(ncomp + 1)).flatten() )
                                print(valueBotBC,valueTopBC)
                                print(buCCBefore_,sum(buCCBefore_))
                                print(buCCAfter_, sum(buCCAfter_))
                                raise Exception
                        
                    except:
                        print( "l",l,"a_in", self.radii[gId] ,"a_out",self.outer_radii[gId] )
                        print("water conservation ",lId, "before", buWBefore ,"buWAfter ",buWAfter , "all cells",buWBefore_,buWAfter_)
                        print("QWflowIn", QflowIn ,"QWflowOut", QflowOut,"qWflowIn", qIn ,"qWflowOut", qOut,"dt", dt)
                        print("theta_old", WaterContentOld)
                        print("head_old", SolutionHeadOld)
                        print("theta_new", cyl.getWaterContent())
                        print("head_new", cyl.getSolutionHead())
                        print("CC conservation ",lId,"before", sum(buTotCBefore) )
                        print("QflowIn", botVal, botVal_mucil ,"QflowOut", topVal, topVal_mucil,"dt", dt)
                        print("concentration", [cyl.getSolution(i+1) for i in range(self.numComp)])
                        myStr = "RhizoMappedSegments.solve: dumux exception with boundaries in flow {:g} cm3/day, out flow {:g} cm3/day, segment radii [{:g}-{:g}] cm"
                        myStr = myStr.format(proposed_inner_fluxes[gId] / (2 * np.pi * self.radii[gId] * l), proposed_outer_fluxes[gId] / (2 * np.pi * self.radii[gId] * l), self.radii[gId], self.outer_radii[gId])

                        
                        raise Exception(myStr)
                
                    # write_file_array("pressureHead_no"+str(i),np.array(cyl.getSolutionHead()).flatten())
                    # write_file_array("coord_no"+str(i), cyl.getDofCoordinates().flatten())
                    # for jj in range(self.numFluidComp):
                        # write_file_array("solute_conc"+str(jj+1)+"_no"+str(i), 
                                        # np.array(cyl.getSolution(jj+1)).flatten()* self.molarDensityWat ) 
                    # for jj in range(self.numFluidComp, self.numComp):
                        # write_file_array("solute_conc"+str(jj+1)+"_no"+str(i), 
                                        # np.array(cyl.getSolution(jj+1)).flatten()* self.bulkDensity_m3 /1e6 ) 
                
                print(rank, "FINISHED cyl no ",lId,"/",len(self.cyls),'gId', isinstance(cyl, AirSegment))
                
                
        elif self.mode == "dumux":
            proposed_inner_fluxes = argv[0]
            proposed_outer_fluxes = argv[1]
            for i, cyl in enumerate(self.cyls):  # run cylindrical models
                j = self.eidx[i]  # for one process j == i
                l = self.seg_length[j]
                if isinstance(cyl, AirSegment):  
                    cyl.setInnerFluxCyl(proposed_inner_fluxes[j])
                else:                
                    cyl.setInnerFluxCyl(proposed_inner_fluxes[j] / (2 * np.pi * self.radii[j] * l))  # [cm3/day] -> [cm /day]
                    cyl.setOuterFluxCyl(proposed_outer_fluxes[j] / (2 * np.pi * self.outer_radii[j] * l))  # [cm3/day] -> [cm /day]
                    try:
                        cyl.solve(dt)
                    except:
                        myStr = "RhizoMappedSegments.solve: dumux exception with boundaries in flow {:g} cm3/day, out flow {:g} cm3/day, segment radii [{:g}-{:g}] cm"
                        myStr = myStr.format(proposed_inner_fluxes[j] / (2 * np.pi * self.radii[j] * l), proposed_outer_fluxes[j] / (2 * np.pi * self.radii[j] * l), self.radii[j], self.outer_radii[j])
                        print("node ", self.nodes[self.segments[j].y])
                        self.plot_cylinder(j)
                        self.plot_cylinders()
                        raise Exception(myStr)
        
        elif ((self.mode == "dumux_dirichlet_nc") or (self.mode == "dumux_dirichlet_2c")
                or (self.mode == "dumux_dirichlet_10c")):            
            rx = argv[0]
            proposed_outer_fluxes = argv[1]
            proposed_inner_fluxes = argv[2]
            for i, cyl in enumerate(self.cyls):  # run cylindrical models
                    
                j = self.eidx[i]  # for one process j == i
                if isinstance(cyl, AirSegment):  
                    cyl.setInnerFluxCyl(proposed_inner_fluxes[j])
                    #cyl.setParameter( "Soil.BC.Bot.C1Value", str(exud)) 
                else:
                    l = self.seg_length[j]
                    cyl.setInnerMatricPotential(rx[i])
                    #inner = bot = plant
                    #outer = top = soil
                    
                    cyl.setParameter( "Soil.BC.Top.C1Value", proposed_outer_fluxes[j] / (2 * np.pi * self.outer_radii[j] * l))  # [cm3/day] -> [cm /day]
                    
                    if not isinstance(argv[2],float): 
                        cyl.setInnerBC_solute("constantFluxCyl", argv[2][i])
                        #setBotBC_solute
                        cyl.setParameter(self.param_group + "BC.Bot.CValue", str(value_bot))
                    else:
                        cyl.setInnerBC_solute("constantFluxCyl", argv[2])
                    cyl.setOuterBC_solute("constantFluxCyl", argv[3])
                    try:
                        cyl.solve(dt)
                    except:
                        # str = "RhizoMappedSegments.solve: dumux exception with boundaries in flow {:g} cm3/day, out flow {:g} cm3/day, segment radii [{:g}-{:g}] cm"
                        # str = str.format(proposed_inner_fluxes[j] / (2 * np.pi * self.radii[j] * l), proposed_outer_fluxes[j] / (2 * np.pi * self.radii[j] * l), self.radii[j], self.outer_radii[j])
                        # print("node ", self.nodes[self.segments[j].y])
                        self.plot_cylinder(j)
                        self.plot_cylinders()
                        self.plot_cylinders_solute()
                        raise Exception
        else:
            print(self.mode)
            raise Exception("RhizoMappedSegments.initialize: unknown solver {}".format(self.mode))
        
    def getTotCContent(self, cyl, l):
        if isinstance(cyl, AirSegment):
            return np.full(self.numComp,0.)
        totC = 0
        vols = cyl.getCellSurfacesCyl()  * l  #cm3 scv
        for i in range(self.numComp):
            isDissolved = (i < 2)
            totC += cyl.getContentCyl(i+1, isDissolved,l)
        C_S_W = np.array(cyl.getSolution(1)).flatten()*self.molarDensityWat_m3
        
        init = (cyl.simTime == 0.)
        if init:
            css1 = self.soilModel.CSSmax * (C_S_W/(C_S_W+ self.soilModel.k_sorp)) * self.soilModel.f_sorp
        else:
            css1 = np.array(cyl.base.getCSS1_out()).flatten()/1e6 #mol C / cm3 scv
            assert css1[-1] == 0.
            css1 = css1[:-1]
        totC += css1 * vols
        assert len(np.array(totC).shape)==1
        return totC
    
    def get_water_volume(self):
        """ returns the water volume of the cylindrical models [cm3] """
        volumes = np.zeros((len(self.cyls),))
        if self.mode.startswith("dumux"):
            for i, cyl in enumerate(self.cyls):  # run cylindrical models
                j = self.eidx[i]  # for one process j == i
                cyl_water = 0.
                cyl_water_content = cyl.getWaterContent()  # getWaterContent() in cpp/pyhton_binding/richards.hh
                nodes = cyl.getPoints()
                for k, wc in enumerate(cyl_water_content):
                    r1 = nodes[k]
                    r2 = nodes[k + 1]
                    cyl_water += np.pi * (r2 * r2 - r1 * r1) * self.seg_length[j] * wc
                volumes[i] = cyl_water
        elif self.mode.startswith("python"):
            for i, cyl in enumerate(self.cyls):  # run cylindrical models
                j = self.eidx[i]  # for one process j == i
                cyl_water = 0.
                cyl_water_content = cyl.get_water_content()  # segment 0
                for k, wc in enumerate(cyl_water_content):
                    r1 = cyl.grid.nodes[k]
                    r2 = cyl.grid.nodes[k + 1]
                    cyl_water += np.pi * (r2 * r2 - r1 * r1) * self.seg_length[j] * wc
                volumes[i] = cyl_water
        else:
            raise Exception("RhizoMappedSegments.get_water_volume: unknown solver {}".format(self.mode))
        return self._map(self._flat0(comm.gather(volumes, root = 0)))  # gathers and maps correctly
    
    def splitSoilVals(self, soilVals, seg_values, troubleShootId = -1):
        """ split soilFlux array according to the values in seg_values """
        verbose = False
        if troubleShootId <0. :
            cellIds = np.fromiter(self.cell2seg.keys(), dtype=int)
            cellIds =  np.array([x for x in cellIds if x >= 0])#take out air segments
        else:
            cellIds = np.array([troubleShootId])
            verbose = True
            
        assert soilVals.shape == ( len(self.soilModel.getCellVolumes()), )
        assert seg_values.shape == (len(self.cyls), )
        # sum seg_values per voxel
        organTypes = np.array(self.organTypes)
        
        if (organTypes != 2).any():
            try:
                assert (seg_values[np.where(organTypes != 2)] == 0.).all()
            except:
                print(seg_values, organTypes)
                print(seg_values[np.where(organTypes != 2)])
                raise Exception
            
        
        seg_values_voxel = np.full(len(self.cyls), 0.)
        splitVals = np.full(len(self.cyls), 0.)
        
        for cellid in cellIds:
            segIds = self.cell2seg[cellid]
            assert (splitVals[segIds] == 0.).all()
            ots = organTypes[segIds]
            if (ots != 2).any():
                try:
                    assert (seg_values[segIds][np.where(ots != 2)] == 0.).all()
                except:
                    print(cellid, segIds, ots, seg_values[segIds])
                    raise Exception
            if verbose:
                print("soilVals[cellid]", cellid, soilVals[cellid], segIds)
            if soilVals[cellid] < 0:# goes away from the 1d models
                seg_values_voxel[segIds] = sum(seg_values[segIds])
                weightVals = seg_values[segIds] / seg_values_voxel[segIds]
                
                if verbose:
                    print("sum(seg_values[segIds])", seg_values_voxel[segIds], weightVals)
            else:# goes toward  the 1d models
                seg_values_voxel[segIds] = sum(1/seg_values[segIds])
                weightVals = (1 / seg_values[segIds]) / seg_values_voxel[segIds]
                if verbose:
                    print("sum(1/seg_values[segIds])", seg_values_voxel[segIds], weightVals)
            
            splitVals[segIds] = weightVals * soilVals[cellid]
            if verbose:
                print("splitVals[segIds]",splitVals[segIds])
                print("sum(weightVals)", sum(weightVals), sum(splitVals[segIds]))
            try:
                assert (sum(weightVals) - 1.) < 1e-13
                assert abs(sum(splitVals[segIds]) - soilVals[cellid]) < 1e-13
            except:
                print(weightVals, sum(weightVals))
                print(splitVals[segIds], soilVals[cellid])
                print(sum(splitVals[segIds]), soilVals[cellid])
                print(sum(splitVals[segIds]) - soilVals[cellid])
                raise Exception
                
        try:
            assert abs((sum(splitVals) - sum(soilVals[cellIds]))/sum(soilVals[cellIds])) < 0.1
        except:
            print(sum(splitVals), sum(soilVals),  sum(soilVals[cellIds]))
            print(splitVals,soilVals )
            print(sum(splitVals) -  sum(soilVals[cellIds]))
            splitVals_ = 0.
            soilVals_ = 0.
            for cellid in cellIds:
                segIds = self.cell2seg[cellid]
                splitVals_ += sum(splitVals[segIds])
                soilVals_ +=  soilVals[cellid]
                print(cellid,segIds,sum(splitVals[segIds]), soilVals[cellid],
                    "current",splitVals_,  soilVals_,splitVals_,-  soilVals_,
                    "to",sum(splitVals),sum(soilVals[cellIds]),sum(splitVals)-sum(soilVals[cellIds]),
                    (sum(splitVals)-sum(soilVals[cellIds]))/sum(soilVals[cellIds]),
                    sum(soilVals))
            raise Exception
            
        try:
            assert (splitVals[np.where(organTypes != 2)] == 0.).all()
        except:
            print("np.where(organTypes != 2)", np.where(organTypes != 2))
            print("splitVals",splitVals,organTypes )
            print(splitVals[np.where(organTypes != 2)] )
            print("seg_values",seg_values)
            print(seg_values[np.where(organTypes != 2)] )
            raise Exception
        return splitVals
    
    def _map(self, x):
        """Converts @param x to a numpy array and maps it to the right indices                 """
        indices = self._flat0(comm.gather(self.eidx, root = 0))  # gather segment indices from all threads
        if len(indices) > 0:  # only for rank 0 it is not empty
            assert len(indices) == len(x), "RhizoMappedSegments._map: indices and values have different length"
            p = np.zeros((len(x),), dtype = np.float64)
            for i in range(0, len(indices)):  #
                p[indices[i]] = x[i]            
            return p
        else:
            return np.array([])

    def _flat0(self, xx):
        """flattens the gathered list in rank 0, empty list for other ranks """
        #print("rhizoMappedPlant::_flat0", rank, xx)
        if rank == 0:
            return [item for sublist in xx for item in sublist]
        else:
            return []

    def plot_cylinder(self, i):
        """ plots a specific cylinder (DUMUX only, TODO) """
        cyl = self.cyls[i]
        x_ = cyl.getDofCoordinates()
        y_ = cyl.getSolutionHead()

        SMALL_SIZE = 22
        MEDIUM_SIZE = 22
        BIGGER_SIZE = 22
        plt.rc('font', size = SMALL_SIZE)  # controls default text sizes
        plt.rc('axes', titlesize = SMALL_SIZE)  # fontsize of the axes title
        plt.rc('axes', labelsize = MEDIUM_SIZE)  # fontsize of the x and y labels
        plt.rc('xtick', labelsize = SMALL_SIZE)  # fontsize of the tick labels
        plt.rc('ytick', labelsize = SMALL_SIZE)  # fontsize of the tick labels
        plt.rc('legend', fontsize = SMALL_SIZE)  # legend fontsize
        plt.rc('figure', titlesize = BIGGER_SIZE)  # fontsize of the figure title
        # plt.xlim([0., 0.2 ])
        plt.plot(x_, y_)
        plt.xlabel("distance [cm]")
        plt.ylabel("matric potential [cm]")
        plt.show()

    def plot_cylinders(self):
        """ plots a specific cylinder (DUMUX only, TODO) """
        inner, outer = [], []
        zz = -self.minBound.z
        for i, cyl in enumerate(self.cyls):
            if not isinstance(cyl, AirSegment):  
                x_ = cyl.getDofCoordinates()
                y_ = cyl.getSolutionHead()
                inner.append(y_[0])
                outer.append(y_[-1])
                j = self.segments[i].y
                z = self.nodes[j].z
                col_i = int(-z / zz * 255.)
                c_ = '#%02x%02x%02x' % (col_i, col_i, 64)
                plt.plot(x_, y_, alpha = 0.1, c = c_)
        plt.xlabel("distance [cm], deeper roots are yellow")
        plt.ylabel("matric potential [cm]")
        # plt.xlim([0.05, 0.6])
        # plt.ylim([-8500, 0. ])
        plt.show()
        return  np.argmin(inner), np.argmax(inner), np.argmin(outer), np.argmax(inner)

    def plot_cylinders_solute(self):
        """ plots a specific cylinder (DUMUX only, TODO) """
        inner, outer = [], []
        zz = -self.minBound.z
        for i, cyl in enumerate(self.cyls):
            if not isinstance(cyl, AirSegment):  
                x_ = cyl.getDofCoordinates()
                y_ = cyl.getSolution(1)
                #y_ = cyl.getWaterContent() works
                inner.append(y_[0])
                outer.append(y_[-1])
                j = self.segments[i].y
                z = self.nodes[j].z
                col_i = int(-z / zz * 255.)
                c_ = '#%02x%02x%02x' % (col_i, col_i, 64)
                plt.plot(x_-x_[0], y_, alpha = 0.1, c = c_)
        plt.xlabel("distance [cm], deeper roots are yellow")
        plt.ylabel('solute concentration (g/cm3)')
        # plt.xlim([0.05, 0.6])
        # plt.ylim([-8500, 0. ])
        plt.show()
        return  np.argmin(inner), np.argmax(inner), np.argmin(outer), np.argmax(inner)


    def map_cylinders_solute(self, XX,YY,ZZ,name = ""):
        """ maps cylinders to soil grid """

        shape = np.shape(XX)
        conc = np.zeros((shape))
        for i, cyl in enumerate(self.cyls):
            if not isinstance(cyl, AirSegment):  
                R_ = cyl.getDofCoordinates()
                vv_ = cyl.getSolution(1)
                
                p0 = np.array(self.nodes[self.segments[i].x])
                p1 = np.array(self.nodes[self.segments[i].y])
                
                v = p1 - p0
                mag = norm(v)
                v = v / mag
                not_v = np.array([1, 0, 0])
                if (v == not_v).all():
                    not_v = np.array([0, 1, 0])
                n1 = np.cross(v, not_v)
                n1 /= norm(n1)
                n2 = np.cross(v, n1)
                t = np.linspace(0, mag, 20)
                theta = np.linspace(0, 2 * np.pi, 20)
                t, theta = np.meshgrid(t, theta)
                
                x = []
                y = []
                z = []
                vv = []
                
                for j in range(0,len(R_)):
                    R = R_[j]
                    X, Y, Z = [p0[i] + v[i] * t + R * np.sin(theta) * n1[i] + R * np.cos(theta) * n2[i] for i in [0, 1, 2]]
                    x.extend(X.flatten())
                    y.extend(Y.flatten())
                    z.extend(Z.flatten())
                    vv.extend(np.ones((len(X.flatten())))*vv_[j])

                # interpolate "data.v" on new grid "inter_mesh"
                V = gd((x,y,z), vv, (XX.flatten(),YY.flatten(),ZZ.flatten()), method='linear')
                V[np.isnan(V)] = 0
                V = np.array(V.reshape(shape))
                conc = conc+V
                print('cyl '+str(i)+' of ' + str(len(self.cyls))+ ' is finished!')

        return conc
        
    def calc_model(number,XX,YY,ZZ,q):
        q.put([(self.calculate_cyls(number,XX,YY,ZZ,))])
        
    def calculate_cyls(i,XX,YY,ZZ): 
        cyl = self.cyls[i]
        R_ = cyl.getDofCoordinates()
        vv_ = cyl.getSolution(1)
        
        p0 = np.array(self.nodes[self.segments[i].x])
        p1 = np.array(self.nodes[self.segments[i].y])
        
        v = p1 - p0
        mag = norm(v)
        v = v / mag
        not_v = np.array([1, 0, 0])
        if (v == not_v).all():
            not_v = np.array([0, 1, 0])
        n1 = np.cross(v, not_v)
        n1 /= norm(n1)
        n2 = np.cross(v, n1)
        t = np.linspace(0, mag, 20)
        theta = np.linspace(0, 2 * np.pi, 20)
        t, theta = np.meshgrid(t, theta)
        
        x = []
        y = []
        z = []
        vv = []
        
        for j in range(0,len(R_)):
            R = R_[j]
            X, Y, Z = [p0[i] + v[i] * t + R * np.sin(theta) * n1[i] + R * np.cos(theta) * n2[i] for i in [0, 1, 2]]
            x.extend(X.flatten())
            y.extend(Y.flatten())
            z.extend(Z.flatten())
            vv.extend(np.ones((len(X.flatten())))*vv_[j])

        # interpolate "data.v" on new grid "inter_mesh"
        V = gd((x,y,z), vv, (XX.flatten(),YY.flatten(),ZZ.flatten()), method='linear')
        V[np.isnan(V)] = 0
        V = np.array(V.reshape(shape))
        return V

    def map_cylinders_solute_parallel(self, XX,YY,ZZ,name = ""):
        """ maps cylinders to soil grid """
        
        shape = np.shape(XX)
        conc = np.zeros((shape))
        q = multiprocessing.Queue()
        for i in range(0,len(self.cyls)):
            if not isinstance(cyl, AirSegment):              
                p = multiprocessing.Process(target=self.calc_model, args=(i,XX,YY,ZZ,q,))
                p.start()
            
        for i in range(0,len(self.cyls)):
            conc = np.add(conc,np.array(q.get()))
        return conc
                  
    def collect_cylinder_solute_data(self):
        """ collects solute data from cylinders  """
        
        l_all = self.segLength()
        a_all = self.radii
        dist = []
        conc = []
        l = []
        a = []
        for i, cyl in enumerate(self.cyls):
            if not isinstance(cyl, AirSegment):  
                x = cyl.getDofCoordinates()
                x_ = x-x[0]
                dist.append(x_)
                conc.append(cyl.getSolution(1))
                l.append(l_all[i]*np.ones((len(cyl.getDofCoordinates()))))
                #print('shape dist', np.shape(dist))

        return  dist, conc, l
    
def plot_transpiration(t, soil_uptake, realised_trans):
    """ plots potential and actual transpiration over time 
    
    depending on discretisation soil_uptake and root_uptake might differ  
    
    @param t                  times [day]
    @param soil_uptake        actual transpiration [cm3/day] of soil 
    @param realised_trans    [cm3/day]
    """
    fig, ax1 = plt.subplots()
    ax1.plot(t, realised_trans, 'k', label = "realised transpiration")  # realised transpiration
    ax1.plot(t, np.array(soil_uptake), 'g', label = "soil uptake")  # actual transpiration  according to soil model
    ax1.set_xlabel("Time [d]")
    ax1.set_ylabel("Transpiration $[cm^3 d^{-1}]$")
    ax2 = ax1.twinx()
    dt = np.diff(t)
    so = np.array(soil_uptake)
    cum_transpiration = np.cumsum(np.multiply(so[:-1], dt))
    ax2.plot(t[1:], cum_transpiration, 'c--', label = "Cumulative soil uptake")  # cumulative transpiration (neumann)
    ax2.set_ylabel("Cumulative soil uptake $[cm^3]$")
    print("Cumulative soil uptake", cum_transpiration[-1], "[cm^3]")
    fig.legend()
    plt.show()


def plot_info(x_, water_collar_cell, water_cyl, collar_sx, min_sx, min_rx, min_rsx, water_uptake, water_domain):
    """ 2x2 plot with additional information """
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
    ax1.set_title("Water amount")
    ax1.plot(x_, np.array(water_collar_cell), label = "water cell")
    ax1.plot(x_, np.array(water_cyl), label = "water cylindric")
    ax1.legend()
    ax1.set_xlabel("Time (days)")
    ax1.set_ylabel("(cm3)")
    ax2.set_title("Pressure")
    ax2.plot(x_, np.array(collar_sx), label = "soil at root collar")
    ax2.plot(x_, np.array(min_sx), label = "min soil")
    ax2.plot(x_, np.array(min_rx), label = "min xylem")
    ax2.plot(x_, np.array(min_rsx), label = "min 1d at root surface")
    ax2.set_ylim([-15000, 0])
    ax2.legend()
    ax2.set_xlabel("Time (days)")
    ax2.set_ylabel("Matric potential (cm)")
    ax3.set_title("Water uptake")
    ax3.plot(x_, -np.array(water_uptake))
    ax3.set_xlabel("Time (days)")
    ax3.set_ylabel("Uptake (cm/day)")
    ax4.set_title("Water in domain")
    ax4.plot(x_, np.array(water_domain))
    ax4.set_xlabel("Time (days)")
    ax4.set_ylabel("cm3")
    plt.show()

# print(len(rs.nodeCTs), len(rs.segments))
# ana2 = pb.SegmentAnalyser(r.rs.nodes, r.rs.segments, r.rs.nodeCTs[1:], r.rs.radii)
# types = np.array(r.rs.subTypes, dtype=np.float64)
# ana2.addData("subType", types)
# ana2.addData("age", r.get_ages())
# pd = vp.segs_to_polydata(ana2, 1., ["radius", "subType", "creationTime", "age"])
# vp.plot_roots(pd, "creationTime")

