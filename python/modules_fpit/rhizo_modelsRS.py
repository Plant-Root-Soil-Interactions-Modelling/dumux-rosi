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
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank(); max_rank = comm.Get_size()
import psutil
from air_modelsPlant import AirSegment
from scipy import sparse
import scipy.sparse.linalg as LA

from scipy.interpolate import PchipInterpolator,  CubicSpline

from scenario_setup import write_file_array, write_file_float
class RhizoMappedSegments(pb.MappedRootSystem):#XylemFluxPython):#
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
                usemoles, seedNum=None,results_dir = './results/' ,limErr1d3dAbs = 1e-11,
                l_ks = "dx_2"):
        """ @param file_name is either a pb.MappedRootSystem, pb.MappedSegments, or a string containing a rsml filename"""
        if not seedNum is None:
            super().__init__(seednum = seedNum)
        else:
            super().__init__()
        self.NC = NC #dof +1
        self.logbase = logbase
        self.mode = mode  # more precise RichardsCylFoam, mode="dumux"
        #print(self.mode)
        self.results_dir = results_dir
        # changes with growing plant (to update)
        # for soil 1ds
        self.nsSoilOld = 0
        self.toAddSoil = np.array([])
        self.repartitionSoilOld = np.array([0 for i in range( max_rank)])
        # for air 1ds
        self.nsAirOld = 0
        self.toAddAir = np.array([])
        self.repartitionAirOld = np.array([0 for i in range( max_rank)])
        
        self.seg_length_ = np.array([])
        self.radii_ = np.array([])
        self.outer_radii_ = np.array([])
        self.rhizoVol_ = np.array([])
        self.seg2cell_old = {}
        self.cell2seg_old = {}
        self.limErr1d3dAbs = limErr1d3dAbs
        self.l_ks = l_ks        
        
        
        self.cyls = []
        self.cylsSoil = []
        self.outer_radii = None
        self.rhizoVol = None
        self.seg_length = np.array([]) 
        self.seg_length_old = np.array([]) 
        self.eidx = np.array([]) #root segments id
        self.cylidx = [] # perirhizal cylinders id
        self.cylSoilidx = [] # perirhizal cylinders id without airSegments
        self.cylSoilidx_all = []
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
        
        # when put some components to 0
        self.canBeNull = np.where(self.soilModel.ICcc == 0.)[0] + 1
        
        # self.setSoilData()# initial soil condition
            
        self.dx2 = []
        self.eidx = np.array([], dtype=np.int64)
        # additional variables
        self.last_dt = 0.

    def reset(self):
        for cyl in self.cyls:
            cyl.reset()
        
    def getDistribution(self, ns:int):      
        dcyl = int(np.floor(ns /  max_rank))
        repartition = np.array([dcyl for i in range( max_rank)])
        residual = ns - max_rank * dcyl 
        repartition[:residual] += 1 # perfect repartition 
        return repartition
    
    def broadcastPlantShape(self):
        #backups
        self.outer_radii_old = self.outer_radii
        self.seg_length_old = self.seg_length
        self.rhizoVol_old = self.rhizoVol
        self.seg2cell_old = self.seg2cell
        self.cell2seg_old = self.cell2seg
        self.airSegs = np.array([])
        if rank == 0:
    
            self.outer_radii = np.array(self.segOuterRadii(type = 0)) 
            aboveGround = np.array([])
            if not (self.cell2seg.get(-1) is None):
                aboveGround = self.cell2seg.get(-1)                
            self.airSegs = np.array(list(set(np.concatenate((aboveGround,
                                                        np.where(np.array(self.organTypes) != 2)[0])) )))
            if len (self.airSegs) > 0:
                self.outer_radii[self.airSegs] = np.array(self.radii)[self.airSegs]*1.1 #dummy outer radii
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
            repartitionSoil = self.getDistribution(len(self.seg_length) - len (self.airSegs))
            self.toAddSoil = repartitionSoil - self.repartitionSoilOld
            
            repartitionAir = self.getDistribution(len (self.airSegs))
            self.toAddAir = repartitionAir - self.repartitionAirOld
            
            #print('parsing plant data_0Soil', rank, self.toAddSoil,repartitionSoil ,
            #      self.repartitionSoilOld,self.nsSoilOld, sum(self.toAddSoil), len(self.radii))
            #print('parsing plant data_0Air', rank, self.toAddAir,repartitionAir ,
            #      self.repartitionAirOld,self.nsAirOld, sum(self.toAddAir), len(self.radii))
            if (min(self.toAddSoil) <0.) or (min(self.toAddAir) <0.) :
                print('min(self.toAddSoil) or min(self.toAddSoil)',min(self.toAddSoil), min(self.toAddAir))
                raise Exception
            self.repartitionSoilOld = repartitionSoil
            self.repartitionAirOld = repartitionAir
            self.organTypes =np.array( self.organTypes)
            assert sum(repartitionSoil) == len(self.seg_length) - len (self.airSegs)
            assert sum(self.toAddSoil) == (len(self.seg_length) - len (self.airSegs) - self.nsSoilOld )
            assert sum(repartitionAir) ==  len (self.airSegs)
            assert sum(self.toAddAir) == ( len (self.airSegs) - self.nsAirOld )
            #for next time step
        
        
        self.airSegs = comm.bcast (self.airSegs, root = 0) 
        self.toAddSoil = comm.bcast(self.toAddSoil, root = 0) 
        self.toAddAir = comm.bcast(self.toAddAir, root = 0) 
        self.seg_length = comm.bcast(self.seg_length, root = 0)  
        self.outer_radii = np.array( comm.bcast(self.outer_radii, root = 0))
        self.radii =np.array( comm.bcast(self.radii, root = 0)    )
        self.seg2cell = comm.bcast(self.seg2cell, root = 0)    
        self.cell2seg = comm.bcast(self.cell2seg, root = 0)         
        self.organTypes =np.array( comm.bcast(self.organTypes, root = 0) ) 
        self.organTypesArray =np.array( comm.bcast(self.organTypes, root = 0)   )
        self.organTypesArray =np.array( self.organTypesArray)
        self.rhizoVol = (self.outer_radii**2 - np.array(self.radii)**2)*np.pi*self.seg_length
        #print('parsing plant data', rank, len(self.seg_length), len(self.radii),
        #      'soil',self.toAddSoil ,  self.nsSoilOld, sum(self.toAddSoil),
        #      'air', self.toAddAir ,  self.nsAirOld, sum(self.toAddAir))
        newLidxSoil = np.array([i for i in range(self.toAddSoil[rank])]) + sum(self.toAddSoil[:rank])+self.nsSoilOld 
        newLidxAir = np.array([i for i in range(self.toAddAir[rank])]) + sum(self.toAddAir[:rank])+self.nsAirOld
        self.rhizoSegsId = np.array([i for i in range(len(self.organTypes)) if i not in self.airSegs])
        newEidxSoil = np.array([], dtype = int)
        if len(newLidxSoil) > 0:
            newEidxSoil = self.rhizoSegsId[newLidxSoil]
        newEidxAir = np.array([], dtype = int)
        if len(newLidxAir) > 0:
            newEidxAir = self.airSegs[newLidxAir]
        self.newEidx = np.concatenate((newEidxSoil, newEidxAir))
        if len(newEidxSoil) > 0:
            assert (np.array(self.organTypes)[newEidxSoil] == 2).all()
        
        #print('parsing plant dataBIS', rank,
        #    np.array( self.organTypesArray),newEidxSoil,newEidxAir,self.newEidx)
        
        
    def update(self):
        """ calls the specific initializer according to @param mode (no parallelisation)  
        @param soil     van genuchten parameters as list        
        @param x        is the solution (or initial condition) of the soil model
        @þaram newEidx  in case of mpi the cylinder indices per process
        """
        
        sizeSoilCell = comm.bcast(self.soilModel.getCellVolumes()/1e6, root = 0) #m3
        try:
            maxlim1d3d = max(self.maxDiff1d3dCW_abs[0]*10,self.limErr1d3dAbs)
        except:
            maxlim1d3d = self.limErr1d3dAbs
        
        print("broadcastPlantShape",rank)
        self.broadcastPlantShape() 
        
        comm.barrier()
        if self.newEidx is None:
            self.newEidx = np.array(range(len(self.eidx), len(self.radii)), np.int64)  # segment indices for process
        self.eidx = np.concatenate((self.eidx,np.array(self.newEidx, dtype = np.int64)), dtype = np.int64)
        
        
        comm.barrier()
        print('checkMassOMoleBalance2 after broadcastPlantShape')
        self.checkMassOMoleBalance2( sourceWat = np.full(len(sizeSoilCell),0.), # cm3/day 
                                     sourceSol = np.full((self.numComp, len(sizeSoilCell)),0.), # mol/day
                                     dt = 0.,        # day    
                                     seg_fluxes = 0.,# [cm3/day]
                                     doWater = True, doSolute = True, doSolid = True,
                                     # useSoilData = True,
                                     diff1d3dCW_abs_lim = maxlim1d3d)
        
        
        self.checkVolumeBalance(finishedUpdate = False)
        self.checkRadii()
        
        self.seg_length_ = self.seg_length[self.eidx]
        self.radii_ = np.array(self.radii)[self.eidx]
        self.outer_radii_ = self.outer_radii[self.eidx]
        self.rhizoVol_ = self.rhizoVol[self.eidx]
        
        cellIds = np.fromiter(self.cell2seg.keys(), dtype=int)
        cellIds =  np.array([x for x in cellIds if x >= 0])#take out air segments
        
        
        print("updateOld", rank)
        for i, cyl in enumerate(self.cyls):
            if not isinstance(cyl, AirSegment) :
                self.updateOld(i, cyl)
                
        if (len(self.seg_length) != len(self.seg_length_old)) or( len(self.cylSoilidx_all) == 0):
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
            #watcont = self.soilModel.getWaterContent()
            # write_file_array('watcontbis',watcont, ',',self.results_dir)
            # if rank == 0:
                # print('XX_leftover[53]',XX_leftover[53], watcont[53])
            # raise Exception
            #raise Exception
            for gId in self.newEidx:#only initialize the new eidx
                self.initialize_(gId,XX_leftover,CC_leftover)
        self.cylidx = self.eidx
        self.cylSoilidx = np.array([gId for lId, gId in enumerate(self.cylidx) if (not isinstance(self.cyls[lId], AirSegment))])
        self.cylsSoil = np.array([cyl for cyl in self.cyls if (not isinstance(cyl, AirSegment))])
        self.cylSoilidx_all = np.concatenate(comm.bcast(comm.gather(np.array(self.cylSoilidx),root = 0),root = 0))# all epreviously existsing segs
        
        self.checkVolumeBalance(finishedUpdate = True)
        print("afterUpdate")
        self.checkMassOMoleBalance2( sourceWat = np.full(len(sizeSoilCell),0.), # cm3/day 
                                     sourceSol = np.full((self.numComp, len(sizeSoilCell)),0.), # mol/day
                                     dt = 0.,        # day    
                                     seg_fluxes = 0.,# [cm3/day]
                                     doWater = True, doSolute = True, doSolid = True,
                                     # useSoilData = True,
                                     diff1d3dCW_abs_lim = max(self.maxDiff1d3dCW_rel[0]*2,self.limErr1d3dAbs) )
        
        self.nsSoilOld = sum(self.repartitionSoilOld )
        self.nsSoilOld = comm.bcast(self.nsSoilOld, root = 0)  
        self.nsAirOld = sum(self.repartitionAirOld )
        self.nsAirOld = comm.bcast(self.nsAirOld, root = 0)  
    
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
        #print('get_XX_leftover_AA', rank, watVol_leftover, volLeftOver, theta_leftOver)
        lowTheta = np.where(theta_leftOver <self.vg_soil.theta_R )
        if len(volLeftOver[lowTheta]) > 0:
            assert max(volLeftOver[lowTheta]) < 1e-13 # low theta, no new segments in the soil voxel
        highTheta = np.where(theta_leftOver >self.vg_soil.theta_S )
        if len(volLeftOver[highTheta]) > 0:
            assert max(volLeftOver[highTheta]) < 1e-13 # low theta, no new segments in the soil voxel
        correctTheta = np.where(((theta_leftOver <= self.vg_soil.theta_S) & (theta_leftOver >= self.vg_soil.theta_R)))
        #print('get_XX_leftover_A', rank, theta_leftOver, correctTheta)
        XX_leftover = np.full(theta_leftOver.shape, np.nan)
        if len(theta_leftOver[correctTheta]) > 0:
            try:
                #print('get_XX_leftover_B', rank,np.array([tlo for tlo in theta_leftOver[correctTheta]]))
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
            vol_total = self.soilModel.getCellVolumes() #cm3
            if rank == 0:
                vol_total = vol_total[cellId] # cm3
                
            vol_total = comm.bcast(vol_total, root = 0) 
            
            idSegs = np.array(self.cell2seg[cellId])#all segments in the cell
            idCylsAll = np.array([ids for ids in idSegs if self.organTypes[ids]==2 ]) 
            if len(idCylsAll)>0:
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
                            print("\n",rank,cellId,"checkRadii error",  vol_total,
                                  'vol_rootRhizo',vol_rootRhizo,sum(vol_rootRhizo))
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
                vol_rhizo_ =self.getVolumesCyl(idCyls,doSum = False)# comm.bcast(self._flat0(comm.gather(vol_rhizo_, root=0)), root=0)
                vol_rhizo = sum(vol_rhizo_)#cm3 sum(vol_rhizo_)
                    
                
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
                                     diff1d3dCW_abs_lim = None):#would need to do it for each cell, not overall?
        verbose = False
        if verbose:
            print("checkMassOMoleBalance2", rank)
        if diff1d3dCW_abs_lim is None:
            diff1d3dCW_abs_lim = self.limErr1d3dAbs
        
        cellIds = np.fromiter(self.cell2seg.keys(), dtype=int)
        cellIds =  np.array([x for x in cellIds if x >= 0])#take out air segments
        if verbose:
            print("cellIds for checkMassOMoleBalance2",rank, cellIds)
        
        cellIds_All = np.array(comm.bcast(comm.gather(cellIds, root=0), root=0))

        try:
            assert (cellIds_All[1:] == cellIds_All[:-1]).all()
        except:
            print('different cell ids', cellIds_All)
            raise Exception
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
        
        
        if useSoilData:
            raise Exception
            wat_total_ = self.soilWatVol_old * 1e6# [cm3].
        else:
            wat_total_ = self.soilModel.getWaterVolumes()# [cm3].
            
        wat_total_ = comm.bcast(wat_total_, root=0)
        
        self.all1d3dDiff =  np.full(len(wat_total_), 0.)

        if useSoilData:
            raise Exception
            mol_total_ = [self.soilContent_old[idComp] for  idComp in range(NC)]# solute content [mol]
        else:
            mol_total_ = [self.soilModel.getContent(idComp, idComp <= 2)  for  idComp in range(1,NC+1) ]# solute content [mol].
        
        mol_total_ = comm.bcast(mol_total_, root=0)
        
        for cellId in cellIds:
            verbose = False
            if verbose:
                print('check cell', cellId,rank)
            comm.barrier()
            wat_total = wat_total_[cellId] # [cm3].
            
            idCylsAll, idCyls = self.getIdCyllMPI(cellId)
            if verbose:
                print("rank",rank,"cell",cellId,  "\tidCylsAll",idCylsAll)
            if len(idCylsAll)> 0:
                localIdCyls = np.array([ np.where(self.eidx == i)[0] for i in idCyls]).flatten()#np.where(idCyls in self.eidx)[0]
                wat_rhizo_ = self.getWaterVolumesCyl(idCyls, doSum = False) #cm3
                wat_rhizo = sum(wat_rhizo_)
                    
                
                
                if rank == 0:
                    if verbose:
                        print("in water",rank, cellId, 'wat_total',wat_total,'wat_rhizo',  wat_rhizo, sourceWat[cellId] * dt , idCyls , self.eidx,idCyls,localIdCyls )    
                    #print(cellId, 0,wat_total , sourceWat[cellId] * dt , wat_rhizo, wat_rhizo_)
                    diff1d3dCW_abs = abs(wat_total+ sourceWat[cellId] * dt - wat_rhizo)
                    self.all1d3dDiff[cellId] = diff1d3dCW_abs
                    diff1d3dCW_rel = abs(diff1d3dCW_abs/(wat_total+ sourceWat[cellId] * dt) *100)
                    self.sumDiff1d3dCW_abs[0] += diff1d3dCW_abs
                    self.sumDiff1d3dCW_rel[0] += diff1d3dCW_rel
                    self.maxDiff1d3dCW_abs[0] = max(self.maxDiff1d3dCW_abs[0], diff1d3dCW_abs)
                    self.maxDiff1d3dCW_rel[0] = max(self.maxDiff1d3dCW_rel[0], diff1d3dCW_rel)
                    if diff1d3dCW_abs > diff1d3dCW_abs_lim:#1e-10 :
                        print("checkMassOMoleBalance2error")
                        print(cellId, 0, 'wat_total',wat_total,'wat_rhizo', wat_rhizo, sourceWat[cellId], dt)
                        print("wat_rhizo_",wat_rhizo_)
                        print("Diff1d3dW",diff1d3dCW_abs, diff1d3dCW_rel, self.sumDiff1d3dCW_abs[0],
                            self.sumDiff1d3dCW_rel[0], self.maxDiff1d3dCW_abs[0], self.maxDiff1d3dCW_rel[0],
                            diff1d3dCW_abs_lim)
                        print(np.array(self.organTypes)[idCylsAll])
                        print(np.array(self.radii)[idCylsAll],np.array(self.outer_radii)[idCylsAll])
                        # raise Exception #turn off for now to make evaluation easier 
                    
                for idComp in range(1, NC +1): 
                    #comm.barrier()
                    isDissolved = (idComp <= self.numFluidComp)
                    if verbose:
                        print("check comp",  idComp,cellId,rank, len(idCylsAll))
                    #comm.barrier()
                        
                    mol_total = mol_total_[idComp -1][cellId] # solute content [mol]
                    
                    #comm.barrier()
                    #print('getting mol total', rank, idComp,cellId)
                    mol_total = comm.bcast(mol_total, root=0)
                    #print('got mol total', rank, idComp,cellId)
                    #comm.barrier()
                    
                    if idComp not in self.canBeNull:
                        assert mol_total > 0.
                    else:
                        assert mol_total >= 0.
                                        
                    #mol_rhizo_ = np.array([sum(self.cyls[i].getContentCyl( idComp, isDissolved, self.seg_length_[i] )) for i in localIdCyls ]) #cm3
                    #print('getting mol rhizo', rank, idComp,cellId)
                    mol_rhizo_ = self.getContentCyl(idCyls, idComp, doSum = False)#sum(comm.bcast(self._flat0(comm.gather(mol_rhizo_, root=0)), root=0))
                    mol_rhizo = sum(mol_rhizo_)
                    #print('got mol rhizo', rank, idComp,cellId)
                        
                            
                    if rank == 0:
                        try:
                            sourceSol_ = sourceSol[idComp-1][cellId]             
                        except:
                            print('error could not get sourcesol_')
                            print(sourceSol, idComp, cellId)
                            raise Exception
                        if (mol_total  + sourceSol_ * dt  == 0) :
                            diff1d3dCW_abs = abs(mol_total  + sourceSol_ * dt)
                            diff1d3dCW_rel = np.nan
                            if verbose:
                                print("content", cellId, idComp, mol_total, mol_rhizo,"sourceSol_ * dt",sourceSol_ * dt, 
                                        diff1d3dCW_abs, diff1d3dCW_rel)
                            #print(cellId, idComp,mol_total , sourceSol[idComp-1][cellId] * dt , mol_rhizo, mol_rhizo_)
                            if ((mol_rhizo > 0)):
                                print('error for component', rank, cellId)
                                print(idComp, mol_total, mol_rhizo)
                                raise Exception
                        else:
                            #print(cellId, idComp,mol_total , sourceSol[idComp-1][cellId] * dt , mol_rhizo, mol_rhizo_)
                            #print(cellId, idComp,mol_total , sourceSol[idComp-1][cellId] * dt)
                            diff1d3dCW_abs = abs(mol_total + sourceSol_ * dt - mol_rhizo)
                            diff1d3dCW_rel = abs(diff1d3dCW_abs/(mol_total +sourceSol_ * dt) *100)
                            if verbose:
                                print("content", cellId, idComp, mol_total, 'mol_rhizo',mol_rhizo,mol_rhizo_,"sourceSol_ * dt",sourceSol_ * dt, 
                                        diff1d3dCW_abs, diff1d3dCW_rel)
                            if diff1d3dCW_abs > diff1d3dCW_abs_lim:
                                print('error for component', rank, cellId)
                                print(cellId, idComp,mol_total , sourceSol_ , dt)
                                print(mol_rhizo, abs((mol_total + sourceSol_ * dt - mol_rhizo)/(mol_total + sourceSol_* dt)*100))                                
                                # raise Exception #turn off for now to make evaluation easier 
                                
                        
                        self.sumDiff1d3dCW_abs[idComp] += diff1d3dCW_abs
                        self.sumDiff1d3dCW_rel[idComp] += diff1d3dCW_rel
                        self.maxDiff1d3dCW_abs[idComp] = max(self.maxDiff1d3dCW_abs[idComp], diff1d3dCW_abs)
                        self.maxDiff1d3dCW_rel[idComp] = max(self.maxDiff1d3dCW_rel[idComp], diff1d3dCW_rel)
                        if verbose:
                            print("check compBis", idComp,cellId, rank, mol_rhizo,mol_total, sourceSol_ * dt,"error",diff1d3dCW_abs, diff1d3dCW_rel)
                    #comm.barrier()
            #print('finished check cell',cellId,rank)
            #comm.barrier()
        #if verbose:
        #("finished checkmass2", rank)
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
        # comm.barrier()
        # if rank == 0:
            # print('\n\t\tgetXcyl', rank ,'X_rhizo',X_rhizo, 'X_rhizo-0.00071296' ,X_rhizo-0.00071296, 
                    # 'sum', sum(X_rhizo) ,"\n")
        # comm.barrier()
        # if rank == 1:
            # print('\t\tgetXcyl', rank ,'X_rhizo',X_rhizo, 'X_rhizo-0.00071296' ,X_rhizo-0.00071296, 
                    # 'sum', sum(X_rhizo) ,"\n")
        # comm.barrier()
        if reOrder:
            idCyll = np.array( comm.bcast(self._flat0(comm.gather(self.eidx, root=0)), root=0))
            try:
                assert X_rhizo.shape == idCyll.shape
            except:
                print(X_rhizo.shape, idCyll.shape)
                raise Exception
            # for some reason, I need to create another variable (X_rhizo[idCyll] = X_rhizo fails)
            X_rhizo2 = np.full(len(X_rhizo), np.nan)
            X_rhizo2[idCyll] = X_rhizo
        else:
            X_rhizo2 = X_rhizo
        if doSum:
            X_rhizo2 = sum(X_rhizo)
        return X_rhizo2
        
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
        #print('lens', len(self.eidx), len(wat_rhizo))
        #print('getWaterVolumesCyl',wat_rhizo,idCyls, np.array(self.organTypes)[np.array(self.eidx)]  )
        return self.getXcyl(wat_rhizo, doSum, reOrder)
        
    def getVolumesCyl(self,idCyls=None, doSum = True, reOrder = False):#cm3
        localIdCyls =   self.getLocalIdCyls(idCyls)                     
        V_rhizo = np.array([sum(self.cyls[i].getCellSurfacesCyl()*self.seg_length_[i]) for i in localIdCyls ]) #cm3
        #print('getVolumesCyl',rank, V_rhizo )
        return self.getXcyl(V_rhizo, doSum, reOrder)
        
       
    def getKrw(self,idCyls=None):#[-] unitless
        doSum = False
        reOrder = True
        localIdCyls =   self.getLocalIdCyls(idCyls)                     
        krws = np.array([self.cyls[i].getKrw()[0] for i in localIdCyls ]) 
        return self.getXcyl(krws, doSum, reOrder)
        
    def getDeltaR(self,idCyls=None):#[-] unitless
        doSum = False
        reOrder = True
        localIdCyls =   self.getLocalIdCyls(idCyls)      
        cc0 = np.array([self.cyls[i].getCellCenters()[0] for i in localIdCyls ]) 
        p0 = np.array([ self.cyls[i].getPoints()[0]  for i in localIdCyls ]) 
        deltaR = cc0 - p0
        return self.getXcyl(deltaR, doSum, reOrder)
    
            
    def getC_rhizo(self,soilShape, idComp = 1, konz = True): # if konz:  mol/m3 wat or mol/m3 scv, else: mol
        """ return component concentration or content, only in voxel with 1D-domains"""
        #print('getC_rhizo', rank, idComp, konz)
        isDissolved = (idComp <= self.numFluidComp)
        
        contentOrKonz = np.full(soilShape, 0.)
        
            
        cellIds = np.fromiter(self.cell2seg.keys(), dtype=int)
        cellIds =  np.array([x for x in cellIds if x >= 0])#take out air segments
        for idCell in cellIds:
            idSegs = self.cell2seg[idCell]#all segments in the cell
            #print('getC', rank , idCell)
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
                
        return contentOrKonz
        
    def checkSoilBC(self, outer_R_bc_wat, #cm3 wat
                    outer_R_bc_sol):#mol
        raise Exception
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
        raise Exception
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
        raise Exception
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
        #print('setSoilData',rank, self.soilModel.getWaterContent())
        #raise Exception
        
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
    
    def updateConcentration(self,totContent, gradient, cellCenters, volumes,isWater):
        """ give the new theta or solute mole fraction to get specific total content and gradient
            volumes: either bulk soil or water volumes
        """
        matrix_size = self.NC -1
        sub_diag_values = -1.
        main_diag_values = 1.
        matrix = np.diag(np.full(matrix_size-1,sub_diag_values), k=1) + np.diag(np.full(matrix_size,main_diag_values), k=0) 
        matrix[-1,] = volumes
        aB = - gradient * (cellCenters[1:] - cellCenters[:-1])
        aB = np.append(aB,totContent)
        SPmatrix = sparse.csc_matrix(sparse.coo_matrix(matrix))
        val_new = LA.spsolve(SPmatrix, aB, use_umfpack = True) #either that or mol fraction
        try:
            assert min(val_new) >= 0. # need to also check that theta >= theta_r and <= theta_s
            assert abs(sum(val_new *volumes ) - totContent) < 1e-13 #check tot content ok
            assert (abs(((val_new[1:] - val_new[:-1] )/ (cellCenters[1:] - cellCenters[:-1])) - gradient) < 1e-13).all() #check gradient ok
        except:
            print(rank, 'updateConcentration error', val_new,totContent, gradient, cellCenters, volumes )
            raise Exception
            
        if isWater and ((min(val_new) <  self.vg_soil.theta_R) or (max(val_new) >  self.vg_soil.theta_S)):
            print(rank, 'updateConcentration: issue with theta concentraiton, adapt gradient', val_new )
            tempTheta = val_new
            if (max(val_new) >  self.vg_soil.theta_S):
                idtemp = np.where(val_new >  self.vg_soil.theta_S)
                extraWater = sum((val_new[idtemp] - self.vg_soil.theta_S) * volumes[idtemp])
                val_new[idtemp] = self.vg_soil.theta_S
                while  extraWater > 0:
                    print(rank, 'updateConcentration: take out extra water', extraWater )
                    canAdd = (self.vg_soil.theta_S - val_new)* volumes 
                    assert sum(canAdd) >= extraWater 
                    idtemp = np.where(canAdd > 0 )
                    val_new[idtemp] = (val_new[idtemp]*volumes[idtemp]+extraWater/len(idtemp))/volumes[idtemp]
                    extraWater -=  extraWater/len(idtemp)
                print(rank, 'updateConcentration: finished taking out extra water', extraWater , val_new)
                
            if (min(val_new) <  self.vg_soil.theta_R):
                idtemp = np.where(val_new <  self.vg_soil.theta_R)
                missingWater = sum((self.vg_soil.theta_R - val_new[idtemp]) * volumes[idtemp])
                val_new[idtemp] = self.vg_soil.theta_R
                while  missingWater > 0:
                    print(rank, 'updateConcentration: add missing water', missingWater )
                    canTake = (val_new - self.vg_soil.theta_R )* volumes 
                    assert sum(canTake) >= missingWater 
                    idtemp = np.where(canTake > 0 )
                    val_new[idtemp] = (val_new[idtemp]*volumes[idtemp] - missingWater/len(idtemp))/volumes[idtemp]
                    missingWater -=  missingWater/len(idtemp)
                print(rank, 'updateConcentration: finished adding missing water', missingWater , val_new)
                
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

            
    def updateOld(self, lId,  cyl):
        """ update distribution of cylinder if it s volume has changed
        """
        gId = self.eidx[lId]
        cellId = self.seg2cell[gId]
        verbose = False # (cellId == 437) or (cellId == 257)
        segsId =  np.array(self.cell2seg[cellId])
        segsId = np.array([ids for ids in segsId if self.organTypes[ids]==2 ])# only root organs
        if verbose:
            print(rank,'updateold', np.array(self.cylSoilidx),self.cylSoilidx_all)
        try:
            hasNewSegs =(np.array([ bool(rootId not in self.cylSoilidx_all) for rootId in segsId]).any())
        except:
            print(rank, 'error with ahsnewsegs', lId, gId,  np.array(self.cell2seg[cellId]), 
                  np.array(self.organTypes)[np.array(self.cell2seg[cellId])])
            print(rank, segsId,self.cylSoilidx, self.cylSoilidx_all)
            raise Exception
        oldPoints = np.array(cyl.getPoints()).flatten() # cm
        lb = self.logbase            
        a_in = self.radii[gId] # cm
        a_out = self.outer_radii[gId]  # cm
        
        if ((self.seg2cell[gId] > 0) and (self.organTypes[gId] == 2)):
            points = np.logspace(np.log(a_in) / np.log(lb), np.log(a_out) / np.log(lb), self.NC, base = lb) # cm
        else:
            points = np.array([a_in,a_out ])
        if (len(points) != len(oldPoints)) or (max(abs(points - oldPoints)) > 1e-13):#new vertex distribution
            if verbose:
                print(rank, 'updateOld', cellId, lId, gId, 'hasNewSegs',hasNewSegs )
            ## old shape:
            centersOld = np.array(cyl.getCellCenters()).flatten();assert len(centersOld)==9      # cm 
            volOld = self.getVolumes(oldPoints, self.seg_length_old[gId] )# cm3
            ## new shape: 
            volNew = self.getVolumes(points, self.seg_length[gId] ) # cm^3
            centersNew = (points[1:] + points[:-1])/2  # cm
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
                print('radii', self.radii[gId], self.outer_radii[gId],self.outer_radii_old[gId])
                raise Exception
            if verbose:
                print('\t',gId,'volNew', volNew,'volOld',volOld)
                print('\t',gId,"points",oldPoints, points)
                print('\t',gId,'lengths', self.seg_length_old[gId],self.seg_length[gId])
                print('\t',gId,'radii', self.radii[gId], self.outer_radii[gId],self.outer_radii_old[gId])
            ##  water:
            theta_old = cyl.getWaterContent() # cm3/cm3
            gradientOld = (theta_old[1:] - theta_old[:-1])/(centersOld[1:] - centersOld[:-1])           
            gradientNew = self.interAndExtraPolation(points[1:-1],oldPoints[1:-1], gradientOld)
            
            wOld = sum(theta_old*volOld)
            changeRatioW = max(min(changeRatio, sum(self.vg_soil.theta_S*volNew)/wOld),sum(self.vg_soil.theta_R*volNew)/wOld)
            try:
                assert ((changeRatioW <= 1.) and (changeRatioW > 0.))
            except:
                print('volNew', volNew,'volOld',volOld)
                print("points",oldPoints, points)
                print('lengths', self.seg_length_old[gId],self.seg_length[gId])
                print('radii', self.radii[gId], self.outer_radii[gId],self.outer_radii_old[gId])
                print('theta', theta_old)
                raise Exception
            theta_new = self.updateConcentration(totContent = wOld*changeRatioW, gradient =gradientNew, cellCenters = centersNew, volumes = volNew,isWater = True)
            newHead = np.array([vg.pressure_head(nt, self.vg_soil) for nt in theta_new])# cm

            if verbose:            
                print('\t',gId,"x_old",theta_old* volOld,volOld )
                print('\t',gId,"xnew", theta_new* volNew,volNew )
                print('\t',gId,"newHead", newHead )
                print('\t',gId,"points",oldPoints, points,centersNew)
                print('\t',gId,"gradient",gradientOld, gradientNew)
                print('\t',gId,"theta",theta_new, theta_old)
                print('\t',gId,"vg param", self.vg_soil.theta_R, self.vg_soil.theta_S)   
                print('\t',gId,"changeRatio",changeRatio, changeRatioW)    
            
                
            ## new contents:  
            molFrOld =np.array( [np.array(cyl.getSolution(nC+1)) for nC in range(self.numComp)])   #mol/mol 
            volWatNew = theta_new *volNew
            molFrNew = []
            for nComp in range(1, self.numComp +1):
                isDissolved = (nComp <= self.numFluidComp)
                if isDissolved: # mol phase = [cm3 phase] * [m3/cm3] * [mol phase /m3 phase] 
                    molarPhaseOld = theta_old*volOld/1e6 * self.phaseDensity(nComp ) 
                    molarPhaseNew = volWatNew/1e6 * self.phaseDensity(nComp ) # 
                else:
                    molarPhaseOld = volOld/1e6 * self.phaseDensity(nComp )
                    molarPhaseNew = volNew/1e6 * self.phaseDensity(nComp ) 
                gradientOld = (molFrOld[nComp -1][1:] - molFrOld[nComp -1][:-1])/(centersOld[1:] - centersOld[:-1])   
                gradientNew = self.interAndExtraPolation(points[1:-1],oldPoints[1:-1], gradientOld)
                cOld = sum(molFrOld[nComp -1] *   molarPhaseOld  )    
                try:
                    assert abs(cOld - sum(cyl.getContentCyl(nComp, isDissolved, self.seg_length_old[gId] )))< 1e-13
                except:
                    print('contentoldA',cyl.getContentCyl(nComp, isDissolved, self.seg_length_old[gId] ), 
                            sum(cyl.getContentCyl(nComp, isDissolved, self.seg_length_old[gId] )))
                    print('contentoldB',molFrOld[nComp -1] *   molarPhaseOld, cOld,'diff', cOld - sum(cyl.getContentCyl(nComp, isDissolved, self.seg_length_old[gId] )))
                    print('watcontA',cyl.getWaterContent(), cyl.getCellSurfacesCyl() / 1e4 * self.seg_length_old[gId] / 100)
                    print('watcontB',theta_old,volOld/1e6, self.phaseDensity(nComp ) , cyl.phaseDensity(nComp ))
                    print('molFrOld',molFrOld[nComp -1])
                    print('molarPhaseOld', molarPhaseOld,np.multiply(cyl.getCellSurfacesCyl() / 1e4 * self.seg_length_old[gId] / 100 , 
                                            cyl.getWaterContent()) *cyl.molarDensityWat_m3)
                    raise Exception
                    
                molFrNew.append(self.updateConcentration(totContent = cOld*changeRatio, gradient =gradientNew, 
                                cellCenters = centersNew, volumes = molarPhaseNew,isWater = False))  
                if verbose:            
                    print('\t',gId,'nComp', nComp, "cOld",cOld,molFrOld[nComp -1],molarPhaseOld )
                    print('\t',gId,"molFrNew", molFrNew[nComp -1], 
                                "cNew", molFrNew[nComp -1]* molarPhaseNew, 
                                sum(molFrNew[nComp -1]* molarPhaseNew) )
                    print('\t',gId,"gradient", gradientOld, gradientNew,'ratio', 
                    sum(molFrNew[nComp -1]* molarPhaseNew)/sum(molFrOld[nComp -1] *   molarPhaseOld)) 
                    print('\t',gId,"error",nComp,  abs(sum(molFrNew[nComp -1]* molarPhaseNew)/sum(molFrOld[nComp -1] *   molarPhaseOld) -changeRatio),
                            abs(sum(molFrNew[nComp -1]* molarPhaseNew)/sum(molFrOld[nComp -1] *   molarPhaseOld) -changeRatio) < 1e-13)
                assert abs(sum(molFrNew[nComp -1]* molarPhaseNew)/sum(molFrOld[nComp -1] *   molarPhaseOld) -changeRatio) < 1e-13
            self.cyls[lId] = self.initialize_dumux_nc_( gId, 
                                                        x = newHead,# cm
                                                        cAll = molFrNew, # mol/mol water or mol/mol scv
                                                        Cells = centersNew) # cm
            
            
            
    def getVolumes(self,vertices, length):
        return   np.pi * (vertices[1:] ** 2 - vertices[:-1]**2) * length
    
    def get_Vol_leftoverI(self, idCell):# m3
        comm.Barrier()
        
        idCylsAll = []
        idSegs = []
        if idCell in self.cell2seg:
            idSegs = self.cell2seg[idCell]#all segments in the cell
            idSegs=np.array( [ids for ids in idSegs if (self.organTypes[ids] == 2)])#only root segments
            idCylsAll, idCyls =   self.getIdCyllMPI(idCell, True)  
        
        newVol = 0
        if len(idSegs) > len(idCylsAll):
            newSegs = np.array([ids for ids in idSegs if ids not in idCylsAll])
            newVol = sum(np.pi* (np.array(self.outer_radii)[newSegs]**2 - np.array(self.radii)[newSegs]**2 )*self.seg_length[newSegs]/1e6)
        return newVol
        
    def get_watVol_leftoverI(self, idCell):# m3
        comm.Barrier()
        wat_total = self.soilModel.getWaterContent() # m3/m3 #self.soilWatVol_old
        sizeSoilCell = self.soilModel.getCellVolumes()/1e6 #m3
        if rank == 0:
            wat_total = wat_total[idCell] * sizeSoilCell[idCell] # [m3] water content before bulk flow
        wat_total = comm.bcast(wat_total,root = 0)
        wat_rhizo = 0
        idCyls = []
        idCylsAll = []
        if idCell in self.cell2seg:
            idSegs = self.cell2seg[idCell]#all segments in the cell
            idSegs=np.array( [ids for ids in idSegs if (self.organTypes[ids] == 2)])#only root segments
            idCylsAll, idCyls =   self.getIdCyllMPI(idCell, True)   
            if len(idCylsAll) > 0: # we have old segments
                wat_rhizo= self.getWaterVolumesCyl(idCyls)/1e6 #m3
                #wat_rhizo = sum(np.array([sum(self.cyls[int(np.where(self.eidx == gId)[0])].getWaterVolumesCyl(self.seg_length[gId]))/1e6 for lId, gId in enumerate(idCyls) ])) #[m3]
        # print('get_watVol_leftoverI_B',rank, idCell,wat_total,wat_rhizo,idCylsAll,idCyls)
        watVol_res = wat_total - wat_rhizo #m3
        if watVol_res < 0.:
            if ((watVol_res > -1e-13) and (len(idSegs) == len(idCylsAll))): # rounding error probably
                watVol_res = 0.
            else:
                print("getXX_leftoverI")
                print(idCell, wat_total,wat_rhizo, wat_total - wat_rhizo)
                print( len(idSegs), len(idCyls),  len(idCylsAll))# how many old and new cylinders?
                raise Exception
        comm.Barrier()        
        return watVol_res
        
    
    def getC_content_leftoverI(self, idCell, idComp):# mol
        comm.Barrier()
        isDissolved = (idComp <= self.numFluidComp)
        mol_total = self.soilModel.getContent(idComp, isDissolved)#mol #self.soilContent_old[idComp]
        if rank == 0:
            mol_total = mol_total[idCell] # mol
        mol_total = comm.bcast(mol_total,root = 0)
        mol_rhizo = 0
        mol_rhizo_ = 0
        idCyls = []
        idCylsAll = []
        if idCell in self.cell2seg:
            idSegs = self.cell2seg[idCell]#all segments in the cell
            idSegs=np.array( [ids for ids in idSegs if (self.organTypes[ids] == 2)])#only root segments
            idCylsAll, idCyls =  self.getIdCyllMPI(idCell, True)   
            
            if len(idCylsAll) > 0:#we have old segments
                mol_rhizo_= self.getContentCyl(idCyls, idComp=idComp, doSum = False)#adapted for MPI
                mol_rhizo = sum(mol_rhizo_)
                #mol_rhizo = sum(np.array([sum(self.cyls[int(np.where(self.eidx == gId)[0])].getContentCyl(idComp, isDissolved, self.seg_length[gId])) for lId, gId in enumerate(idCyls) ])) #[m3]

        res_CC = mol_total - mol_rhizo
        #if rank==0:
        #    print('getC_content_leftoverI',rank, idCell,idComp,res_CC,'mol_total',mol_total, 'mol_rhizo',mol_rhizo,mol_rhizo_,idCylsAll,idCyls)
        if res_CC < 0.:
            if ((res_CC > -1e-13) and (len(idSegs) == len(idCylsAll))): # rounding error probably 
                res_CC = 0.
            else:
                print("getC_content_leftoverI", idCell)
                print("res_CC = ",res_CC," < 0, idComp:",idComp,'mol_total',mol_total ,'mol_rhizo', mol_rhizo,"res_CC",res_CC ) 
                print("mol_rhizo_", mol_rhizo_, "idCylsAll",idCylsAll)
                print( len(idSegs), len(idCyls))# how many old and new cylinders?
                raise Exception
        comm.Barrier()
        return res_CC

                
    def initialize_(self,gId,x,cc ):
        if ((self.seg2cell[gId] >= 0) and (self.organTypes[gId] == 2)):
            self.cyls.append(self.initialize_dumux_nc_(gId, x[self.seg2cell[gId]], cc[self.seg2cell[gId]]))
        else:
            a_in = self.radii[gId]#get Perimeter instead? not used for now anyway
            a_out = self.outer_radii[gId]
            self.cyls.append(AirSegment(a_in, a_out)) #psi_air computed by the photosynthesis module.
    

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
        #print('to wrap', rank,gId, a_in,a_out )
        if a_in < a_out:
            cyl = RichardsNoMPIWrapper(Richards10CCylFoam(), self.useMoles)  # only works for RichardsCylFoam compiled without MPI
            cyl.initialize(verbose = False)
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
            
            if self.l_ks == "dx_2":
                cyl.setParameter("Soil.BC.dzScaling", "1")
            elif self.l_ks == "dx":
                cyl.setParameter("Soil.BC.dzScaling", "2")
            else:
                raise Exception
                
            cyl.setParameter("Problem.verbose", "0")
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

            cyl.setParameter("Soil.k_sorp", str(self.soilModel.k_sorp)) # mol / cm3
            cyl.setParameter("Soil.f_sorp", str(self.soilModel.f_sorp)) #[-]
            cyl.setParameter("Soil.CSSmax", str(self.soilModel.CSSmax)) #[mol/cm3 scv]
            cyl.setParameter("Soil.alpha", str(self.soilModel.alpha)) #[1/d]


            cyl.setParameter("1.Component.LiquidDiffusionCoefficient", str(self.soilModel.Ds)) #m^2/s

            cyl.setParameter("2.Component.LiquidDiffusionCoefficient", str(self.soilModel.Dl)) #m^2/s

            
            for j in range( 1, self.numComp+1):
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
                        
            cyl.setParameter("Newton.MaxRelativeShift",str(self.soilModel.MaxRelativeShift))
            cyl.initializeProblem()
            cyl.setCriticalPressure(self.soilModel.wilting_point)  # cm pressure head
            cyl.bulkDensity_m3 = self.soilModel.bulkDensity_m3
            cyl.solidDensity =self.soilModel.solidDensity 
            cyl.solidMolarMass =self.soilModel.solidMolarMass
            cyl.solidMolDensity =self.soilModel.solidMolDensity         
            cyl.ddt = 1.e-5
            cyl.gId = gId    
            ThetaCyl = cyl.getWaterContent()
            assert (ThetaCyl >= self.vg_soil.theta_R).all()
            assert (ThetaCyl <= self.vg_soil.theta_S).all()
            #print(cyl.getSolutionHead(),x, gId, self.seg2cell[gId])#cyl.getSolutionHead_(),
            #raise Exception
            assert len(cyl.getSolutionHead()) == (self.NC - 1)
            return cyl
        else:
            print("RhizoMappedSegments.initialize_dumux_: Warning, segment {:g} might not be in domain, radii [{:g}, {:g}] cm".format(i, a_in, a_out))
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
        isAirSegs = self.airSegs 
        #np.array(list(set(np.concatenate((self.cell2seg.get(-1),
        #                                              np.where(np.array(self.organTypes) != 2)[0])) )))
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

    def solve(self, dt, n_iter, *argv):#dt, seg_rx, proposed_outer_fluxes,kex,proposed_outer_sol_fluxes
        """ set bc for cyl and solves it """
        self.last_dt = dt
        verbose = False
        #inner = bot = plant
        #outer = top = soil
        proposed_inner_fluxes = argv[0]
        proposed_outer_fluxes = argv[1]

        self.seg_fluxes_limited = np.full(len(self.cyls), np.nan) # store for post processing
        
        for lId, cyl in enumerate(self.cyls):  # run cylindrical models
            gId = self.eidx[lId]  # for one process gId == lId
            if verbose:
                print(rank, "cyl no ",lId+1,"/",len(self.cyls),'gId')
            if isinstance(cyl, AirSegment):  
                cyl.setInnerFluxCyl(proposed_inner_fluxes[gId])
                self.seg_fluxes_limited[lId] = proposed_inner_fluxes[gId]
            else:
                #if n_iter > 0:
                #    cyl.reset()
                    
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
                    botVal = argv[2][gId]
                else:
                    botVal = argv[2]
                if not isinstance(argv[3],float): 
                    topVal = argv[3][gId]
                else:
                    topVal = argv[3]
                if not isinstance(argv[4],float): 
                    botVal_mucil = argv[4][gId]
                else:
                    botVal_mucil = argv[4]
                if not isinstance(argv[5],float): 
                    topVal_mucil = argv[5][gId]
                else:
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
                
                QflowIn_temp = QflowIn
                QflowOut_temp = QflowOut
                maxDt_temp = 2500/(24*3600)
                maxRelShift = self.soilModel.MaxRelativeShift
                try:
                    if verbose:
                        print(rank, lId,gId, 'start solve','buWBefore',buWBefore,'Qin',QflowIn,
                              proposed_inner_fluxes[gId],qIn, 'QflowOut',QflowOut,'dt',dt )
                                        
                    redoSolve = True
                    n_iter_solve = 0
                    while redoSolve:
                        try:
                            cyl.ddt = 1.e-5 #do I need to reset it each time?
                            cyl.solve(dt, maxDt = maxDt_temp)
                                             
                            # newton parameters are re-read at each 'solve()' calls
                            neumanns = cyl.getAllNeumann(0)    
                            QflowIn_limited = neumanns[0] * (2 * np.pi * self.radii[gId] * l)
                            QflowOut_limited = neumanns[self.NC -2 ] * (2 * np.pi *np.array( self.outer_radii)[gId] * l)
                            
                            buWAfter_ =  cyl.getWaterVolumesCyl(l)
                            buWAfter = sum(buWAfter_ )
                            buCCAfter_ = cyl.getContentCyl(1, True, l)
                            buTotCAfter = self.getTotCContent(cyl,l)
                            diffWproposed = buWAfter - ( buWBefore + (QflowIn + QflowOut) * dt)
                            diffWtemp = buWAfter - ( buWBefore + (QflowIn_temp + QflowOut_temp) * dt)
                            diffWlimited = buWAfter - ( buWBefore + (QflowIn_limited + QflowOut_limited) * dt)
                            if verbose:
                                print(rank, lId,gId,n_iter_solve,dt, 'end solve QflowIn',QflowIn,
                                      'QflowIn_limited',QflowIn_limited,'diff',QflowIn - QflowIn_limited,
                                     'QflowOut',QflowOut, 'QflowOut_limited',QflowOut_limited,
                                      'diff',QflowOut - QflowOut_limited, 'buWAfter',buWAfter,
                                      'buWBefore',buWBefore,'diffWproposed',diffWproposed,
                                      'diffWtemp',diffWtemp,'diffWlimited',diffWlimited,'neumanns',neumanns,
                                      'with qflowIn',qIn,' and qflowout',qOut,cyl.getCellCenters())


                            cyl.setParameter("Newton.MaxRelativeShift", str(self.soilModel.MaxRelativeShift))
                            redoSolve = False
                        except:     
                            if n_iter_solve < 5:
                                print(rank, lId,gId,' with qflowIn',qIn,' or qflowout',qOut,'maxDt',maxDt_temp,'dt',dt,
                                      '. Increase NewtonMaxRelativeShift from',maxRelShift,'to',maxRelShift*10.)
                                maxRelShift *= 10.
                                # newton parameters are re-read at each 'solve()' calls
                                cyl.setParameter("Newton.MaxRelativeShift", str(maxRelShift))
                                
                                
                                
                            elif n_iter_solve > 50:
                                raise Exception
                            elif (QflowIn_temp < 0 or QflowOut_temp < 0):
                                divVale = 1/0.9
                                print(rank, lId,gId,'qflowIn',qIn,' or qflowout',qOut,
                                      'too low, and/or','maxDt',maxDt_temp,'too high.',
                                      'Decrease manually the lowest and maxDt_temp by',divVale)
                                if QflowIn_temp <= QflowOut_temp:
                                    QflowIn_temp /=divVale
                                    qIn = QflowIn_temp/ (2 * np.pi * self.radii[gId] * l) # [cm3/day] -> [cm /day]
                                    cyl.setInnerFluxCyl(qIn) 
                                else:
                                    QflowOut_temp /= divVale
                                    qOut = QflowOut_temp/(2 * np.pi * self.outer_radii[gId] * l)# [cm3/day] -> [cm /day]
                                    cyl.setOuterFluxCyl(qOut)
                                maxDt_temp /= divVale
                                maxDt_temp = max(maxDt_temp, cyl.ddt)
                            else:
                                raise Exception
                            cyl.reset()
                            n_iter_solve += 1
                            
                            
                    buWAfter_ =  cyl.getWaterVolumesCyl(l)
                    buWAfter = sum(buWAfter_ )
                    buCCAfter_ = cyl.getContentCyl(1, True, l)
                    buTotCAfter = self.getTotCContent(cyl,l)
                    diffW = buWAfter - ( buWBefore + (QflowIn_limited + QflowOut_limited) * dt)
                                                                                                                                   
                                                                                                                   
                    
                    
                    if abs(diffW) > 1e-13:
                        print("abs(diffW) > 1e-13",rank, lId,gId,n_iter_solve,dt, 'end solve QflowIn',QflowIn, 
                              'QflowIn_limited',QflowIn_limited,'diff',QflowIn - QflowIn_limited,
                             'QflowOut',QflowOut, 'QflowOut_limited',QflowOut_limited,
                              'diff',QflowOut - QflowOut_limited, 'buWAfter',buWAfter,
                              'buWBefore',buWBefore, 'diff',diffW)
                        # raise Exception only raise exception at the end of the iteration loop I suppose.


                    self.seg_fluxes_limited[lId] = QflowIn_limited
                    if ( sum(buTotCBefore) + (botVal + topVal + botVal_mucil+ topVal_mucil) * dt) > 0:
                        try:#check solute mass balance
                            assert abs((sum(buTotCAfter) - \
                                        ( sum(buTotCBefore) + (botVal + topVal + botVal_mucil+ \
                                                               topVal_mucil) * dt))/sum(buTotCAfter)*100) < 1.
                        except:
                            print("error mass C")
                            print("before",buTotCBefore )
                            print("after",buTotCAfter )
                            print("CC conservation ",lId, sum(buTotCAfter) ,"before", sum(buTotCBefore) )
                            print("QCflowIn", botVal, botVal_mucil ,"QCflowOut", topVal, topVal_mucil,"dt", dt)
                            print("qCIn",valueBotBC,"qCOut",valueTopBC)
                            print( "l",l,"a_in", self.radii[gId] ,"a_out",self.outer_radii[gId] )
                            print("diff",( sum(buTotCAfter) -( sum(buTotCBefore) + \
                                                              (botVal + topVal + botVal_mucil+ topVal_mucil) * dt))/sum(buTotCAfter)*100)
                            for ncomp in range(self.numComp):
                                print("ncomp_before", ncomp +1, solsBefore[ncomp])
                                print("ncomp_after",ncomp +1 , np.array(cyl.getSolution(ncomp + 1)).flatten() )
                            print("css1 before", css1_before)
                            print("css1 after",  np.array(cyl.base.getCSS1_out()).flatten()/1e6 )
                            raise Exception
                        
                    for ncomp in range(self.numComp):
                        assert (np.array(cyl.getSolution(ncomp + 1)).flatten() >= 0).all()
                    
                except:
                    print( "gId",gId,"l",l,"a_in", self.radii[gId] ,
                          "a_out",self.outer_radii[gId],'maxRelShift',maxRelShift )
                    print("water conservation ",lId, "before", buWBefore , "all cells",buWBefore_)
                    print("QWflowIn", QflowIn ,"QWflowOut", QflowOut,"qWflowIn", qIn ,"qWflowOut", qOut,"dt", dt)
                    print("theta_old", WaterContentOld)
                    print("head_old", SolutionHeadOld)
                    print("theta_new", cyl.getWaterContent())
                    print("head_new", cyl.getSolutionHead())
                    print("CC conservation ",lId,"before", sum(buTotCBefore) )
                    print("QflowInC", botVal, botVal_mucil ,"QflowOutC", topVal, topVal_mucil,"dt", dt)
                    print("concentration", [cyl.getSolution(i+1) for i in range(self.numComp)])
                    print("typeBC, valueBotBC", typeBC, valueBotBC,valueTopBC)
                    myStr = "RhizoMappedSegments.solve: dumux exception with boundaries in flow {:g} cm3/day, out flow {:g} cm3/day, segment radii [{:g}-{:g}] cm"
                    myStr = myStr.format(qIn,qOut, self.radii[gId], self.outer_radii[gId])
                    raise Exception(myStr)
            if verbose:
                print(rank, "FINISHED cyl no ",lId,"/",len(self.cyls),'gId', isinstance(cyl, AirSegment))
            
        
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
            cellIds =  np.array([x for x in cellIds if x >= 0])#take out air segments, cyl
        else:
            cellIds = np.array([troubleShootId])
            verbose = True
            
        
        assert seg_values.shape == (len(self.organTypes), )
        # sum seg_values per voxel
        organTypes = np.array(self.organTypes)
        
        if (organTypes != 2).any():
            try:
                assert (seg_values[np.where(organTypes != 2)] == 0.).all()
            except:
                print(seg_values, organTypes)
                print(seg_values[np.where(organTypes != 2)])
                raise Exception
            # print(seg_values, organTypes)
            # print(seg_values[np.where(organTypes != 2)])
            # raise Exception
            
        
        seg_values_voxel = np.full(len(self.organTypes), 0.)
        splitVals = np.full(len(self.organTypes), 0.)
        
        for cellid in cellIds:
            segIds = self.cell2seg[cellid]
            if soilVals[cellid] != 0:#
                assert (splitVals[segIds] == 0.).all()
                ots = organTypes[segIds]
                rootIds = np.array([sid for sid in segIds if (organTypes[sid] == 2)])
                if len(rootIds) > 0:
                    weightVals = np.full(len(seg_values[segIds]), 0.)

                    if (ots != 2).any():
                        try:
                            assert (seg_values[segIds][np.where(ots != 2)] == 0.).all()
                        except:
                            print(cellid, segIds, ots, seg_values[segIds])
                            raise Exception
                    if verbose:
                        print("soilVals[cellid]", cellid, soilVals[cellid], segIds)
                    if soilVals[cellid] < 0:# goes away from the 1d models
                        seg_values_voxel[rootIds] = sum(seg_values[rootIds])
                        weightVals[np.where(ots == 2)] = seg_values[rootIds] / seg_values_voxel[rootIds]

                        if verbose:
                            print("sum(seg_values[segIds])", seg_values_voxel[segIds], weightVals)
                    else:# goes toward  the 1d models
                        seg_values_voxel[rootIds] = sum(1/seg_values[rootIds])
                        weightVals[np.where(ots == 2)] = (1 / seg_values[rootIds]) / seg_values_voxel[rootIds]
                        if verbose:
                            print("sum(1/seg_values[segIds])", seg_values_voxel[segIds], 
                                  weightVals, sum(1/seg_values[rootIds]))

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

