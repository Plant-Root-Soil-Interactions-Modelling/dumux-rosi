# import sys; sys.path.append("../modules/"); sys.path.append("../modules/fv/"); sys.path.append("../../../CPlantBox/"); sys.path.append("../../../CPlantBox/src");
# sys.path.append("../../build-cmake/cpp/python_binding/")

import plantbox as pb
import functional.xylem_flux as xylem_flux
import sys
from fv.fv_grid import *
import fv.fv_richards as rich  # local pure Python cylindrical models
import functional.van_genuchten as vg

import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import norm
from scipy.interpolate import griddata as gd
from mpl_toolkits.mplot3d import Axes3D
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()
import multiprocessing
from multiprocessing import Process, active_children
import psutil
from threading import Thread 
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank(); max_rank = comm.Get_size()

#TODO: adapt to represent better air domain (air flux, co2, h2o, T gradient...
# soil-air analogy:
# rhizosphere => leaf boundary layer
# bulk soil: canopy space + above-canopy space
# outer BC of bulk soil: value at point of  measurment (Xm aboveground)

class AirSegment():#solve later the atmospheric flow also via dumux?
    """
        dummy class to use insteade of the dumux wrapper for shoot segments and     
        root segments aboveground. 
        Same function names as the ones in rhichards but they don t do anything
        used to avoid getting erros when iterating through the segments
        in the rhizo_modelsPlant class
    """
    def __init__(self, a_in:float = 1., a_out:float = 1.1, length = 1.):
        
        self.a_in = a_in #cm
        
        self.a_out = a_out
        self.points = np.array([a_in,self.a_out ]) #need at least 2 nodes == 1 cell
        self.grid = FVGrid1Dcyl(self.points)
        self.length = length
        
        #to use?
        self.n = self.grid.n_cells
        self.x0 = np.zeros((self.n,))  # solution of last time step [cm]
        self.sim_time = 0.
        self.bc = { }  # boundary conditions, map with key (cell_id, face_id) containing a list of values
        self.sources = np.zeros((self.n,))  # [cm3 / cm3]
        
    def solve(self, dt, maxDt = None):
        print("dummy solve function for airSegments. does nothing (for now?)")
    def reset(self):
        pass
    def resetManual(self):
        pass
    def saveManual(self):
        pass
    def save(self):
        pass
    #water
    #get
    def getInnerHead(self,shift=0, val:float=-90000):  # [cm]
        return self.get_inner_head(val)
    def get_inner_head(self, val):
        return val
    def getInnerFlux(self,val=0):
        return self.get_inner_flux()
    def get_inner_flux(self):#return val equal to the "proposed flux" for that segment
        return self.innerFlux
    def getWaterVolume(self):
        return 0
    def getWaterContent(self):
        return np.array([0])
    def getSolutionHead(self):
        return np.array([0])
    def getWaterVolumesCyl(self ):
        return np.array([0])
    def getKrw(self):# could return here the air resistance instead of having it in the photosynthesis class
        return np.array([np.Inf])
        
    #set
    def setInnerFluxCyl(self,val) :#flux computed by the plant (transpiration for leaf, 0 else)
        self.innerFlux = val
    def setOuterFluxCyl(self,val) :
        pass
    def setInnerMatricPotential(self,val):
        pass
    def setRootSystemBC(self,*arg):
        pass
    
    #solute
    #get
    def getInnerSolutes(self,shift=0, compId = 1):
        return 0
    def get_inner_concentration(self):
        return 0
    def getCSS1_out_real(self):
        return np.array([0])
    def getCSS1_out_th(self):
        return np.array([0])
    def getSolution_(self,val):
        return np.array([0])
    def getSolution(self,val):
        return np.array([0])
    def getContent(self,idComp):
        return np.array([0])
    #set
    def setOuterBC_solute(self,*arg):
        pass
    def setInnerBC_solute(self,*arg):
        pass
    
    #all
    #def solve(self,*arg):
    #    print("solve aire segment")
        
    #shape
    def getDofCoordinates(self):
        raise Exception
        return np.array([0])
    def getCellCenters(self):
        return (self.points[1:] +  self.points[:-1])/2
    def getPoints(self):
        return self.points
    def getCellSurfacesCyl(self):
        """nompi version of  """
        return np.array([np.pi * (self.a_out*self.a_out - self.a_in*self.a_in)])  # cm2
    def getCellSurfacesCyl_(self):
        """nompi version of  """
        return np.array([np.pi * (self.a_out*self.a_out - self.a_in*self.a_in)])  # cm2
    
    def getCellVolumes(self):
        return self.getCellSurfacesCyl() * self.length
        



