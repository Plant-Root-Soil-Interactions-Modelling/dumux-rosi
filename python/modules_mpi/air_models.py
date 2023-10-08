# import sys; sys.path.append("../modules/"); sys.path.append("../modules/fv/"); sys.path.append("../../../CPlantBox/"); sys.path.append("../../../CPlantBox/src");
# sys.path.append("../../build-cmake/cpp/python_binding/")

import plantbox as pb
import functional.xylem_flux as xylem_flux
import sys
from rosi_richards_cyl import RichardsCylFoam  # C++ part (Dumux binding) of cylindrcial model
from rosi_richardsnc_cyl import RichardsNCCylFoam  # C++ part (Dumux binding)
from richards_no_mpi import RichardsNoMPIWrapper  # Python part of cylindrcial model (a single cylindrical model is not allowed to run in parallel)
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

#TODO: adapt to represent better air domain (air flux, co2, h2o, T gradient...
# soil-air analogy:
# rhizosphere => leaf boundary layer
# bulk soil: canopy space + above-canopy space
# outer BC of bulk soil: value at point of  measurment (Xm aboveground)

class AirSegment():#solve later the atmospheric flow also via dumux?
    def __init__(self, psi_air, a_in):
        self.psi_air = psi_air
        points = np.array([a_in,a_in+a_in*0.1]) #need at least 2 nodes == 1 cell
        self.grid = FVGrid1Dcyl(points)
        
        #to use?
        self.n = self.grid.n_cells
        self.x0 = np.zeros((self.n,))  # solution of last time step [cm]
        self.sim_time = 0.
        self.bc = { }  # boundary conditions, map with key (cell_id, face_id) containing a list of values
        self.sources = np.zeros((self.n,))  # [cm3 / cm3]
        
    #water
    #get
    def getInnerHead(self,shift=0):  # [cm]
        return self.get_inner_head()
    def get_inner_head(self):
        return self.psi_air
    def getInnerFlux(self,val=0):
        return self.get_inner_flux()
    def get_inner_flux(self):
        return 0
    def getWaterVolume(self):
        return 0
    def getWaterContent(self):
        return np.array([0])
    def getSolutionHead(self):
        return np.array([0])
        
    #set
    def setInnerFluxCyl(self,val) :
        pass
    def setOuterFluxCyl(self,val) :
        pass
    def setInnerMatricPotential(self,val):
        pass
    def setRootSystemBC(self,*arg):
        pass
    
    #solute
    #get
    def getInnerSolutes(self,shift=0):
        return 0
    def get_inner_concentration(self):
        return 0
    def getSolution_(self,val):
        return np.array([0])
    #set
    def setOuterBC_solute(self,*arg):
        pass
    def setInnerBC_solute(self,*arg):
        pass
    
    #all
    def solve(self,*arg):
        pass
    def getDofCoordinates(self):
        return np.array([0])
    def getPoints(self):
        return self.grid.nodes



