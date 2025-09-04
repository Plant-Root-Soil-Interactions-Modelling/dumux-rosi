from solverbase import SolverWrapper
from richards_no_mpi import RichardsNoMPIWrapper

import numpy as np


class RichardsNoMPIFlatWrapper(RichardsNoMPIWrapper):
    """ 
    get the outputs as flattened arrays for 1D models, without re-writing RichardsNoMPIWrapper
    """

    def __init__(self, base):
        super().__init__(base)

    def getSolutionHead(self, eqIdx=0):
        """Gathers the current solution into rank 0, and converts it into a numpy array (Ndof,), 
        model dependent units, [Pa, ...]"""
        return super().getSolutionHead(eqIdx).flatten()
        
    def getSolution(self, eqIdx=0):
        """Gathers the current solution into rank 0, and converts it into a numpy array (Ndof,), 
        model dependent units"""
        return super().getSolution(eqIdx).flatten()

    def getWaterContent(self):
        """Gathers the current solution's saturation into rank 0, and converts it into a numpy array (Nc,) [1]"""
        return super().getWaterContent().flatten()
    
    def getPoints(self):
        """Gathers vertices into rank 0, and converts it into numpy array (Np,) [cm]"""
        return super().getPoints().flatten() 

    def getCellCenters(self):
        """Gathers cell centers into rank 0, and converts it into numpy array (Nc,) [cm]"""
        return super().getCellCenters().flatten()

    def getDofCoordinates(self):
        """Gathers dof coorinates into rank 0, and converts it into numpy array (Ndof,) [cm]"""
        return super().getDofCoordinates().flatten() 

    def getCellVolumes(self):
        return super().getCellVolumes().flatten() 

    def getSource(self, eqIdx = 0):
        return super().getSource(eqIdx).flatten() 
        
    def getContent(self,eqIdx):
        return super().getContent(eqIdx).flatten() 