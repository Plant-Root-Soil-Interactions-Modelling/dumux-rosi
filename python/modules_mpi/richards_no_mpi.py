from solverbase import SolverWrapper
from richards import RichardsWrapper

import numpy as np


class RichardsNoMPIWrapper(RichardsWrapper):
    """ 
    rewrites all methods using MPI to single process ones
    """

    def __init__(self, base, usemoles):
        super().__init__(base, usemoles)
        self.useMoles = usemoles

    def initialize(self, args_ = [""], verbose = True, doMPI = False):
        """ Writes the Dumux welcome message, and creates the global Dumux parameter tree """
        #print('solverbase_no_mpi:init')
        self.base.initialize(args_, verbose, doMPI)
        #print('solverbase_no_mpi:init_end')
        
    def solve(self, dt:float, maxDt = -1.):
        """ Simulates the problem, the internal Dumux time step ddt is taken from the last time step 
        @param dt      time span [days] 
        @param mxDt    maximal time step [days] 
        """
        #print("solveNoMPI")
        self.base.solveNoMPI(dt * 24.*3600., maxDt * 24.*3600.)  # days -> s

    def getSolutionHead(self, eqIdx=0):
        """Gathers the current solution into rank 0, and converts it into a numpy array (Ndof, neq), 
        model dependent units, [Pa, ...]"""
        self.checkInitialized()
        return self._map((self.base.getSolutionHead(eqIdx)), 0)

    def getWaterContent(self):
        """Gathers the current solution's saturation into rank 0, and converts it into a numpy array (Nc, 1) [1]"""
        self.checkInitialized()
        return self._map((self.base.getWaterContent()), 2)
    
    def getPoints(self):
        """Gathers vertices into rank 0, and converts it into numpy array (Np, 3) [cm]"""
        self.checkInitialized()
        return self._map((self.base.getPoints()), 1) * 100.  # m -> cm

    def getCellCenters(self):
        """Gathers cell centers into rank 0, and converts it into numpy array (Nc, 3) [cm]"""
        self.checkInitialized()
        return self._map((self.base.getCellCenters()), 2) * 100.  # m -> cm

    def getDofCoordinates(self):
        """Gathers dof coorinates into rank 0, and converts it into numpy array (Ndof, 3) [cm]"""
        self.checkInitialized()
        return self._map((self.base.getDofCoordinates()), 0) * 100.  # m -> cm

    def getCells(self):
        """ Gathers dune elements (vtk cells) as list of list of vertex indices (vtk points) (Nc, Number of corners per cell) [1]"""
        return self._map((self.base.getCells()), 2, np.int64)

    def getCellVolumes(self):
        """ Gathers element volumes (Nc, 1) [cm3] """
        return self._map((self.base.getCellVolumes()), 2) * 1.e6  # m3 -> cm3

    def getCellVolumesCyl(self):
        """ Gathers element volumes (Nc, 1) [cm3] """
        return self._map((self.base.getCellVolumesCyl()), 2) * 1.e6  # m3 -> cm3

    def getDofIndices(self):
        """Gathers dof indicds into rank 0, and converts it into numpy array (dof, 1)"""
        self.checkInitialized()
        return  (self.base.getDofIndices())

    def getCellSurfacesCyl(self):
        """ Gathers element volumes (Nc, 1) [cm3] """
        return self._map((self.base.getCellSurfacesCyl()), 2) * 1.e4  # m3 -> cm3

    def getCellVolumesCyl(self):
        """ Gathers element volumes (Nc, 1) [cm3] """
        return self._map((self.base.getCellVolumesCyl()), 2) * 1.e6  # m3 -> cm3

    def getCellCenters(self):
        """Gathers cell centers into rank 0, and converts it into numpy array (Nc, 3) [cm]"""
        self.checkInitialized()
        return self._map((self.base.getCellCenters()), 2) * 100.  # m -> cm

    def getSolution(self, eqIdx=0):
        """Gathers the current solution into rank 0, and converts it into a numpy array (dof, neq), 
        model dependent units [Pa, ...]"""
        self.checkInitialized()
        return self._map((self.base.getSolution(eqIdx)), 0)

    def getAllNeumann(self, eqIdx=0):
        """ Gathers the neuman fluxes into rank 0 as a map with global index as key [cm / day]"""
        assert not self.useMoles
        dics = self.base.getAllNeumann(eqIdx)
        flat_dic = {}
        for d in dics:
            flat_dic.update(d)
        for key, value in flat_dic:
            flat_dic[key] = value / 1000 * 24 * 3600 * 100.  # [kg m-2 s-1] / rho = [m s-1] -> cm / day
        return flat_dic

    def getNetFlux(self, eqIdx=0):
        """ Gathers the net fluxes fir each cell into rank 0 as a map with global index as key [cm3 / day]"""
        assert not self.useMoles
        self.checkInitialized()
        return self._map((self.base.getNetFlux(eqIdx)), 0) * 1000. *24 * 3600  # kg/s -> cm3/day

    def _map(self, x, idType, dtype=np.float64):
        """Converts rows of x to numpy array and maps it to the right indices         
        @param idType 0 dof indices, 1 point (vertex) indices, 2 cell (element) indices   
        """
        if idType == 0:  # auto (dof)
            indices = (self.base.getDofIndices())
        elif idType == 1:  # points
            indices = (self.base.getPointIndices())
        elif idType == 2:  # cells
            indices = (self.base.getCellIndices())
        else:
            raise Exception('PySolverBase._map: idType must be 0, 1, or 2.')
        if len(indices) >0:  # only for rank 0 not empty
            try:
                assert len(indices) == len(x), "_map: indices and values have different length"
            except:
                print(len(indices) , len(x), indices)
                raise Exception
            ndof = max(indices) + 1
            try:
                if isinstance(x[0], (list,type(np.array([])))) :
                    m = len(x[0])
                    p = np.zeros((ndof, m), dtype = dtype)
                    for i in range(0, len(indices)):  #
                        p[indices[i],:] = np.array(x[i], dtype = dtype)
                    if m == 1:
                        p = p.flatten()
                else:
                    p = np.zeros(ndof, dtype = dtype)
                    for i in range(0, len(indices)):  #
                        p[indices[i]] = np.array(x[i], dtype = dtype)
            except:
                print('indices',indices, 'x',x)
                print('x[0]',x[0])
                raise Exception
            return p
        else:
            return np.array([])
