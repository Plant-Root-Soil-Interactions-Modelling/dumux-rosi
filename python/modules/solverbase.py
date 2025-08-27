import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.interpolate import griddata
from mpi4py import MPI; comm = MPI.COMM_WORLD; size = comm.Get_size(); rank = comm.Get_rank()
import helpful


class SolverWrapper():
    """ Additional functionality to the C++ Python binding. 
    
        The C++ class derived from SolverBase is wrapped, 
        SolverBase functions are passed to C++ 
        
        Additionally, contains mainly methods that are easier to write in Python, 
        e.g. MPI communication, writeVTK, interpolate        
    """

    def __init__(self, base):
        """ @param base is the C++ base class that is wrapped. """
        self.base = base
        self.molarMassWat = 18. # [g/mol]
        self.densityWat = 1. #[g/cm3]
        self.molarDensityWat =  self.densityWat / self.molarMassWat # [mol wat/cm3 wat] 
        #only needed if we have element in the solid phase. no effect on result if soil is static
        self.molarBulkDensity = 1. # [mol soil minerals / cm3 bulk soil density]
        self.bulkDensity = 1. # [g soil minerals / cm3 bulk soil density]
        self.length = 1. # [cm] optional parameter for the 1D axisymmetric models to go from 2d to 3d
        self.results_dir = "./"

    def initialize(self, args_ = [""], verbose = False, doMPI_ = True):
        """ Creates the global Dumux parameter tree 
        
            @params verbose: level of verbosity of the dumux object [bool]
            @params doMPI_: do we use multi processes in the dumux object [bool]
            If the python code is run in parallel, we can use doMPI_=False for the
            1d models only (FoamGrid). Does not work for the other (3d) grids.
        """
        with helpful.StdoutRedirector() as redirector:
            try:
                self.base.initialize(args_, verbose, doMPI=doMPI_)
            except Exception as e:
                target_filepath = self.results_dir + 'stdcout_cpp_'+str(rank)+'.txt'
                with open(target_filepath, 'w') as f:
                    f.write(redirector.buffer)
                raise Exception
            
    def setMaxTimeStepSize(self, maxDt):
        """
            change the maximum inner time step which can be tested by dumux
        """
        self.base.setMaxTimeStepSize(maxDt)
        
    def createGridFromInput(self, modelParamGroup = ""):
        """ Creates the Grid and gridGeometry from the global DuMux parameter tree """
        self.base.createGrid(modelParamGroup)

    def is_periodic(self):
        """ returns if periodic in X and Y direction """
        return self.base.periodic

    def createGrid(self, boundsMin, boundsMax, numberOfCells, periodic = False):
        """ Creates a rectangular grid with a given resolution             
            @param boundsMin        domain corner [cm]
            @param boundsMax        domain corner [cm]        
            @param numberOfCells    resoultion [1]
            @param periodic         If true, the domain is periodic in x and y, not in z 
        """
        self.base.createGrid(np.array(boundsMin) / 100., np.array(boundsMax) / 100., np.array(numberOfCells), periodic)  # cm -> m
        
        self.numberOfCellsTot = np.prod(self.numberOfCells)
        if self.dimWorld == 3:
            self.numberOfFacesTot = self.numberOfCellsTot * 6
        elif self.dimWorld == 1:
            self.numberOfFacesTot = self.numberOfCellsTot * 2
        else:
            raise Exception
            

    def createGrid1d(self, points, length = 1.):
        """ Creates a 1D grid with a given resolution             
            @param points           1D coordinates of the vertices [cm]
            @param length           length of the domain, perpendicular to the 1D grid. 
        """
        p = []
        for v in points:
            p.append([v / 100.])  # cm -> m
        
        with helpful.StdoutRedirector() as redirector:
            try:
                self.base.createGrid1d(p)
            except Exception as e:
                target_filepath = self.results_dir + 'stdcout_cpp_'+str(rank)+'.txt'
                with open(target_filepath, 'w') as f:
                    f.write(redirector.buffer)
                raise Exception
            
        self.length = length

        self.numberOfCellsTot = len(points) -1
        self.numberOfFacesTot = self.numberOfCellsTot * 2
        
#     def createGrid3d(self, points, p0):
#         """ todo
#         """
#         p, p0_ = [], []
#         for v in points:
#             p.append(list(v / 100.))  # cm -> m
#         for v in p0:
#             p0_.append(list(v / 100.))  # cm -> m
#         self.base.createGrid3d(p, p0_)

    def readGrid(self, file:str):
        """ Creates a grid from a file (e.g. dgf or msh)"""
        self.base.readGrid(file)

    def getGridBounds(self):
        """  Returns a rectangular bounding box around the grid geometry [cm] """
        return np.array(self.base.getGridBounds()) * 100.  # m -> cm

    def setParameter(self, key:str, value:str):
        """ Writes a parameter into the global Dumux parameter map """
        self.base.setParameter(key, value)

    def getParameter(self, key:str):
        """ Reads a parameter from the global Dumux parameter map, returns an empty string if value is not set """
        return self.base.getParameter(key)

    def initializeProblem(self, maxDt = -1.):
        """ After the grid is created, the problem can be initialized 
        @param maxDt    maximal time step [days] 
        """
        with helpful.StdoutRedirector() as redirector:
            try:
                self.base.initializeProblem(maxDt * 24.*3600.)
            except Exception as e:
                target_filepath = self.results_dir + 'stdcout_cpp_'+str(rank)+'.txt'
                with open(target_filepath, 'w') as f:
                    f.write(redirector.buffer)
                raise Exception

        ##       saves shape data of the grid  (to limit the communication needed between threads)
        # local indices
        self.dofIndices_   = self.base.getDofIndices()
        self.pointIndices_ = self.base.getPointIndices()
        self.cellIndices_  = self.base.getCellIndices()
        # all indices
        self.dofIndices   = np.asarray(self._flat0(self.gather(self.dofIndices_)), np.int64)
        self.pointIndices = np.asarray(self._flat0(self.gather(self.pointIndices_)), np.int64) 
        self.cellIndices  = np.asarray(self._flat0(self.gather(self.cellIndices_)), np.int64) 
        
        if self.dimWorld == 3:
            self.cellsVertex = self._map(self._flat0(self.gather(self.base.getCells())), 2, np.int64)
        else:
            self.cellsVertex = [None]
            
        # volumes, surface
        self.cellVolumes_ =np.array( self.base.getCellVolumes()) * 1.e6  # m3 -> cm3
        self.cellVolumes = self._map(self._flat0(self.gather(self.cellVolumes_)), 2)   # m2 -> cm2
        
        # coordinates
        self.pointCoords = self._map(self._flat0(self.gather(self.base.getPoints())), 1) * 100.  # m -> cm
        self.cellCenters = self._map(self._flat0(self.gather(self.base.getCellCenters())), 2) * 100.  # m -> cm
        self.dofCoordinates = self._map(self._flat0(self.gather(self.base.getDofCoordinates())), 0) * 100.  # m -> cm
        if rank == 0:
            self.lowerBoundary = np.array([self.base.onLowerBoundary( cc ) for cc in self.cellCenters])
            self.upperBoundary = np.array([self.base.onUpperBoundary( cc ) for cc in self.cellCenters])
        else:    
            self.lowerBoundary = []
            self.upperBoundary = []

    def setInitialCondition(self, ic, eqIdx = 0):
        """ Sets the initial conditions for all global elements, processes take from the shared @param ic """
        self.base.setInitialCondition(ic, eqIdx)

    def setInitialConditionHead(self, ic):
        """ Sets the initial conditions for all global elements, processes take from the shared @param ic """
        self.base.setInitialConditionHead(ic)

    def solve(self, dt:float, doMPIsolve_ = True, saveInnerFluxes_ = False):
        """ Simulates the problem, the internal Dumux time step ddt is taken from the last time step 
        @param dt      time span [days] 
        
            # get error:
            # Successive accumulation of data on coarse levels only works with ParMETIS installed.  
            # Fell back to accumulation to one domain on coarsest level
            # see matrixhierarchy.hh l695
            # TODO: suppress message? change input parameters?
        """
        self.base.solve(dt * 24.*3600., doMPIsolve = doMPIsolve_, saveInnerDumuxValues = saveInnerFluxes_)  # days -> s
                
    def solveSteadyState(self):
        """ Finds the steady state of the problem """
        self.base.solveSteadyState()

    def getLowerBoundary(self):
        return self.lowerBoundary
        
    def getUpperBoundary(self):
        return self.upperBoundary

    def getPoints(self):
        """Gathers vertices into rank 0, and converts it into numpy array (Np, 3) [cm]"""
        self.checkGridInitialized()
        return self.pointCoords

    def getPoints_(self):
        """nompi version of """
        self.checkGridInitialized()
        if rank > 0:
            return []
        else:
            return np.array(self.base.getPoints()) * 100.  # m -> cm

    def getCellCenters(self):
        """Gathers cell centers into rank 0, and converts it into numpy array (Nc, 3) [cm]"""
        self.checkGridInitialized()
        return self.cellCenters

    def getCellCenters_(self):
        """nompi version of """
        self.checkGridInitialized()
        return np.array(self.base.getCellCenters()) * 100.  # m -> cm

    def getDofCoordinates(self):
        """ dof coorinates into rank 0, and converts it into numpy array (Ndof, 3) [cm]"""
        self.checkGridInitialized()
        return self.dofCoordinates

    def getDofCoordinates_(self):
        """nompi version of """
        self.checkGridInitialized()
        return np.array(self.base.getDofCoordinates()) * 100.  # m -> cm

    def getCells(self):
        """ dune elements (vtk cells) as list of list of vertex indices (vtk points) (Nc, Number of corners per cell) [1]"""
        return self.cellsVertex

    def getCells_(self):
        """nompi version of """
        return np.array(self.base.getCells(), dtype = np.int64)

    def getCellSurfacesCyl(self):
        """ Gathers element volumes (Nc, 1) [cm3] """
        return self._map(self._flat0(self.gather(self.base.getCellSurfacesCyl(), root = 0)), 2) * 1.e4  # m2 -> cm2
        
    def getCellSurfacesCyl_(self):
        """nompi version of  """
        return np.array(self.base.getCellSurfacesCyl()) * 1.e4  # m2 -> cm2
        
    def getCellVolumes(self):
        """ Gathers element volumes (Nc, 1) [cm3] """
        if self.dimWorld == 3 :
            return self.cellVolumes
        elif self.dimWorld == 1 :
            return self.getCellVolumesCyl()

    def getCellVolumes_(self):
        """nompi version of  """
        if self.dimWorld == 3 :
            return self.cellVolumes_ 
        elif self.dimWorld == 1 :
            return self.getCellVolumesCyl_()
        
    def getCellVolumesCyl(self):
        """ Gathers element volumes (Nc, 1) [cm3] """
        return self._map(self._flat0(self.gather(self.getCellVolumesCyl_(), root = 0)), 2) 

    def getCellVolumesCyl_(self):
        """ Gathers element volumes (Nc, 1) [cm3] 
            uses the optional input parameter of @see solverbase::createGrid1d
        """
        return self.getCellSurfacesCyl_() * self.length # [m2 -> cm2] * cm

    def getWaterVolumes(self):
        """Returns total water volume of the domain [cm3]"""
        self.checkGridInitialized()
        vols = self.getCellVolumes() #cm3 scv
        watCont = self.getWaterContent()# # cm3 wat/cm3 scv        
        return np.multiply(vols , watCont  )  
        
    # def quad, int or something (over all domain)

    def getDofIndices(self):
        """Gathers dof indicds into rank 0, and converts it into numpy array (dof, 1)"""
        self.checkGridInitialized()
        if rank > 0:
            return []
        else:
            return self.dofIndices

    def getDofIndices_(self):
        """nompi version of  """
        self.checkGridInitialized()
        return self.dofIndices_

    def getPointIndices(self):
        """Gathers dof indicds into rank 0, and converts it into numpy array (dof, 1)"""
        if rank > 0:
            return []
        else:
            return self.pointIndices
        
    def getPointIndices_(self):
        """nompi version of  """
        self.checkGridInitialized()
        return self.pointIndices_
        
    def getCellIndices(self):
        """Gathers dof indicds into rank 0, and converts it into numpy array (dof, 1)"""
        if rank > 0:
            return []
        else:
            return self.cellIndices
        
    def getSolution(self, eqIdx = 0):
        """Gathers the current solution into rank 0, and converts it into a numpy array (dof, neq), 
        model dependent units [Pa, ...]"""
        self.checkGridInitialized()
        return self._map(self._flat0(self.gather(self.base.getSolution(eqIdx), root = 0)), 0)

    def getSolution_(self, eqIdx = 0):
        """nompi version of  """
        self.checkGridInitialized()
        return np.array(self.base.getSolution(eqIdx))

    def setSolution(self, x, eqIdx = 0):
        """Gathers the current solution into rank 0, and converts it into a numpy array (dof, neq), 
        model dependent units [Pa, ...]"""
        self.checkGridInitialized()
        self.base.setSolution(x, eqIdx)

    def getSolutionAt(self, gIdx, eqIdx = 0):
        """Returns the current solution at a cell index, model dependent units [Pa, ...]"""
        return self.base.getSolutionAt(gIdx, eqIdx)

    def getNeumann(self, gIdx, eqIdx = 0):
        """ Gathers the neuman fluxes into rank 0 as a map with global index as key [cm / day]"""
        # print("gIdx, eqIdx",gIdx, eqIdx)
        return self.base.getNeumann(gIdx, eqIdx) / 1000 * 24 * 3600 * 100.  # [kg m-2 s-1] / rho = [m s-1] -> cm / day

    def getAllNeumann(self, eqIdx = 0):
        """ Gathers the neuman fluxes into rank 0 as a map with global index as key [cm / day]"""
        dics = self.gather(self.base.getAllNeumann(eqIdx), root = 0)
        flat_dic = {}
        for d in dics:
            flat_dic.update(d)
        for key, value in flat_dic:
            flat_dic[key] = value / 1000 * 24 * 3600 * 100.  # [kg m-2 s-1] / rho = [m s-1] -> cm / day
        return flat_dic

    def getAllNeumann_(self, eqIdx = 0):
        """nompi version of (TODO is that working?)"""
        dics = self.base.getAllNeumann(eqIdx)
        flat_dic = {}
        for d in dics:
            flat_dic.update(d)
        for key, value in flat_dic:
            flat_dic[key] = value / 1000 * 24 * 3600 * 100.  # [kg m-2 s-1] / rho = [m s-1] -> cm / day
        return flat_dic

    def getVelocities(self, eqIdx = 0):
        """ Gathers the net fluxes fir each cell into rank 0 as a map with global index as key [cm3 / day]"""
        self.checkGridInitialized()
        return self._map(self._flat0(self.gather(self.base.getVelocities(eqIdx), root = 0)), 0) * 100. *24 * 3600  # m/s -> cm/day

    def getVelocities_(self, eqIdx = 0):
        """nompi version of """
        self.checkGridInitialized()
        return np.array(self.base.getVelocities(eqIdx)) * 100. *24 * 3600  # cm/s -> cm/day
        
    def getFaceSurfaces_(self, length = None):
        if self.dimWorld == 1:
            if length is None:
                raise Exception('getFaceSurfaces_: length parameter required to get flux of 1d domain')
            return np.array(self.base.getCylFaceCoordinates()) * 100 * 2 * np.pi * length
        else:
            return np.array( self.base.getFaceSurfaces()) * 1e4 # cm2
        
    def getFace2CellIdx_(self):
        face2CellIdx = np.array(self.base.getFace2CellIdx())
        notGhost = face2CellIdx >= 0
        return face2CellIdx[notGhost]
        
    def getFlowsPerCell(self, eqIdx = 0, length = None):
        """ Gathers the net sources fir each cell into rank 0 as a map with global index as key [cm3/ day or kg/ day or mol / day]"""
        return self._map(self._flat0(self.gather(self.getFlowsPerCell_(eqIdx, length), root = 0)), 2)
 
    def getFlowsPerCell_(self, eqIdx = 0, length = None):
        """nompi version of
           [cm3/day] or [g/day] or [mol/day]
        """
        self.checkGridInitialized()
        scvfFlows = np.array(self.getFlowsPerFace_(eqIdx, length)  ) # [cm3/day] 
        face2CellIdx = np.array(self.getFace2CellIdx_())
        localcellsnumber = len(self.cellIndices_)
        scvFlows = np.zeros(localcellsnumber)
        scvFlows[np.unique(face2CellIdx)] = np.array([
             sum(scvfFlows[face2CellIdx == cellindx]) for cellindx in np.unique(face2CellIdx) 
            ])
        return scvFlows
            
    def getFlowsPerFace_(self, eqIdx = 0, length = None):
        """nompi version of 
            [cm3/day] or [g/day] or [mol/day]
        """
        self.checkGridInitialized() 
        scvfFlows = self.getBoundaryFlowsPerFace_(eqIdx, length)      
        scvfFlows += self.getCell2CellFlowsPerFace_(eqIdx, length)     
        return scvfFlows 
        
    def getUpperBoundaryFlowPerCell(self, eqIdx = 0, length = None):
        """nompi version of
           [cm3/day] or [g/day] or [mol/day]
        """
        scvFlows = self.getBoundaryFlowsPerCell(eqIdx, length)
        upperBoundary = self.getUpperBoundary()
        if rank == 0:
            scvFlows[~ upperBoundary] = 0.
        return scvFlows
        
        
    def getLowerBoundaryFlowsPerCell(self, eqIdx = 0, length = None):
        """nompi version of
           [cm3/day] or [g/day] or [mol/day]
        """
        scvFlows = self.getBoundaryFlowsPerCell(eqIdx, length)
        lowerBoundary = self.getLowerBoundary()
        if rank == 0:
            scvFlows[~ lowerBoundary] = 0.
        return scvFlows
        
    def getLowerBoundaryFluxesPerCell(self, eqIdx = 0, length = None):
        """nompi version of
           [cm3/day] or [g/day] or [mol/day]
        """
        scvFluxes = self.getBoundaryFluxesPerCell(eqIdx, length)
        lowerBoundary = self.getLowerBoundary()
        if rank == 0:
            scvFluxes[~ lowerBoundary] = 0.
        return scvFluxes
        
    def getBoundaryFluxesPerCell(self, eqIdx = 0, length = None):
        """nompi version of
           [cm3/day] or [g/day] or [mol/day]
        """
        return self._map(self._flat0(self.gather(self.getBoundaryFluxesPerCell_(eqIdx, length), root = 0)), 2)
 
    def getBoundaryFluxesPerCell_(self, eqIdx = 0, length = None):
        """nompi version of
           [cm3/day] or [g/day] or [mol/day]
        """
        self.checkGridInitialized()
        scvfFluxes = np.array(self.getBoundaryFluxesPerFace_(eqIdx, length)  ) # [cm3/day] 
        face2CellIdx = np.array(self.getFace2CellIdx_())
        localcellsnumber = len(self.cellIndices_)
        scvFluxes = np.zeros(localcellsnumber)
        scvFluxes[np.unique(face2CellIdx)] = np.array([
             sum(scvfFluxes[face2CellIdx == cellindx]) for cellindx in np.unique(face2CellIdx) 
            ])
        return scvFluxes
        
    def getBoundaryFlowsPerCell(self, eqIdx = 0, length = None):
        """nompi version of
           [cm3/day] or [g/day] or [mol/day]
        """
        return self._map(self._flat0(self.gather(self.getBoundaryFlowsPerCell_(eqIdx, length), root = 0)), 2)
        
    def getBoundaryFlowsPerCell_(self, eqIdx = 0, length = None):
        """nompi version of
           [cm3/day] or [g/day] or [mol/day]
        """
        self.checkGridInitialized()
        scvfFlows = np.array(self.getBoundaryFlowsPerFace_(eqIdx, length)  ) # [cm3/day] 
        face2CellIdx = np.array(self.getFace2CellIdx_())
        localcellsnumber = len(self.cellIndices_)
        scvFlows = np.zeros(localcellsnumber)
        scvFlows[np.unique(face2CellIdx)] = np.array([
             sum(scvfFlows[face2CellIdx == cellindx]) for cellindx in np.unique(face2CellIdx) 
            ])
        return scvFlows
        
    def getBoundaryFluxesPerFace_(self, eqIdx = 0, length = None):
        """  [cm3/cm2/day] """
        self.checkGridInitialized()
        unitChange = 24 * 3600 / 1e4 # [ kgOrmol/m2/s -> kgOrmol/cm2/day]
        if eqIdx == 0:
            if self.useMoles:
                unitChange *=  18.068 # [mol/cm2/day -> cm3/cm2/day]
            else:
                unitChange *= 1e3  # [kg/cm2/day] -> [cm3/cm2/day] 
        elif not self.useMoles:
            unitChange *= 1e3  # [kg/cm2/day] -> [g/cm2/day]
            
        surfaces = self.getFaceSurfaces_(length)  # cm2 
        notGhost = surfaces > 0
        return (np.array(self.base.getScvfBoundaryFluxes()[eqIdx]))[notGhost] * unitChange
    
    def getBoundaryFlowsPerFace_(self, eqIdx = 0, length = None):
        """  [cm3/day] """
        self.checkGridInitialized()
        scvfFluxes = self.getBoundaryFluxesPerFace_(eqIdx, length)
        surfaces = self.getFaceSurfaces_(length)  # cm2 
        notGhost = surfaces > 0
        surfaces = surfaces[notGhost]
        return scvfFluxes * surfaces
            
    def getCell2CellFlowsPerFace_(self, eqIdx = 0, length = None):
        """  [cm3/day] """
        self.checkGridInitialized()
        unitChange = 24 * 3600 # [ kgOrmol/s -> kgOrmol/day]
        if eqIdx == 0:
            if self.useMoles:
                unitChange *=  18.068 # [mol/day -> cm3/day]
            else:
                unitChange *= 1e3  # [kg/cm2/day] -> [cm3/cm2/day] 
        elif not self.useMoles:
            unitChange *= 1e3  # [kg/cm2/day] -> [g/cm2/day]
            
        surfaces = self.getFaceSurfaces_(length)        
        notGhost = surfaces > 0
        return np.array(self.base.getScvfInnerFluxes()[eqIdx])[notGhost] * unitChange
        
    def pickCell(self, pos):
        """ Picks a cell and returns its global element cell index """
        return self.base.pickCell(np.array(pos) / 100.)  # cm -> m

    def pick(self, x):
        """ Picks a cell and returns its global element cell index """
        return self.base.pick(np.array(x) / 100.)  # cm -> m

    def __str__(self):
        """ Solver representation as string """
        return str(self.base)

    def checkGridInitialized(self):
        """ Checks if the grid was created, i.e. createGrid() or createGrid1d() was called """
        return self.base.checkGridInitialized()

    def checkProblemInitialized(self, throwError = True):
        """ Checks if the problem was created, i.e. initializeProblem() was called """
        if self.base.checkProblemInitialized():
            return True
        elif throwError:
            raise Exception('problem was not created, call initializeProblem()')
        else:
            return False

    @property
    def dimWorld(self):
        """ dimention (1 or 3 """
        if (self.base.dimWorld != 1) and (self.base.dimWorld != 3):
            raise Exception(f'recieved unexpected dimention: {self.base.dimWorld}')
        return self.base.dimWorld 
        
    @property
    def useMoles(self):
        """ dumux units in moles [True] or kg [False] """
        return self.base.useMoles() 
        
    @property
    def simTime(self):
        """ Current simulation time (read only) [days]"""
        return self.base.simTime / 24. / 3600.  # s->days

    @property
    def rank(self):
        """ MPI rank (read only) """
        return self.base.rank

    @property
    def maxRank(self):
        """ number of MPI processes (read only) """
        return self.base.maxRank

    @property
    def ddt(self):
        """ last internal time step, i.e. Dumux time step [days]"""
        return self.base.ddt / 24. / 3600.  # s-> days

    @ddt.setter
    def ddt(self, value):
        """ sets internal time step, i.e. Dumux time step [days]"""
        self.base.ddt = value * 24.*3600.  # days -> s

    def interpolate(self, xi, eq = 0):
        """ interpolates the solution at position ix [cm],
        model dependent units
        todo: test"""
        self.checkGridInitialized()
        points = self.getDofCoordinates()
        solution = self.getSolution()
        if rank == 0:
            yi = np.zeros((xi.shape[0]))
            for i in range(0, xi.shape[0]):
                yi[i] = griddata(points, solution[:, eq], xi / 100., method = 'linear', rescale = True)  # cm -> m
            return yi
        else:
            return []

    def interpolate_(self, xi, points, solution, eq = None):
        """ interpolates the solution at position ix [cm] (same es interpolate_ but passes points and solution),
        model dependent units """
        self.checkGridInitialized()
        if rank == 0:
            yi = np.zeros((xi.shape[0]))
            for i in range(0, xi.shape[0]):
                if eq:
                    yi[i] = griddata(points, solution[:, eq], xi / 100., method = 'linear')  # cm -> m
                else:
                    yi[i] = griddata(points, solution[:], xi / 100., method = 'linear')  # cm -> m
            return yi
        else:
            return []

    def interpolateNN(self, xi, eq = None):
        """ solution at the points xi (todo currently works only for CCTpfa)"""
        self.checkGridInitialized()
        solution = self.getSolution()
        if rank == 0:
            y = np.zeros((xi.shape[0]))
        else:
            y = []
        for i in range(0, xi.shape[0]):
            idx = self.pickCell(xi[i,:])
            if eq and rank == 0:
                y[i] = solution[idx, eq]
            elif rank == 0:
                y[i] = solution[idx]
        return y

    def writeVTK(self, file:str, small:bool = False):
        """writes a vtk file (todo additional fields) 
        @param file 
        @param small Determines if data are compressed and stroed binary  TODO
        """
        points = self.getDofCoordinates()
        # cells = self.getCells()
        if self.rank == 0:
            pd = vtk.vtkPolyData()
            pd.SetPoints(_vtkPoints(points))

            writer = vtk.vtkXMLPolyDataWriter()
            writer.SetFileName(file)
            writer.SetInputData(pd)
            writer.SetDataModeToBinary()
            writer.SetCompressorTypeToZLib()
            writer.Write()

    def _map(self, x, type_, dtype = np.float64):
        """Converts rows of x to numpy array and maps it to the right indices         
        @param type_ 0 dof indices, 1 point (vertex) indices, 2 cell (element) indices   
        """
        if type_ == 0:  # auto (dof)
            indices = self._flat0(self.gather(self.base.getDofIndices(), root = 0))
        elif type_ == 1:  # points
            indices = self._flat0(self.gather(self.base.getPointIndices(), root = 0))
        elif type_ == 2:  # cells
            indices = self._flat0(self.gather(self.base.getCellIndices(), root = 0))
        else:
            raise Exception('PySolverBase._map: type_ must be 0, 1, or 2.')
        if len(indices) > 0:  # only for rank 0 not empty
            try:
                assert len(indices) == len(x), "_map: indices and values have different length"
            except:
                print(len(indices) , len(x), indices)
                raise Exception
            ndof = max(indices) + 1
            if isinstance(x[0], (list, type(np.array([])))):
                m = len(x[0])
                p = np.zeros((ndof, m), dtype = dtype)
                for i in range(0, len(indices)):  #
                    p[indices[i],:] = np.array(x[i], dtype = dtype)
            else:  # option to get array of shape (N,) instead of (N,1)
                p = np.zeros(ndof, dtype = dtype)
                for i in range(0, len(indices)):  #
                    p[indices[i]] = np.array(x[i], dtype = dtype)
            return p
        else:
            return 0

    def _flat0(self, xx, dtype_=None):
        """ flattens the gathered list in rank 0, empty list for other ranks """
        if rank == 0:
            return np.array([item for sublist in xx for item in sublist],dtype = dtype_)
        else:
            return np.array([])

    def gather(self, data2gather, root = 0):
        """ to make the overloading of functions by richards_no_mpi easier
            @see richards_no_mpi::gather()
        """
        return comm.gather(data2gather, root) 
        
    
    def allgatherv(self,data2share, keepShape = False, data2share_type_default = float): 
        """ own allgather() function for vectors
            use it if the vectors exchanged between the threads might be of
            different size
            @param data2share: (local) vector of the thread to share
            @param keepShape: the output vector has to be of the same shape as data2share
            @data2share_type_default: type of the data in data2share
        """
        try:
            assert isinstance(data2share, (list, type(np.array([]))))
        except:
            error = ('on rank {}: '.format(rank)+
                'solverbase::allgatherv() the data to share should be of type list or np.array.'+
                r'not {}.'.format(type(data2share))     
            ) 
            print(error)
            raise Exception
            
        data2share = np.array(data2share)
            
        # get shape and size of local array
        if len((data2share).shape) == 2:
            local_size = (data2share).shape[0] * (data2share).shape[1]
            shape0 = (data2share).shape[0]
            shape1 = (data2share).shape[1]
            data2share_type = type(data2share[0][0])
        elif len((data2share).shape) == 1:
            local_size = (data2share).shape[0] 
            shape0 = (data2share).shape[0]
            shape1 = 0
            if len(data2share) > 0:
                data2share_type = type(data2share[0])
            else:
                data2share_type = data2share_type_default
        else:
            raise Exception  
        
        # data2share needs to use floats for Allgatherv with MPI.DOUBLE to work.
        data2share = np.array(data2share, dtype = np.float64)
        
        # other data needed by comm.Allgatherv
        all_sizes = np.array(comm.allgather(local_size))
        work_size = sum(all_sizes)
        all_data2share = np.zeros(work_size)

        offsets = np.zeros(len(all_sizes), dtype=np.int64)
        offsets[1:]=np.cumsum(all_sizes)[:-1]
        all_sizes =tuple(all_sizes)
        offsets =tuple( offsets)      
        
        # share the vectors
        comm.Allgatherv( [data2share.reshape(-1), MPI.DOUBLE],[all_data2share,all_sizes,offsets,MPI.DOUBLE])
        
        all_data2share = np.array(all_data2share, dtype =data2share_type)
        
        if keepShape:
            if shape1 > 0:
                all_data2share = all_data2share.reshape( size,shape0,shape1)
            else:
                all_data2share = all_data2share.reshape(size, shape0)
        else:
            if shape1 > 0:
                all_data2share = all_data2share.reshape(-1,shape1)
            else:
                all_data2share = all_data2share.reshape(-1)
        
        return all_data2share
        
    def createLinearSolver(self):
        """
            manually (re)create nonlinear solver. 
            useful to implement new solver parameters
        """
        self.base.createLinearSolver()
    
    def createNewtonSolver(self):
        """
            manually (re)create nonlinear solver. 
            useful to implement new solver parameters
        """
        self.base.createNewtonSolver()
            
    @staticmethod
    def _vtkPoints(p):
        """ Creates vtkPoints from an numpy array"""
        assert(p.shape(1) == 3)
        da = vtk.vtkDataArray.CreateDataArray(vtk.VTK_DOUBLE)
        da.SetNumberOfComponents(3)  # vtk point dimension is always 3
        da.SetNumberOfTuples(p.shape[0])
        for i in range(0, p.shape[0]):
            da.InsertTuple3(i, p[i, 0], p[i, 1], p[i, 2])
        points = vtk.vtkPoints()
        points.SetData(da)
        return points

    @staticmethod
    def _vtkCells(t):
        """ Creates vtkCells from an numpy array"""
        cellArray = vtk.vtkCellArray()
        if t.shape[1] == 2:
            Simplex = vtk.vtkLine
        elif t.shape[1] == 4:
            Simplex = vtk.vtkTetra
        elif t.shape[1] == 8:
            Simplex = vtk.Hexaedron
        else:
            raise Exception('PySolverBase._vtkCells: do not know what to do with {} vertices'.format(t.shape[1]))
        for vert in t:
            tetra = Simplex()
            for i, v in enumerate(vert):
                tetra.GetPointIds().SetId(i, int(v))
            cellArray.InsertNextCell(tetra)
        return cellArray

    @staticmethod
    def _vtkData(data):
        """ todo """
        pass

    @staticmethod
    def to_pa(ph):
        """ Converts cm pressure head to Pascal [kg/ (m s^2)] """
        return 1.e5 + ph / 100 * 1000. * 9.81;

    @staticmethod
    def to_head(p):
        """ Converts Pascal [kg/ (m s^2)] to cm pressure head """
        return (p - 1.e5) * 100. / 1000. / 9.81;
        
    @property
    def numComp(self):
        """Get the number of components evaluated (including water) ."""
        return self.base.numComp()  
        
    @property
    def numSoluteComp(self):
        """Get the number of components evaluated (not counting water)."""
        return self.base.numComp() - 1 
    
    @property
    def numFluidComp(self):
        """Get the number of components in the fluid phase (including water)."""
        return self.base.numFluidComp()
        
    @property
    def numDissolvedSoluteComp(self):
        """Get the number of disolved solutes."""
        return self.base.numFluidComp() - 1
