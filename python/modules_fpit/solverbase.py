import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.interpolate import griddata
from mpi4py import MPI; comm = MPI.COMM_WORLD; size = comm.Get_size(); rank = comm.Get_rank()


class SolverWrapper():
    """ Additional functionality to the C++ Python binding. 
    
        The C++ class derived from SolverBase is wrapped, 
        SolverBase functions are passed to C++ 
        
        Additionally, contains mainly methods that are easier to write in Python, 
        e.g. MPI communication, writeVTK, interpolate        
    """

    def __init__(self, base, usemoles):
        """ @param base is the C++ base class that is wrapped. """
        self.base = base
        self.useMoles = usemoles
        self.solidDensity = 0
        self.solidMolarMass=0
        self.solidMolDensity=0
        self.bulkDensity_m3 =0 #mol / m3 bulk soil
        self.molarMassWat = 18. # [g/mol]
        self.densityWat_m3 = 1e6 #[g/m3]
        # [mol/m3] = [g/m3] /  [g/mol] 
        self.molarDensityWat_m3 =  self.densityWat_m3 / self.molarMassWat # [mol wat/m3 wat] 
        

    def initialize(self, args_ = [""], verbose = True):
        """ Writes the Dumux welcome message, and creates the global Dumux parameter tree """
        self.base.initialize(args_, verbose)
        
    def reset(self):
        """ reset solution vector to value before solve function """
        self.base.reset()
        
    def setVerbose(self, verbose:int):
        """ set verbose level """
        self.base.setVerbose(verbose)
        

    def createGridFromInput(self, modelParamGroup = ""):
        """ Creates the Grid and gridGeometry from the global DuMux parameter tree """
        self.base.createGrid(modelParamGroup)

    def createGrid(self, boundsMin, boundsMax, numberOfCells, periodic = False):
        """ Creates a rectangular grid with a given resolution             
            @param boundsMin        domain corner [cm]
            @param boundsMax        domain corner [cm]        
            @param numberOfCells    resoultion [1]
            @param periodic         If true, the domain is periodic in x and y, not in z 
        """
        self.base.createGrid(np.array(boundsMin) / 100., np.array(boundsMax) / 100., np.array(numberOfCells), periodic)  # cm -> m

    def createGrid1d(self, points):
        """ todo
        """
        p = []
        for v in points:
            p.append([v / 100.])  # cm -> m
        self.base.createGrid1d(p)

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

    def initializeProblem(self, rank_ = 0):
        """ After the grid is created, the problem can be initialized """
        self.base.initializeProblem()
        cell = self.getCellCenters()
        if len(cell) > 0: # cell float if rank >0
            cell0 = cell[0]
            if (isinstance(cell0, float) or isinstance(cell0, np.float) ):
                self.dimWorld = 1
            else:
                self.dimWorld = len(cell0)        
            assert ((self.dimWorld  == 1) or ( self.dimWorld  == 3))
        else:
            self.dimWorld = 0


    def setInitialCondition(self, ic, eqIdx = 0):
        """ Sets the initial conditions for all global elements, processes take from the shared @param ic """
        self.base.setInitialCondition(ic, eqIdx)

    def setInitialConditionHead(self, ic):
        """ Sets the initial conditions for all global elements, processes take from the shared @param ic """
        self.base.setInitialConditionHead(ic)

    def solve(self, dt:float, maxDt = -1.):
        """ Simulates the problem, the internal Dumux time step ddt is taken from the last time step 
        @param dt      time span [days] 
        @param mxDt    maximal time step [days] 
        """
        self.base.solve(dt * 24.*3600., maxDt * 24.*3600.)  # days -> s

    def solveSteadyState(self):
        """ Finds the steady state of the problem """
        self.base.solveSteadyState()

    def getPoints(self):
        """Gathers vertices into rank 0, and converts it into numpy array (Np, 3) [cm]"""
        self.checkInitialized()
        return self._map(self._flat0(comm.gather(self.base.getPoints(), root = 0)), 1) * 100.  # m -> cm

    def getPoints_(self):
        """nompi version of """
        self.checkInitialized()
        assert max_rank == 1
        return np.array(self.base.getPoints()) * 100.  # m -> cm

    def getCellCenters(self):
        """Gathers cell centers into rank 0, and converts it into numpy array (Nc, 3) [cm]"""
        self.checkInitialized()
        return self._map(self._flat0(comm.gather(self.base.getCellCenters(), root = 0)), 2) * 100.  # m -> cm

    def getCellCenters_(self):
        """nompi version of """
        self.checkInitialized()
        #assert max_rank == 1
        return np.array(self.base.getCellCenters()) * 100.  # m -> cm

    def getDofCoordinates(self):
        """Gathers dof coorinates into rank 0, and converts it into numpy array (Ndof, 3) [cm]"""
        self.checkInitialized()
        return self._map(self._flat0(comm.gather(self.base.getDofCoordinates(), root = 0)), 0) * 100.  # m -> cm

    def getDofCoordinates_(self):
        """nompi version of """
        self.checkInitialized()
        assert max_rank == 1
        return np.array(self.base.getDofCoordinates()) * 100.  # m -> cm

    def getCells(self):
        """ Gathers dune elements (vtk cells) as list of list of vertex indices (vtk points) (Nc, Number of corners per cell) [1]"""
        return self._map(self._flat0(comm.gather(self.base.getCells(), root = 0)), 2, np.int64)

    def getCells_(self):
        """nompi version of """
        assert max_rank == 1
        return np.array(self.base.getCells(), dtype = np.int64)

    def getCellSurfacesCyl(self):
        """ Gathers element volumes (Nc, 1) [cm3] """
        return self._map(self._flat0(comm.gather(self.base.getCellSurfacesCyl(), root = 0)), 2) * 1.e4  # m3 -> cm3

    def getCellSurfacesCyl_(self):
        """nompi version of  """
        assert max_rank == 1
        return np.array(self.base.getCellSurfacesCyl()) * 1.e4  # m2 -> cm2
        
    def getCellVolumes(self):
        """ Gathers element volumes (Nc, 1) [cm3] """
        return self._map(self._flat0(comm.gather(self.base.getCellVolumes(), root = 0)), 2) * 1.e6  # m2 -> cm2

    def getCellVolumes_(self):
        """nompi version of  """
        return np.array(self.base.getCellVolumes()) * 1.e6  # m3 -> cm3

    def getCellVolumesCyl(self):
        """ Gathers element volumes (Nc, 1) [cm3] """
        return self._map(self._flat0(comm.gather(self.base.getCellVolumesCyl(), root = 0)), 2) * 1.e6  # m3 -> cm3

    def getCellVolumesCyl_(self):
        """nompi version of  """
        return np.array(self.base.getCellVolumesCyl()) * 1.e6  # m3 -> cm3

    # def quad, int or something (over all domain)

    def getDofIndices(self):
        """Gathers dof indicds into rank 0, and converts it into numpy array (dof, 1)"""
        self.checkInitialized()
        return self._flat0(comm.gather(self.base.getDofIndices(), root = 0))

    def getDofIndices_(self):
        """nompi version of  """
        self.checkInitialized()
        return np.array(self.base.getDofIndices(), dtype = np.int64)

    def getSolution(self, eqIdx = 0):
        """Gathers the current solution into rank 0, and converts it into a numpy array (dof, neq), 
        model dependent units [Pa, ...]"""
        self.checkInitialized()
        return self._map(self._flat0(comm.gather(self.base.getSolution(eqIdx), root = 0)), 0)

    def getAvgDensity_(self):
        """nompi version of  """
        self.checkInitialized()
        return np.array(self.base.getAvgDensity())
        
    def getSolution_(self, eqIdx = 0):
        """nompi version of  """
        self.checkInitialized()
        return np.array(self.base.getSolution(eqIdx))

    def getSolutionAt(self, gIdx, eqIdx = 0):
        """Returns the current solution at a cell index, model dependent units [Pa, ...]"""
        return self.base.getSolutionAt(gIdx, eqIdx)

    def getNeumann(self, gIdx, eqIdx = 0):
        """ Gathers the neuman fluxes into rank 0 as a map with global index as key [cm / day]"""
        assert not self.useMoles
        return self.base.getNeumann(gIdx, eqIdx) / 1000 * 24 * 3600 * 100.  # [kg m-2 s-1] / rho = [m s-1] -> cm / day

    def getAllNeumann(self, eqIdx = 0):
        """ Gathers the neuman fluxes into rank 0 as a map with global index as key [cm / day]"""
        assert not self.useMoles
        dics = comm.gather(self.base.getAllNeumann(eqIdx), root = 0)
        flat_dic = {}
        for d in dics:
            flat_dic.update(d)
        for key, value in flat_dic:
            flat_dic[key] = value / 1000 * 24 * 3600 * 100.  # [kg m-2 s-1] / rho = [m s-1] -> cm / day
        return flat_dic

    def getAllNeumann_(self, eqIdx = 0):
        """nompi version of (TODO is that working?)"""
        assert not self.useMoles
        dics = self.base.getAllNeumann(eqIdx)
        flat_dic = {}
        for d in dics:
            flat_dic.update(d)
        for key, value in flat_dic:
            flat_dic[key] = value / 1000 * 24 * 3600 * 100.  # [kg m-2 s-1] / rho = [m s-1] -> cm / day
        return flat_dic

    def getNetFlux(self, eqIdx = 0):
        """ Gathers the net fluxes fir each cell into rank 0 as a map with global index as key [cm3 / day]"""
        assert not self.useMoles
        self.checkInitialized()
        return self._map(self._flat0(comm.gather(self.base.getNetFlux(eqIdx), root = 0)), 0) * 1000. *24 * 3600  # kg/s -> cm3/day

    def getNetFlux_(self, eqIdx = 0):
        """nompi version of """
        assert not self.useMoles
        self.checkInitialized()
        return np.array(self.base.getNetFlux(eqIdx)) * 1000. *24 * 3600  # kg/s -> cm3/day

    def pickCell(self, pos):
        """ Picks a cell and returns its global element cell index """
        return self.base.pickCell(np.array(pos) / 100.)  # cm -> m

    def pick(self, x):
        """ Picks a cell and returns its global element cell index """
        return self.base.pick(np.array(x) / 100.)  # cm -> m
        
        
    def pick_(self,coordCell):
        bounds = self.getGridBounds(); min_b = bounds[:3]; max_b = bounds[3:]
        cell_number_ = self.numberOfCells
        ratioDist = (coordCell - min_b)/(max_b - min_b)
        if ((ratioDist > 1.) |(ratioDist < 0.) ).any():#not in the domain
            return -1
        id_rows = ratioDist*cell_number_
        onCellOuterFace = np.where((id_rows == np.round(id_rows)) & (id_rows > 0) )[0]
        id_rows[onCellOuterFace] -= 1
        id_rows= np.floor(id_rows)
        id_cell = id_rows[2]*(cell_number_[0]*cell_number_[1])+id_rows[1]*cell_number_[0]+id_rows[0]
        return int(id_cell) 

    def __str__(self):
        """ Solver representation as string """
        return str(self.base)

    def checkInitialized(self):
        """ Checks if the problem was initialized, i.e. initializeProblem() was called """
        return self.base.checkInitialized()

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
        self.checkInitialized()
        points = self.getDofCoordinates()
        solution = self.getSolution()
        if rank == 0:
            yi = np.zeros((xi.shape[0]))
            for i in range(0, xi.shape[0]):
                yi[i] = griddata(points, solution[:, eq], xi / 100., method = 'linear', rescale = True)  # cm -> m
            return yi
        else:
            return []

    def interpolate_(self, xi, points, solution, eq = 0):
        """ interpolates the solution at position ix [cm] (same es interpolate_ but passes points and solution),
        model dependent units """
        self.checkInitialized()
        if rank == 0:
            yi = np.zeros((xi.shape[0]))
            for i in range(0, xi.shape[0]):
                yi[i] = griddata(points, solution[:, eq], xi / 100., method = 'linear')  # cm -> m
            return yi
        else:
            return []

    def interpolateNN(self, xi, eq = 0):
        """ solution at the points xi (todo currently works only for CCTpfa)"""
        self.checkInitialized()
        solution = self.getSolution()
        if rank == 0:
            y = np.zeros((xi.shape[0]))
        else:
            y = []
        for i in range(0, xi.shape[0]):
            idx = self.pickCell(xi[i,:])  # cm -> m
            if rank == 0:
                y[i] = solution[idx, eq]
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
            indices = self._flat0(comm.gather(self.base.getDofIndices(), root = 0))
        elif type_ == 1:  # points
            indices = self._flat0(comm.gather(self.base.getPointIndices(), root = 0))
        elif type_ == 2:  # cells
            indices = self._flat0(comm.gather(self.base.getCellIndices(), root = 0))
        else:
            raise Exception('PySolverBase._map: type_ must be 0, 1, or 2.')
        if len(indices) >0:  # only for rank 0 not empty
            try:
                assert len(indices) == len(x), "_map: indices and values have different length"
            except:
                print(len(indices) , len(x), indices)
                raise Exception
            ndof = max(indices) + 1
            if isinstance(x[0], (list,type(np.array([])))) :
                m = len(x[0])
                p = np.zeros((ndof, m), dtype = dtype)
                for i in range(0, len(indices)):  #
                    p[indices[i],:] = np.array(x[i], dtype = dtype)
            else:
                p = np.zeros(ndof, dtype = dtype)
                for i in range(0, len(indices)):  #
                    p[indices[i]] = np.array(x[i], dtype = dtype)
            return p#.flatten()
        else:
            #return array instead of float to be able to have same object type no matter what
            return np.array([])

    def _flat0(self, xx):
        """flattens the gathered list in rank 0, empty list for other ranks """
        #print("solverbase::_flat0", rank, xx)
        if rank == 0:
            return np.array([item for sublist in xx for item in sublist])
        else:
            return []

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