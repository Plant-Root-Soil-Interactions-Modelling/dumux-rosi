import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.interpolate import griddata
import vtk

from mpi4py import MPI

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()


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

    def initialize(self, args = [""]):
        """ Writes the Dumux welcome message, and creates the global Dumux parameter tree """
        self.base.initialize(args)

    def createGrid(self, modelParamGroup = ""):
        """ Creates the Grid and gridGeometry from the global DuMux parameter tree """
        self.base.createGrid(modelParamGroup)

    def createGrid(self, boundsMin, boundsMax, numberOfCells, periodic = "false false false"):
        """ Creates a rectangular grid with a given resolution """
        self.base.createGrid(boundsMin, boundsMax, numberOfCells, periodic)

    def readGrid(self, file :str):
        """ Creates a grid from a file """
        self.base.readGrid(file)

    def getGridBounds(self):
        """  Returns a rectangular bounding box around the grid geometry """
        return self.base.getGridBounds()

    def setParameter(self, key :str, value :str):
        """ Writes a parameter into the global Dumux parameter map """
        self.base.setParameter(key, value)

    def getParameter(self, key :str):
        """ Reads a parameter from the global Dumux parameter map, returns an empty string if value is not set """
        return self.base.getParameter(key)

    def initializeProblem(self):
        """ After the grid is created, the problem can be initialized """
        self.base.initializeProblem()

    def simulate(self, dt :float, maxDt = -1.):
        """ Simulates the problem for time span dt, with initial time step ddt """
        self.base.simulate(dt, maxDt)

    def getPoints(self):
        """Gathers vertices into rank 0, and converts it into numpy array (dof, 3)"""
        self.checkInitialized()
        return self._map(self._flat0(MPI.COMM_WORLD.gather(self.base.getPoints(), root = 0)), 1)

    def getCellCenters(self):
        """Gathers cell centers into rank 0, and converts it into numpy array (dof, 3)"""
        self.checkInitialized()
        return self._map(self._flat0(MPI.COMM_WORLD.gather(self.base.getCellCenters(), root = 0)), 2)

    def getDofCoordinates(self):
        """Gathers dof coorinates into rank 0, and converts it into numpy array (dof, 3)"""
        self.checkInitialized()
        return self._map(self._flat0(MPI.COMM_WORLD.gather(self.base.getDofCoordinates(), root = 0)))

    # def getCells # TODO

    def getDofIndices(self):
        """Gathers dof indicds into rank 0, and converts it into numpy array (dof, 1)"""
        self.checkInitialized()
        return self._flat0(MPI.COMM_WORLD.gather(self.base.getDofIndices(), root = 0))

    def getSolution(self):
        """Gathers the current solution into rank 0, and converts it into a numpy array (dof, neq) """
        self.checkInitialized()
        return self._map(self._flat0(MPI.COMM_WORLD.gather(self.base.getSolution(), root = 0)))

    def pickCell(self, pos :tuple):
        """ Picks a cell and returns its global element cell index """
        return self.base.pickCell(pos)

    def __str__(self):
        """ Solver representation as string """
        return str(self.base)

    def checkInitialized(self):
        """ Checks if the problem was initialized, i.e. initializeProblem() was called """
        return self.base.checkInitialized()

    @property
    def simTime(self):
        """ Current simulation time (read only)"""
        return self.base.simTime

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
        """ last internal time step, i.e. Dumux time step """
        return self.base.ddt

    @ddt.setter
    def ddt(self, value):
        self.base.ddt = value

    def interpolate(self, xi, eq = 0):
        """ interpolates the solution at position x todo: test"""
        self.checkInitialized()
        points = self.getDofCoordinates()
        solution = self.getSolution()
        if rank == 0:
            return griddata(points, solution[:, eq], xi, method = 'linear')
        else:
            return []

    def writeVTK(self, file :str, small :bool = False):
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

    def _map(self, x, type = 0):
        """Converts rows of x to numpy array and maps it to the right indices         
        @param type 0 dof indices, 1 point (vertex) indices, 2 cell (element) indices   
        """
        if type == 0:  # auto (dof)
            indices = self._flat0(MPI.COMM_WORLD.gather(self.base.getDofIndices(), root = 0))
        elif type == 1:  # points
            indices = self._flat0(MPI.COMM_WORLD.gather(self.base.getPointIndices(), root = 0))
        elif type == 2:  # cells
            indices = self._flat0(MPI.COMM_WORLD.gather(self.base.getCellIndices(), root = 0))
        else:
            raise Exception('PySolverBase._map: type must be 0, 1, or 2.')
        if indices:  # only for rank 0 not empty
            assert(len(indices) == len(x))
            ndof = max(indices) + 1
            m = len(x[0])
            p = np.zeros((ndof, m))
            for i in range(0, len(indices)):  #
                p[indices[i], :] = np.array(x[i])
            return p
        else:
            return []

    def _flat0(self, xx):
        """flattens the gathered list in rank 0, empty list for other ranks """
        if rank == 0:
            return [item for sublist in xx for item in sublist]
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
