import os
import sys
sys.path.append("../../../build-cmake/rosi_benchmarking/python_solver/")
import richards_sp_solver as solver

import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.interpolate import griddata
import vtk

from mpi4py import MPI

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

CSolverBase = solver.RichardsYaspSolver


class PySolverBase(CSolverBase):
    """ Additional functionality to the C++ Python binding. 
    
        The class defined in the Binding acts as Base class 
        
        Contains mainly methods that are easier to write in Python, 
        e.g. MPI communication, writeVTK, interpolate        
    """

    def getDofIndices(self):
        """Gathers dof indicds into rank 0, and converts it into numpy array (dof, 1)"""
        self.checkInitialized()
        return self._flat0(MPI.COMM_WORLD.gather(super().getDofIndices(), root = 0))

    def getPoints(self):
        """Gathers dof coorinates into rank 0, and converts it into numpy array (dof, 3)"""
        self.checkInitialized()
        return self._map(self._flat0(MPI.COMM_WORLD.gather(super().getPoints(), root = 0)), 1)

    def getCellCenters(self):
        """Gathers dof coorinates into rank 0, and converts it into numpy array (dof, 3)"""
        self.checkInitialized()
        return self._map(self._flat0(MPI.COMM_WORLD.gather(super().getCellCenters(), root = 0)), 2)

    def getDofCoordinates(self):
        """Gathers dof coorinates into rank 0, and converts it into numpy array (dof, 3)"""
        self.checkInitialized()
        return self._map(self._flat0(MPI.COMM_WORLD.gather(super().getDofCoordinates(), root = 0)))

    def getSolution(self):
        """Gathers the current solution into rank 0, and converts it into a numpy array (dof, neq) """
        self.checkInitialized()
        return self._map(self._flat0(MPI.COMM_WORLD.gather(super().getSolution(), root = 0)))

    def _map(self, x, type = 0):
        """Converts rows of x to numpy array and maps it to the right indices         
        @param type 0 dof indices, 1 point (vertex) indices, 2 cell (element) indices   
        """
        if type == 0:  # auto (dof)
            indices = self._flat0(MPI.COMM_WORLD.gather(super().getDofIndices(), root = 0))
        elif type == 1:  # points
            indices = self._flat0(MPI.COMM_WORLD.gather(super().getPointIndices(), root = 0))
        elif type == 2:  # cells
            indices = self._flat0(MPI.COMM_WORLD.gather(super().getCellIndices(), root = 0))
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
        @param small Determines if data are compressed and stroed binary  
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
