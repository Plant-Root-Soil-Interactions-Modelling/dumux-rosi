import os
import sys
sys.path.append("../../../build-cmake/rosi_benchmarking/python_solver/")
import richards_yasp_solver as solver

import numpy as np
import scipy

from mpi4py import MPI

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()


class PySolverBase(solver.RichardsYaspSolver):

    # todo getPoints, getCellCenters
    
    def getDofIndices(self):
        """Gathers dof indicds into rank 0, and converts it into numpy array (dof, 3)"""
        self.checkInitialized()
        comm = MPI.COMM_WORLD
        indices2 = comm.gather(super().getDofIndices(), root = 0)
        if rank == 0: 
            for i,p in enumerate(indices2):
                print("getDofIndices() rank", i, ":", len(p))
        if rank == 0:
            indices = [item for sublist in indices2 for item in sublist]
        else: 
            indices = None
        return indices   
    
    def getDofCorrdinates(self):
        """Gathers dof coorinates into rank 0, and converts it into numpy array (dof, 3)"""
        self.checkInitialized()
        comm = MPI.COMM_WORLD
        points2 = comm.gather(super().getDofCoordinates(), root = 0)
        points = [item for sublist in points2 for item in sublist]
        return self._map(points)

    def getSolution(self):
        """Gathers the current solution into rank 0, and converts it into a numpy array (dof, neq) """
        self.checkInitialized()
        comm = MPI.COMM_WORLD
        solution2 = comm.gather(super().getSolution(), root = 0)
        solution = [item for sublist in solution2 for item in sublist]
        return self._map(solution)

    def _map(self, x):
        """ maps to the right indices """
        comm = MPI.COMM_WORLD
        indices2 = comm.gather(super().getDofIndices(), root = 0)
        indices = [item for sublist in indices2 for item in sublist]
        ndof = max(indices) + 1
        m = len(x[0])
        p = np.zeros((ndof + 1, m))
        for i in range(0, ndof):
            p[indices[i], :] = np.array(x[i])
        return p

    def writeVTK(self):
        pass

    def interpolate(self, xi, eq = 0):
        """ interpolates the solution at position x """
        points = self.getDof()  # todo convert to numpy
        return scipy.interpolate.griddata(points, solution[eq], method = 'linear')

