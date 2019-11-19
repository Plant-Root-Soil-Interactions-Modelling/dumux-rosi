import os
import sys
sys.path.append("../../../build-cmake/rosi_benchmarking/python_solver/")
import richards_yasp_solver as solver

import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.interpolate import griddata

from mpi4py import MPI

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()


class PySolverBase(solver.RichardsYaspSolver):

    # todo getPoints, getCellCenters

    def getDofIndices(self):
        """Gathers dof indicds into rank 0, and converts it into numpy array (dof, 3)"""
        self.checkInitialized()
        indices2 = MPI.COMM_WORLD.gather(super().getDofIndices(), root = 0)
        if rank == 0:
            for i, p in enumerate(indices2):
                print("getDofIndices() rank", i, ":", len(p))
        if rank == 0:
            indices = [item for sublist in indices2 for item in sublist]
        else:
            indices = None
        MPI.COMM_WORLD.Barrier()
        return indices

    def getDofCorrdinates(self):
        """Gathers dof coorinates into rank 0, and converts it into numpy array (dof, 3)"""
        self.checkInitialized()
        points = self._flat0(MPI.COMM_WORLD.gather(super().getDofCoordinates(), root = 0))
        MPI.COMM_WORLD.Barrier()
        return self._map(points)

    def getSolution(self):
        """Gathers the current solution into rank 0, and converts it into a numpy array (dof, neq) """
        self.checkInitialized()
        solution = self._flat0(MPI.COMM_WORLD.gather(super().getSolution(), root = 0))
        MPI.COMM_WORLD.Barrier()
        return self._map(solution)

    def _map(self, x):
        """ converts rows of x to numpy array and maps it to the right indices """
        indices = self._flat0(MPI.COMM_WORLD.gather(super().getDofIndices(), root = 0))
        if indices:  # only for rank 0 not empty
            ndof = max(indices) + 1
            m = len(x[0])
            p = np.zeros((ndof, m))
            for i in range(0, ndof):
                p[indices[i], :] = np.array(x[i])
            return p
        else:
            return []

    def _flat0(self, xx):
        """flattens the gathered list in rank 0, empty list for other ranks """
        if rank == 0:
            x = [item for sublist in xx for item in sublist]
        else:
            x = []
        return x

    def writeVTK(self):
        """ todo """
        pass

    def interpolate(self, xi, eq = 0):
        """ interpolates the solution at position x """
        self.checkInitialized()
        points = self.getDofCoordinates()  # todo convert to numpy
        MPI.COMM_WORLD.Barrier()
        solution = self.getSolution()
        if rank == 0:
            return griddata(points, solution[:, eq], xi, method = 'linear')
        else:
            return []

    def plot_interpolated_Z(self, eq = 0, conv = lambda x: x, xlabel = "Pressure", xy = (0, 0)):
        """ plots the current solution along the z axis (not working, todo) """
        self.checkInitialized()
        bounds = self.getGridBounds()
        z_ = np.linspace(bounds[2], bounds[5], 200)
        y_ = np.ones(z_.shape) * xy[1]
        x_ = np.ones(z_.shape) * xy[0]
        xi = np.hstack((x_, y_, z_))
        sol = self.interpolate(xi, eq)  # for this all processes are needed
        if rank == 0:
            plt.plot(conv(sol), z_ * 100, "*")  #
            plt.xlabel(xlabel)
            plt.ylabel("Depth (cm)")
            plt.show()

    def plotZ(self, x, xlabel = "Pressure"):
        """ plots x along the z axis """
        self.checkInitialized()
        points = np.array(self.getDofCorrdinates())
        if rank == 0:
            plt.plot(x, points[:, 2] * 100, "*")
            plt.xlabel(xlabel)
            plt.ylabel("Depth (cm)")
            plt.show()

