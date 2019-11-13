import os
import sys
sys.path.append("../../../build-cmake/rosi_benchmarking/soil_richards/")
import richards_yasp_solver as solver

import scipy


class RichardsYaspSolver(solver.RichardsYaspSolver):

    def setVanGenuchtenParameter(self):
        pass

    def setBCtop(self):
        pass

    def setBCbot(self):
        pass

    def writeVTK(self):
        pass

    def interpolate(self, xi, eq = 0):
        """ interpolates the solution at position x """
        points = self.getDof()  # todo convert to numpy
        return scipy.interpolate.griddata(points, solution[eq], method = 'linear')

