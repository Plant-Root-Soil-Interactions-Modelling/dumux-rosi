import os
import sys
sys.path.append("../../../build-cmake/rosi_benchmarking/python_solver/")
import solverbase as solver

import scipy


class RichardsYaspSolver(solver.PySolverBase):

    # todo should belong to solverbase.py in

    def setVanGenuchtenParameter(self):
        pass

    def setBCtop(self):
        pass

    def setBCbot(self):
        pass

