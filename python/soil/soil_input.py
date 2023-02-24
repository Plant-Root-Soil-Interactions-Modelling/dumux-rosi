import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src")

from rosi_richards import RichardsSP  # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part
from analytic_b2 import *

import matplotlib.pyplot as plt; from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()

""" 
Same as soil_b2.py, but using an input file instead.

Steady state evaporation with a 3D SPGrid but no resolution in x and y 

everything scripted, no input file needed

also works parallel with mpiexec
"""
s = RichardsWrapper(RichardsSP())
s.initialize(["", "b2_3d.input"])
s.createGridFromInput("Soil")
s.initializeProblem()

s.solveSteadyState()

z = s.getDofCoordinates()
x = s.getSolutionHead()

if rank == 0:
    plt.plot(x, z[:, 2], "r*")
    plt.show()

