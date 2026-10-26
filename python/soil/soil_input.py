"""
Same as soil_b2.py, but using an input file instead.

Steady state evaporation with a 3D SPGrid but no resolution in x and y

everything scripted, no input file needed

also works parallel with mpiexec
"""

import matplotlib.pyplot as plt
from analytic_b2 import *
from mpi4py import MPI
from rosi.richards import RichardsWrapper  # Python part
from rosi.rosi_richards import RichardsSP  # C++ part (Dumux binding)

comm = MPI.COMM_WORLD
rank = comm.Get_rank()


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
