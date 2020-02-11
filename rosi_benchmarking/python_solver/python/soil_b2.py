import sys
sys.path.append("../../../build-cmake/rosi_benchmarking/python_solver/")

from dumux_rosi import RichardsSP  # C++ part (Dumux binding)
from solver.richards import RichardsWrapper  # Python part

import matplotlib.pyplot as plt

from analytic_b2 import *

from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

""" 
Steady state evaporation with a 3D SPGrid but no resolution in x and y (for speed)

everything scripted, no input file needed

also works parallel with mpiexec
"""

cpp_base = RichardsSP()
s = RichardsWrapper(cpp_base)
s.initialize()

N = 52  # resolution
s.createGrid([-5., -5., -53.], [5., 5., 0.], [1, 1, N])  # [cm]
s.setHomogeneousIC(-53, True)  # inital guess in hydraulic equilibrium
s.setTopBC("constantFlux", -0.5)  #  [cm/day]
s.setBotBC("constantPressure", 0.)  # cm pressure head
loam = [0.08, 0.43, 0.04, 1.6, 50]
s.setVGParameters([loam])
s.initializeProblem()

s.solveSteadyState()

z = s.getDofCoordinates()
x = s.getSolution()

if rank == 0:
    print("\n", s)  # some info
    plt.plot(RichardsWrapper.toHead(x), z[:, 2], "r*")
    plt.show()

