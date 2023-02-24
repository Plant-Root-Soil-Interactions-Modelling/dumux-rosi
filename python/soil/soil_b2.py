import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src")

from rosi_richards import RichardsSP  # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part
from analytic_b2 import *

import matplotlib.pyplot as plt; from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()

""" 
Steady state evaporation with a 3D SPGrid but no resolution in x and y 

everything scripted, no input file needed

also works parallel with mpiexec
"""
s = RichardsWrapper(RichardsSP())
s.initialize()

N = 53  # resolution
s.createGrid([-5., -5., -53.], [5., 5., 0.], [1, 1, N])  # [cm]

s.setHomogeneousIC(-53, True)  # initial guess in hydraulic equilibrium
s.setTopBC("constantFlux", -0.5)  #  [cm/day]
s.setBotBC("constantPressure", 0.)  # cm pressure head
loam = [0.08, 0.43, 0.04, 1.6, 50]
s.setVGParameters([loam])
s.initializeProblem()

s.solveSteadyState()

z = s.getDofCoordinates()
x = s.getSolutionHead()

if rank == 0:
    plt.plot(x, z[:, 2], "r*")
    plt.show()

