import sys
sys.path.append("../../../build-cmake/rosi_benchmarking/python_solver/")

from dumux_rosi_cyl import RichardsCylFoam  # C++ part (Dumux binding)
from solver.richards import RichardsWrapper  # Python part

import matplotlib.pyplot as plt
import numpy as np
import os

from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

""" 
Cylindrical 1D model 

everything scripted, no input file needed

also works parallel with mpiexec
""" 
    
cpp_base = RichardsCylFoam()
s = RichardsWrapper(cpp_base)
s.initialize()

loam = [0.045, 0.43, 0.04, 1.6, 50]

s.createGrid([0.02], [0.6], [10])  # [cm]

s.setHomogeneousIC(-100.)  # cm pressure head
s.setTopBC("noFlux")  #  [cm/day]
s.setBotBC("constantFlux", -0.00002) # [cm/day] * m radius TODO
s.setVGParameters([loam])
s.initializeProblem()
 
if rank == 0:
    print(s)
 
times = np.array([0, 288000, 576000, 864000]) / (3600 * 24)  
s.ddt = 1.e-5  # days
 
for dt in np.diff(times):
    if rank == 0:
        print("*****", "external time step", dt, " d, simulation time", s.simTime, "d, internal time step", s.ddt, "d")
    s.solve(dt)
 
points = s.getDofCoordinates()
x = s.getSolution()
plt.plot(points[:], RichardsWrapper.toHead(x), "r*", label = "dumux")

os.chdir("../../../build-cmake/rosi_benchmarking/soil_richards/python")
data = np.loadtxt("cylinder_1d_Comsol.txt")
z_comsol = data[:,0]
h_comsol = data[:,-1]
plt.plot(z_comsol+0.02, h_comsol, "b", label = "comsol")
plt.legend()

plt.show()
