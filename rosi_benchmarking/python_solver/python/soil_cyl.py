import sys
sys.path.append("../../../build-cmake/rosi_benchmarking/python_solver/")

from rosi_richards_cyl import RichardsCylFoam  # C++ part (Dumux binding)
from solver.richards import RichardsWrapper  # Python part
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()  # MPI

import matplotlib.pyplot as plt
import numpy as np
import os

""" 
Cylindrical 1D model, water movement only (DuMux)

DuMux everything scripted, no input file needed, also works parallel with mpiexec
"""

cpp_base = RichardsCylFoam()
s = RichardsWrapper(cpp_base)
s.initialize()

loam = [0.045, 0.43, 0.04, 1.6, 50]

s.createGrid([0.02], [0.6], [100])  # [cm]

s.setHomogeneousIC(-100.)  # cm pressure head
s.setOuterBC("fluxCyl", 0.)  #  [cm/day]
s.setInnerBC("fluxCyl", -0.1)  # [cm/day]
s.setVGParameters([loam])
s.initializeProblem()
s.setCriticalPressure(-15000)  # cm pressure head

if rank == 0:
    print(s)

times = [0., 10, 20]  # days
s.ddt = 1.e-5

col = ["r*", "b*", "g*", "c*", "m*", "y*", ]

for i, dt in enumerate(np.diff(times)):

    if rank == 0:
        print("*****", "external time step", dt, " d, simulation time", s.simTime, "d, internal time step", s.ddt, "d")

    s.solve(dt)
    points = s.getDofCoordinates()
    x = s.getSolutionHead()
    plt.plot(points[:], x, col[i % len(col)], label = "dumux {} days".format(s.simTime))

os.chdir("../../../build-cmake/rosi_benchmarking/soil_richards/python")
data = np.loadtxt("cylinder_1d_Comsol.txt", skiprows = 8)
z_comsol = data[:, 0]
plt.plot(z_comsol + 0.02, data[:, 25], "k", label = "comsol 10 days")
plt.plot(z_comsol + 0.02, data[:, -1], "k:", label = "comsol 20 days")

plt.xlabel("distance from root axis (cm)")
plt.xlabel("soil matric potential (cm)")
plt.legend()
plt.show()
