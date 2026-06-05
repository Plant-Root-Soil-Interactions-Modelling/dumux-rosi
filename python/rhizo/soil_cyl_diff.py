import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src")

from rosi_richards_cyl import RichardsCylFoam  # C++ part (Dumux binding)
from rosi_richards_cyl import RichardsCylFoamAna  # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part

import matplotlib.pyplot as plt
import numpy as np
import os
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()

""" 
Cylindrical 1D model, diffusion only (DuMux), Michaelis Menten

everything scripted, no input file needed, also works parallel with mpiexec
"""
def setupModel(s):
    s.initialize()

    loam = [0.045, 0.43, 0.04, 1.6, 50]

    s.createGrid([0.02], [0.6], [500])  # [cm]

    s.setHomogeneousIC(-100.)  # cm pressure head
    s.setOuterBC("noflux")  #  [cm/day]
    s.setInnerBC("fluxCyl", -0.1)  # [cm/day] -0.1

    s.setVGParameters([loam])
    s.initializeProblem()
    s.setCriticalPressure(-15000)  # cm pressure head

    s.ddt = 1.e-5
    return s


s1 = setupModel(RichardsWrapper(RichardsCylFoam()))
s2 = setupModel(RichardsWrapper(RichardsCylFoamAna()))


fig, (ax1, ax2) = plt.subplots(1, 2)

times = [0., 10., 20.]  # days

col = ["r*", "b*", "g*", "c*", "m*", "y*", ]

for i, dt in enumerate(np.diff(times)):

    if rank == 0:
        print("*****", "external time step", dt, " d, simulation time", s1.simTime, s2.simTime, "d, internal time step", s1.ddt, s2.ddt, "d")
    print('solve s2')
    s2.solve(dt)
    print('solve s1')
    s1.solve(dt)

    points = s.getDofCoordinates()

    x1 = s1.getSolutionHead()
    x2 = s2.getSolutionHead()

    ax1.plot(points[:], x1, col[i % len(col)], label = "dumux num {} days".format(s1.simTime))
    ax1.plot(points[:], x2, col[i % len(col)], label = "dumux ana {} days".format(s2.simTime))


ax1.set_xlabel('distance from the root axis (cm)')
ax1.set_ylabel('psi (hPa or cm)')

ax1.legend()

plt.show()
