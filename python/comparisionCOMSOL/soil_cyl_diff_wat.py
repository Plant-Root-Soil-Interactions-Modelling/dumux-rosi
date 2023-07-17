import sys; sys.path.append("../modules/"); sys.path.append("../../../CPlantBox/");  sys.path.append("../../build-cmake/cpp/python_binding/")
sys.path.append("../../../CPlantBox/src/python_modules")

from rosi_richards_cyl import RichardsCylFoam  # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part

import matplotlib.pyplot as plt
import numpy as np
import os

""" 
Cylindrical 1D model, diffusion only (DuMux), Michaelis Menten

everything scripted, no input file needed, also works parallel with mpiexec
"""

s = RichardsWrapper(RichardsCylFoam())
s.initialize()

#alpha, theta_s, theta_r, n, Ks?
loam = [0.04, 0.43, 0.045, 1.6, 50]

#s.createGrid([0.02], [0.6], [500])  # [cm]
# with open("Hp_results.txt") as f:
    # list2 = [row.split(',')[0] for row in f]
# list2 = list2[1:]
# list2 = [float(ll) for ll in list2]

# s.createGrid1d(np.array(list2) + 0.02)


s.createGrid([0.02], [0.6], [500]) 
s.setHomogeneousIC(-100.)  # cm pressure head
s.setOuterBC("noflux")  #  [cm/day]
s.setInnerBC("constantFlux", -0.5)  #  [cm/day]

s.setVGParameters([loam])
s.initializeProblem()
s.setCriticalPressure(-15000)  # cm pressure head


fig, (ax1, ax2) = plt.subplots(1, 2)

times = [0., 0.5/24., 1./24.]  # days
s.ddt = 1.e-5

col = ["r*", "b*", "g*", "c*", "m*", "y*", ]

x = np.array(s.getSolutionHead())
print(x.flatten())

for i, dt in enumerate(np.diff(times)):

    print("*****", "external time step", dt, " d, simulation time", s.simTime, "d, internal time step", s.ddt, "d")

    s.solve(dt)

    points = s.getDofCoordinates()

    x = np.array(s.getSolutionHead())
    print(points.flatten())
    print(x.flatten())


