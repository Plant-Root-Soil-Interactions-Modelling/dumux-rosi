import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src")

from rosi_richards_cyl import RichardsCylFoam  # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()  # MPI

import matplotlib.pyplot as plt
import numpy as np
import os
import time

import functional.van_genuchten as vg

""" 
Cylindrical 1D model (DuMux) Water movement only 

DuMux everything scripted, no input file needed, also works parallel with mpiexec
"""

s = RichardsWrapper(RichardsCylFoam())
s.initialize()

# loam = [0.045, 0.43, 0.04, 1.6, 50]
loam = [0.03, 0.345, 0.01, 2.5, 28.6]  #  qr, qs, alpha, n, ks

N = 100
s.createGrid([0.02], [0.6], [N])  # [cm]

s.setHomogeneousIC(-100.)  # cm pressure head
s.setOuterBC("fluxCyl", 0.)  #  [cm/day]
s.setInnerBC("fluxCyl", -0.1)  # [cm/day]
s.setVGParameters([loam])
s.initializeProblem()
s.setCriticalPressure(-15000)  # cm pressure head

if rank == 0:
    print(s)

times = [0., 10, 20]  # days
# times = np.linspace(0, 20, 200)
s.ddt = 1.e-5

col = ["r*", "b*", "g*", "c*", "m*", "y*", ]

t = time.time()
v = []

points = s.getDofCoordinates()
inner = points - np.ones(points.shape) * ((0.6 - 0.02) / N / 2)
outer = points + np.ones(points.shape) * ((0.6 - 0.02) / N / 2)
areas = np.pi * (np.multiply(outer, outer) - np.multiply(inner, inner))
print("grid volume", np.sum(areas))
water_vol = np.multiply(s.getWaterContent(), areas)
v.append(np.sum(water_vol))

# print(((0.6 - 0.02) / (N + 1) / 2))
# print(points[0:5])
# print(inner[0:5])
# # print(outer)
# input()

for i, dt in enumerate(np.diff(times)):

    if rank == 0:
        print("*****", "external time step", dt, " d, simulation time", s.simTime, "d, internal time step", s.ddt, "d")

    s.solve(dt)

    x = s.getSolutionHead()
    plt.plot(points[:], x, col[i % len(col)], label = "dumux {} days".format(s.simTime))

    water_vol = np.multiply(s.getWaterContent(), areas)
    v.append(np.sum(water_vol))

print("elapsed time", time.time() - t)

data = np.loadtxt("bauw2020_pressure.txt", skiprows = 8)
z_comsol = data[:, 0]
plt.plot(z_comsol + 0.02, data[:, 25], "k", label = "comsol 10 days")
plt.plot(z_comsol + 0.02, data[:, -1], "k:", label = "comsol 20 days")
plt.xlabel("distance from root axis (cm)")
plt.ylabel("soil matric potential (cm)")
plt.legend()
plt.show()

plt.plot(times, v)
plt.xlabel("time (days)")
plt.ylabel("water volume cm3")

# soil = vg.Parameters(loam)
# theta = vg.water_content(-100., soil)
# v_ = np.pi * (0.6 * 0.6 - 0.02 * 0.02)  # cm2
# print("volume", v_, "cm3")
# print("water content", theta)
# print("initial water volume", v_ * theta, "cm3 = ", v[0], "cm3")
# vol10 = v[0] - 0.1 * 2 * 0.02 * np.pi * 10
# plt.plot([10], [vol10], "r*")
# plt.show()

