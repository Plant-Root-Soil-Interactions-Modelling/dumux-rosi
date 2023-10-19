import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src")

from rosi_richardsnc_cyl import RichardsNCCylFoam  # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part

import matplotlib.pyplot as plt
import numpy as np
import os
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()

""" 
Cylindrical 1D model, diffusion only (DuMux), Michaelis Menten

everything scripted, no input file needed, also works parallel with mpiexec
"""

s = RichardsWrapper(RichardsNCCylFoam())
s.initialize()

loam = [0.045, 0.43, 0.04, 1.6, 50]

s.createGrid([0.02], [0.6], [500])  # [cm]

s.setHomogeneousIC(-100.)  # cm pressure head
s.setOuterBC("noflux")  #  [cm/day]
s.setInnerBC("noflux")  #  [cm/day]

s.setParameter("Soil.BC.Bot.SType", "8")  # michaelisMenten (SType = Solute Type)
s.setParameter("Component.MolarMass", "1.8e-2")  # TODO no idea, where this is neeeded, i don't want to use moles ever
s.setParameter("Component.LiquidDiffusionCoefficient", "1.e-9")

s.setParameter("Soil.IC.C", "0.01")  # g / cm3  # TODO specialised setter?
s.setParameter("Component.BufferPower", "140")  # buffer power = \rho * Kd [1]
s.setParameter("RootSystem.Uptake.Vmax", s.dumux_str(3.26e-6 * 24 * 3600 * 1.e4))  # g /cm^2 / s - > g / m^2 / d
s.setParameter("RootSystem.Uptake.Km", s.dumux_str(5.8e-3 * 1.e4))  # g / cm3 todo setter

s.setVGParameters([loam])
s.initializeProblem()
s.setCriticalPressure(-15000)  # cm pressure head

if rank == 0:
    print(s)  # TODO make super nice

fig, (ax1, ax2) = plt.subplots(1, 2)

times = [0., 10., 20.]  # days
s.ddt = 1.e-5

col = ["r*", "b*", "g*", "c*", "m*", "y*", ]

for i, dt in enumerate(np.diff(times)):

    if rank == 0:
        print("*****", "external time step", dt, " d, simulation time", s.simTime, "d, internal time step", s.ddt, "d")

    s.solve(dt)

    points = s.getDofCoordinates()

    x = s.getSolutionHead()
    y = s.getSolution(1)  # solute concentration

    ax1.plot(points[:], x, col[i % len(col)], label = "dumux {} days".format(s.simTime))
    ax2.plot(points[:], y, col[i % len(col)], label = "dumux {} days".format(s.simTime))

# os.chdir("../../../build-cmake/rosi_benchmarking/soil_richards/python")
# data = np.loadtxt("cylinder_1d_Comsol_water.txt", skiprows=8)
# z_comsol = data[:, 0]
# ax1.plot(z_comsol + 0.02, data[:, 25], "k", label="comsol 10 days")
# ax1.plot(z_comsol + 0.02, data[:, -1], "k:", label="comsol 20 days")
# ax1.set_xlabel("distance from root axis (cm)")
# ax1.set_ylabel("soil matric potential (cm)")
# ax1.legend()

data = np.loadtxt("comsol_c_diff.txt", skiprows = 8)
z_comsol = data[:, 0]
ax2.plot(z_comsol + 0.02, data[:, 25], "k", label = "comsol 10 days")
ax2.plot(z_comsol + 0.02, data[:, -1], "k:", label = "comsol 20 days")
ax2.set_xlabel('distance from the root axis (cm)')
ax2.set_ylabel('solute concentration (g/cm3)')

ax2.legend()

plt.show()
