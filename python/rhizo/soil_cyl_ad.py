import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src")

from rosi_richardsnc_cyl import RichardsNCCylFoam  # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part

import matplotlib.pyplot as plt
import numpy as np
import os
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()

""" 
Cylindrical 1D model, (DuMux solver) Advection Diffusion  

everything scripted, no input file needed, also works parallel with mpiexec

works fine at moderate resolution (100), overestimates water content in final step
grid DOF 11 logspace is good 

RuntimeError: Dune::NotImplemented [evalFlux:/home/daniel/workspace/DUMUX/dumux/dumux/assembly/cclocalresidual.hh:109]: Mixed boundary conditions. Use pure boundary conditions by converting Dirichlet BCs to Robin BCs
"""

s = RichardsWrapper(RichardsNCCylFoam())
s.initialize()

N = 200
loam = [0.045, 0.43, 0.04, 1.6, 50]
# points = np.linspace(0.02, 0.6, N)
points = np.logspace(np.log10(0.02), np.log10(0.6), N)
s.createGrid1d(points)  # [cm]

s.setHomogeneousIC(-100.)  # cm pressure head
s.setOuterBC("noflux")  #  [cm/day]
s.setInnerBC("fluxCyl", -0.1)  # [cm/day] -0.1

s.setParameter("Component.MolarMass", "1.8e-2")  # TODO no idea, where this is neeeded, i don't want to use moles ever
s.setParameter("Component.LiquidDiffusionCoefficient", "1.e-9")  # m^2 s-1
# s.setParameter("SpatialParams.Tortuosity", "0.129")

s.setParameter("Component.BufferPower", "140.")  # buffer power = \rho * Kd [1]
s.setParameter("Soil.IC.C", "0.01")  # (mol)g / cm3  # TODO specialised setter?

s.setParameter("Soil.BC.Top.SType", "1")  # michaelisMenten=8 (SType = Solute Type)
s.setParameter("Soil.BC.Top.CValue", "0.")  # michaelisMenten=8 (SType = Solute Type)
# s.setParameter("Soil.BC.Top.SType", "1")  # michaelisMenten=8 (SType = Solute Type)
# s.setParameter("Soil.BC.Top.CValue", "0.007")  # michaelisMenten=8 (SType = Solute Type)

s.setParameter("Soil.BC.Bot.SType", "1")  # michaelisMenten=8 (SType = Solute Type)
s.setParameter("Soil.BC.Bot.CValue", "0.")
# s.setParameter("TimeLoop.MaxTimeStepSize", "60")  # seconds
# s.setParameter("Soil.BC.Bot.SType", "8")  # michaelisMenten (SType = Solute Type)
# s.setParameter("RootSystem.Uptake.Vmax", s.dumux_str(3.26e-6 * 24 * 3600))  # (mol)g /cm^2 / s - > (mol)g /cm^2 / day
# s.setParameter("RootSystem.Uptake.Km", s.dumux_str(5.8e-3))  # (mol)g / cm3

s.setParameter("Newton.EnableAbsoluteResidualCriterion", "True")
s.setParameter("Newton.MaxAbsoluteResidual", "1.e-10")

s.setVGParameters([loam])
s.initializeProblem()
s.setCriticalPressure(-15000)  # cm pressure head

if rank == 0:
    print(s)  # TODO make super nice

fig, (ax1, ax2) = plt.subplots(1, 2)

days_ = np.linspace(0, 20, 49)  # COMSOL time steps
times = [0., days_[1], days_[23], days_[48]]  # days   , 25, 30
print("times", times, "days")
s.ddt = 1.e-5

col = ["r*", "g*", "b*", "c*", "m*", "y*", ]

idx = s.pick([0.02])
print("boundary element ", idx)
f = []
maxDt = 1

for i, dt in enumerate(np.diff(times)):

    if rank == 0:
        print("*****", "external time step", dt, " d, simulation time", s.simTime, "d, internal time step", s.ddt, "d")

    s.ddt = 1.e-3
    s.solve(dt, maxDt)

    points = s.getDofCoordinates()

    x = s.getSolutionHead()
    y = s.getSolution(1) * 1  # [kg/kg] -> 1/1000 [kg/m3] [] -> 1 [g/cm3] # solute concentration

    ax1.plot(points[:], x, col[i % len(col)], label = "dumux {:g} days".format(s.simTime))
    ax2.plot(points[:], y, col[i % len(col)], label = "dumux {:g} days".format(s.simTime))

    f.append(s.getNeumann(idx))

data = np.loadtxt("../../cpp/soil_richards/python/cylinder_1d_Comsol.txt", skiprows = 8)
z_comsol = data[:, 0]
ax1.plot(z_comsol + 0.02, data[:, 24], "k", label = "comsol 10 days")
ax1.plot(z_comsol + 0.02, data[:, 49], "k:", label = "comsol 20 days")
ax1.set_xlabel("distance from root axis (cm)")
ax1.set_ylabel("soil matric potential (cm)")
ax1.legend()

data = np.loadtxt("../../cpp/soil_richardsnc/python/c_results.txt", skiprows = 8)  # buffer power = 100
# data = np.loadtxt("c_results_nob.txt", skiprows=8) # buffer power = 0

z_comsol = data[:, 0]
ax2.plot(z_comsol + 0.02, data[:, 2], "r", label = "~0.42 days")  # indices = days_ indicdes +1 (radii)
ax2.plot(z_comsol + 0.02, data[:, 24], "g", label = "comsol ~9.2 days")
ax2.plot(z_comsol + 0.02, data[:, 49], "b:", label = "comsol 20 days")
ax2.set_xlabel('distance from the root axis (cm)')
ax2.set_ylabel('solute concentration (g/cm3)')
ax2.legend()

plt.show()
