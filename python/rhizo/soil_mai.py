import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src")

from rosi_richardsnc_cyl import RichardsNCCylFoam  # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part

import matplotlib.pyplot as plt
import numpy as np
import os
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()

""" 
Cylindrical 1D model from Bauw et al. 2020

REVISE
RuntimeError: Dune::NotImplemented [evalFlux:/home/daniel/workspace/DUMUX/dumux/dumux/assembly/cclocalresidual.hh:109]: Mixed boundary conditions. Use pure boundary conditions by converting Dirichlet BCs to Robin BCs

"""

SMALL_SIZE = 22
MEDIUM_SIZE = 22
BIGGER_SIZE = 22
plt.rc('font', size = SMALL_SIZE)  # controls default text sizes
plt.rc('axes', titlesize = SMALL_SIZE)  # fontsize of the axes title
plt.rc('axes', labelsize = MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize = SMALL_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize = SMALL_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize = SMALL_SIZE)  # legend fontsize
plt.rc('figure', titlesize = BIGGER_SIZE)  # fontsize of the figure title

s = RichardsWrapper(RichardsNCCylFoam())
s.initialize()

N = 100

# loam = [0.045, 0.43, 0.04, 1.6, 50]
loam = [0.03, 0.345, 0.01, 2.5, 28.6]  #  qr, qs, alpha, n, ks

points = np.linspace(0.02, 0.6, N)
# points = np.logspace(np.log10(0.02), np.log10(0.6), N)
s.createGrid1d(points)  # [cm]

s.setHomogeneousIC(-100.)  # cm pressure head
s.setOuterBC("noflux")  #  [cm/day]
s.setInnerBC("fluxCyl", -0.1)  # [cm/day] -0.1

# from publication branch input file... (for P)
# MolarMass = 31e-3 #[kg/mol]
# liquidDiffCoeff = 0.6e-9 # [m2 s-1] # http://www.aqion.de/site/194
s.setParameter("Component.MolarMass", "1.801528e-2")  # TODO no idea, where this is neeeded, i don't want to use moles ever
s.setParameter("Component.LiquidDiffusionCoefficient", "6.e-10")  # m^2 s-1

s.setParameter("Component.FreundlichN", ".1")
s.setParameter("Component.FreundlichK", "0.0001")  # 124 -> 0.0311 [g^0.6 * g^-0.4 * cm^3^0.4]
# s.setParameter("Component.BufferPower", "600")  # buffer power = \rho * Kd [1]

s.setParameter("Soil.IC.C", "0.02")  # (mol)g / cm3  # TODO specialised setter?

s.setParameter("Soil.BC.Top.SType", "2")  # michaelisMenten=8 (SType = Solute Type)
s.setParameter("Soil.BC.Top.CValue", "0.")  # michaelisMenten=8 (SType = Solute Type)
# s.setParameter("Soil.BC.Top.SType", "1")  # michaelisMenten=8 (SType = Solute Type)
# s.setParameter("Soil.BC.Top.CValue", "0.007")  # michaelisMenten=8 (SType = Solute Type)

s.setParameter("Soil.BC.Bot.SType", "1")  # michaelisMenten=8 (SType = Solute Type)
s.setParameter("Soil.BC.Bot.CValue", "0.")
# s.setParameter("Soil.BC.Bot.SType", "8")  # michaelisMenten (SType = Solute Type)
# s.setParameter("RootSystem.Uptake.Vmax", s.dumux_str(0.01 * 3.26e-6 * 24 * 3600))  # (mol)g /cm^2 / s - > (mol)g /cm^2 / day
# s.setParameter("RootSystem.Uptake.Km", s.dumux_str(5.8e-3))  # (mol)g / cm3

s.setParameter("TimeLoop.MaxTimeStepSize", "60")  # seconds
s.setParameter("Newton.EnableAbsoluteResidualCriterion", "True")
s.setParameter("Newton.MaxAbsoluteResidual", "1.e-10")

s.setVGParameters([loam])
s.initializeProblem()
s.setCriticalPressure(-15000)  # cm pressure head

if rank == 0:
    print(s)  # TODO make super nice

fig, (ax1, ax2) = plt.subplots(1, 2)

days_ = np.linspace(0, 20, 49)  # COMSOL time steps
times = [0., 0.1, 1, 2, 5, days_[24], days_[48]]  # days   , 25, 30 days_[1],
print("times", times, "days")
s.ddt = 1.e-5

col = ["r*", "g*", "b*", "c*", "m*", "y*", ]

idx = s.pick([0.02])
print("boundary element ", idx)
f = []
maxDt = 0.1

for i, dt in enumerate(np.diff(times)):

    if rank == 0:
        print("*****", "external time step", dt, " d, simulation time", s.simTime, "d, internal time step", s.ddt, "d")

    s.ddt = 1.e-3
    s.solve(dt)

    points = s.getDofCoordinates()

    x = s.getSolutionHead()
    y = s.getSolution(1) * 1.  # [kg/kg] -> 1/1000 [kg/m3] [] -> 1 [g/cm3] # solute concentration

    ax1.plot(points[:], x, col[i % len(col)], label = "dumux {:g} days".format(s.simTime))
    ax2.plot(points[:], y, col[i % len(col)], label = "dumux {:g} days".format(s.simTime))

    f.append(s.getNeumann(idx))

ax1.set_xlabel("distance from root axis (cm)")
ax1.set_ylabel("soil matric potential (cm)")
ax1.legend()
ax2.set_xlabel('distance from the root axis (cm)')
ax2.set_ylabel('solute concentration (g/cm3)')
ax2.legend()

# data = np.loadtxt("bauw2020_pressure.txt", skiprows=8)
data = np.loadtxt("pressure.txt", skiprows = 8)
print(data.shape)
z_comsol = data[:, 0]
ax1.plot(z_comsol + 0.02, data[:, 25], "k--")
ax1.plot(z_comsol + 0.02, data[:, 49], "k--")
ax1.set_xlabel("distance from root axis (cm)")
ax1.set_ylabel("soil matric potential (cm)")
ax1.legend()

# data = np.loadtxt("bauw2020_concentration.txt", skiprows=8)  # buffer power = 100
data = np.loadtxt("concentration.txt", skiprows = 8)  # buffer power = 100
z_comsol = data[:, 0]
# ax2.plot(z_comsol + 0.02, data[:, 25], "k--")  # indices = days_ indicdes +1 (radii)
# ax2.plot(z_comsol + 0.02, data[:, 49], "k--")
ax2.set_xlabel('distance from the root axis (cm)')
ax2.set_ylabel('solute concentration (g/cm3)')
ax2.legend()

plt.show()
