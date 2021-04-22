import sys; sys.path.append("../modules/"); sys.path.append("../../../CPlantBox/");  sys.path.append("../../build-cmake/cpp/python_binding/")
sys.path.append("../../../CPlantBox/src/python_modules")

from rosi_richardsnc_cyl import RichardsNCCylFoam  # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part

import matplotlib.pyplot as plt
import numpy as np
import os
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()

""" 
Cylindrical 1D model from Bauw et al. 2020
"""

N = 1000
# loam = [0.045, 0.43, 0.04, 1.6, 50]
loam = [0.03, 0.345, 0.01, 2.5, 28.6]  #  qr, qs, alpha, n, ks

buffer_powers = [0, 200, 600, 1200]
fluxes = [-0.1, -0.01]
times = [0, 1, 10, 20, 30]  

s = [[ [] for _ in range(len(buffer_powers)) ], [[] for _ in range(len(buffer_powers))]]

for i, b in enumerate(buffer_powers):
    for j, f in enumerate(fluxes):
        
        s[j][i] = RichardsWrapper(RichardsNCCylFoam())
        s[j][i].initialize()
        
        points = np.linspace(0.02, 0.6, N)
        # points = np.logspace(np.log10(0.02), np.log10(0.6), N)
        s[j][i].createGrid1d(points)  # [cm]
        
        s[j][i].setHomogeneousIC(-100.)  # cm pressure head
        s[j][i].setParameter("Soil.IC.C", "0.02")  # (mol)g / cm3  # TODO specialised setter?
         
        s[j][i].setParameter("Component.MolarMass", "1.801528e-2")  # TODO no idea, where this is neeeded, i don't want to use moles ever
        s[j][i].setParameter("Component.LiquidDiffusionCoefficient", "6.e-10")  # m^2 s-1
        s[j][i].setParameter("Component.BufferPower", str(buffer_powers[i]))  # buffer power = \rho * Kd [1]

        s[j][i].setOuterBC("noflux")  #  [cm/day]
        s[j][i].setInnerBC("fluxCyl", str(fluxes[j]))  # [cm/day] -0.1         
        s[j][i].setParameter("Soil.BC.Top.SType", "2")  # michaelisMenten=8 (SType = Solute Type)
        s[j][i].setParameter("Soil.BC.Top.CValue", "0.")  # michaelisMenten=8 (SType = Solute Type)
        s[j][i].setParameter("Soil.BC.Bot.SType", "1")  # michaelisMenten=8 (SType = Solute Type)
        s[j][i].setParameter("Soil.BC.Bot.CValue", "0.")
        
        s[j][i].setParameter("TimeLoop.MaxTimeStepSize", "60")  # seconds
        s[j][i].setParameter("Newton.EnableAbsoluteResidualCriterion", "True")
        s[j][i].setParameter("Newton.MaxAbsoluteResidual", "1.e-10")
         
        s[j][i].setVGParameters([loam])
        s[j][i].initializeProblem()
        s[j][i].setCriticalPressure(-15000)  # cm pressure head
        
        points = s[j][i].getDofCoordinates()
        
fig, axes = plt.subplots(len(buffer_powers), len(fluxes))

col = ["r", "g:", "b--", "c", "m:", "y--", ]

maxDt = 0.1

for i, b in enumerate(buffer_powers):
    for j, f in enumerate(fluxes):

        for k, dt in enumerate(np.diff(times)):
        
            print("Calculating ", i, j, k)
        
            s[j][i].ddt = 1.e-3
            s[j][i].solve(dt, maxDt)
        
            x = s[j][i].getSolutionHead()
            y = s[j][i].getSolution(1) * 1.  # [kg/kg] -> 1/1000 [kg/m3] [] -> 1 [g/cm3] # solute concentration
            
            axes[i, j].plot(points[:], y, col[k % len(col)], label="dumux {:g} days".format(s[j][i].simTime))

            if i == len(buffer_powers) - 1:
                axes[i, j].set_xlabel('distance from the root axis (cm)')
            axes[i, j].set_ylabel('solute concentration (g/cm3)')
            axes[i, j].set_ylim([0, 0.021])
            axes[i, j].set_xlim([0.02, 0.3])
            if i == 0:
                axes[i, j].legend()
                
            axes[i, j].set_title("buffer power {:g}, root influx {:g} cm/s".format(buffer_powers[i], fluxes[j]))

# plt.tight_layout()
plt.show()
