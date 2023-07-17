import sys; sys.path.append("../modules/"); sys.path.append("../../../CPlantBox/");  sys.path.append("../../build-cmake/cpp/python_binding/")
sys.path.append("../../../CPlantBox/src/python_modules")

from rosi_richards3c_cyl import Richards3CCylFoam  # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part

import matplotlib.pyplot as plt
import numpy as np
import os
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()

def write_file_array(name, data):
    name2 = './results/'+ name+ '.txt'
    with open(name2, 'a') as log:
        log.write(','.join([num for num in map(str, data)])  +'\n')
""" 
Cylindrical 1D model, diffusion only (DuMux), Michaelis Menten

everything scripted, no input file needed, also works parallel with mpiexec
"""

results_dir="./results/"
if not os.path.exists(results_dir):
    os.makedirs(results_dir)
else:
    test = os.listdir(results_dir)
    for item in test:
        try:
            os.remove(results_dir+item)
        except:
            pass


s = RichardsWrapper(Richards3CCylFoam())
s.initialize()

#alpha, theta_s, theta_r, n, Ks
loam = [0.045, 0.43, 0.04, 1.6, 50]

s.createGrid([0.02], [0.6], [500])  # [cm]

s.setHomogeneousIC(-100.)  # cm pressure head
s.setOuterBC("constantFluxCyl", 0)  #  [cm/day]
s.setInnerBC("constantFluxCyl", -1)  #  [cm/day]

s.setICZ_solute(0.)  # [kg/m2] 
s.setParameter( "BC.Bot.C1Type", str(3))
s.setParameter( "BC.Bot.C2Type", str(3))
s.setParameter( "BC.Top.C1Type", str(3))
s.setParameter( "BC.Top.C2Type", str(3))
MolarMass = 1.8e-2
exud = 0.01* MolarMass # [mol/cm2/d] * [kg/mol] =  [kg/m2/s]
s.setParameter( "BC.Bot.C1Value", str(exud )) 
s.setParameter( "BC.Bot.C2Value", str(0 )) 
s.setParameter( "BC.Top.C1Value", str(0 )) 
s.setParameter( "BC.Top.C2Value", str(0 ))  

s.setParameter("1.Component.MolarMass", str(MolarMass)) 
s.setParameter("1.Component.LiquidDiffusionCoefficient", "1.e-9") #m^2/s
s.setParameter("2.Component.MolarMass", str(MolarMass)) 
s.setParameter("2.Component.LiquidDiffusionCoefficient", "1.e-9") #m^2/s

s.setParameter("Soil.IC.C", "0.0")  # g / cm3  # TODO specialised setter?
s.setParameter("Component.BufferPower", "0")  # buffer power = \rho * Kd [1]
# s.setParameter("RootSystem.Uptake.Vmax", s.dumux_str(3.26e-6 * 24 * 3600 * 1.e4))  # g /cm^2 / s - > g / m^2 / d
# s.setParameter("RootSystem.Uptake.Km", s.dumux_str(5.8e-3 * 1.e4))  # g / cm3 todo setter

s.setVGParameters([loam])
s.initializeProblem()
s.setCriticalPressure(-15000)  # cm pressure head

if rank == 0:
    print(s)  # TODO make super nice

fig, (ax1, ax2) = plt.subplots(1, 2)

times = [0., 10., 20.]  # days
s.ddt = 1.e-5

col = ["r*", "b*", "g*", "c*", "m*", "y*", ]

x = np.array(s.getSolutionHead())
write_file_array("pressureHead",x.flatten())
#write_file_array("coord",points.flatten())
write_file_array("theta",np.array(s.getWaterContent()).flatten())
write_file_array("getSaturation",np.array(s.getSaturation()).flatten())
write_file_array("krs",np.array(s.getKrw()).flatten())

for i, dt in enumerate(np.diff(times)):

    if rank == 0:
        print("*****", "external time step", dt, " d, simulation time", s.simTime, "d, internal time step", s.ddt, "d")

    s.solve(dt)

    x = np.array(s.getSolutionHead())
    write_file_array("pressureHead",x.flatten())
    #write_file_array("coord",points.flatten())
    write_file_array("theta",np.array(s.getWaterContent()).flatten())
    write_file_array("getSaturation",np.array(s.getSaturation()).flatten())
    write_file_array("krs",np.array(s.getKrw()).flatten())

