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

s.createGrid([0.02], [0.6], [300])  # [cm]

s.setHomogeneousIC(-100.)  # cm pressure head
s.setOuterBC("constantFluxCyl", 0)  #  [cm/day]
s.setInnerBC("constantFluxCyl", 0)  #  [cm/day]

#s.setICZ_solute(0.)  # [kg/m2] 
s.setParameter( "Soil.BC.Bot.C1Type", str(3))
s.setParameter( "Soil.BC.Bot.C2Type", str(3))
s.setParameter( "Soil.BC.Top.C1Type", str(3))
s.setParameter( "Soil.BC.Top.C2Type", str(3))
MolarMass = 1.8e-2 #[kg/mol]
exud = 1. # [g/cm2/d]#1.0* MolarMass *1000# [mol/cm2/d] * [kg/mol] * [g/kg] =  [g/cm2/d]
s.setParameter( "Soil.BC.Bot.C1Value", str(exud )) 
s.setParameter( "Soil.BC.Bot.C2Value", str(0 )) 
s.setParameter( "Soil.BC.Top.C1Value", str(0 )) 
s.setParameter( "Soil.BC.Top.C2Value", str(0 ))  

s.setParameter("1.Component.MolarMass", str(MolarMass)) 
s.setParameter("1.Component.LiquidDiffusionCoefficient", "1.e-9") #m^2/s
s.setParameter("2.Component.MolarMass", str(MolarMass)) 
s.setParameter("2.Component.LiquidDiffusionCoefficient", "1.e-9") #m^2/s

ci = 1
#s.setParameter("Soil.IC.C1Z", "0.0002 0.0003 0.0004 0.0005 0.006") 
#s.setParameter("Soil.IC.C1", "0.02 0.03 0.04 0.05 0.6") 
s.setParameter("Component.BufferPower", "0")  
           

s.setVGParameters([loam])
s.initializeProblem()
s.setCriticalPressure(-15000)  # cm pressure head


s.setParameter("Newton.EnableAbsoluteResidualCriterion", "True")
s.setParameter("Newton.MaxAbsoluteResidual", "1.e-10")
s.setParameter("Newton.EnableChop", "false")
s.setParameter("Newton.EnableResidualCriterion", "true")
s.setParameter("Newton.EnableShiftCriterion", "true")
s.setParameter("Newton.MaxAbsoluteResidual", "1e-10")
s.setParameter("Newton.MaxRelativeShift", "1e-10")
s.setParameter("Newton.MaxSteps", "30")
s.setParameter("Newton.ResidualReduction", "1e-5")
s.setParameter("Newton.SatisfyResidualAndShiftCriterion", "true")
s.setParameter("Newton.TargetSteps", "10")
s.setParameter("Newton.UseLineSearch", "false")
s.setParameter("Newton.EnablePartialReassembly", "false")
s.setParameter("Grid.Overlap", "0")  #no effec5

if rank == 0:
    print(s)  # TODO make super nice

fig, (ax1, ax2) = plt.subplots(1, 2)

times = [0., 5./24, 10./24.]  # days
s.ddt = 1.e-5

col = ["r*", "b*", "g*", "c*", "m*", "y*", ]
points = s.getDofCoordinates()

x = np.array(s.getSolutionHead())
write_file_array("pressureHead",x.flatten())
write_file_array("coord",points.flatten())
write_file_array("theta",np.array(s.getWaterContent()).flatten())
write_file_array("getSaturation",np.array(s.getSaturation()).flatten())
write_file_array("krs",np.array(s.getKrw()).flatten())
# [g/cm3] * [mol/kg] * [kg/g] = [mol/cm3]
write_file_array("solute_conc", np.array(s.getSolution_(1)).flatten()) 

for i, dt in enumerate(np.diff(times)):

    if rank == 0:
        print("*****", "external time step", dt, " d, simulation time", s.simTime, "d, internal time step", s.ddt, "d")

    s.solve(dt)

    x = np.array(s.getSolutionHead())
    write_file_array("pressureHead",x.flatten())
    write_file_array("coord",points.flatten())
    write_file_array("theta",np.array(s.getWaterContent()).flatten())
    write_file_array("getSaturation",np.array(s.getSaturation()).flatten())
    write_file_array("krs",np.array(s.getKrw()).flatten())
    write_file_array("solute_conc", np.array(s.getSolution_(1)).flatten()) 

