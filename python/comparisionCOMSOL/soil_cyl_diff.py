import sys; sys.path.append("../modules/"); sys.path.append("../../../CPlantBox/");  sys.path.append("../../build-cmake/cpp/python_binding/")
sys.path.append("../../../CPlantBox/src/python_modules")

from rosi_richards2c_cyl import RichardsNCCylFoam  # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part

import matplotlib.pyplot as plt
import numpy as np
import os
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()

def write_file_array(name, data):
    name2 = './results2c/'+ name+ '.txt'
    with open(name2, 'a') as log:
        log.write(','.join([num for num in map(str, data)])  +'\n')

""" 
Cylindrical 1D model, diffusion only (DuMux), Michaelis Menten

everything scripted, no input file needed, also works parallel with mpiexec
"""

results_dir="./results2c/"
if not os.path.exists(results_dir):
    os.makedirs(results_dir)
else:
    test = os.listdir(results_dir)
    for item in test:
        try:
            os.remove(results_dir+item)
        except:
            pass
            
s = RichardsWrapper(RichardsNCCylFoam())
s.initialize()

loam = [0.045, 0.43, 0.04, 1.6, 50]

s.createGrid([0.02], [0.6], [500])  # [cm]

s.setHomogeneousIC(-100.)  # cm pressure head
s.setOuterBC("noflux")  #  [cm/day]
s.setInnerBC("constantFluxCyl", -0.26) #  [cm/day]

s.setParameter("Soil.BC.Bot.SType", "8")  # michaelisMenten (SType = Solute Type)
s.setParameter("Component.MolarMass", "1.8e-2")  # TODO no idea, where this is neeeded, i don't want to use moles ever
s.setParameter("Component.LiquidDiffusionCoefficient", "1.e-9")

s.setParameter("Soil.IC.C", "0.01")  # g / cm3  # TODO specialised setter?
s.setParameter("0.Component.BufferPower", "140")  
s.setParameter("1.Component.BufferPower", "140")  
#s.setParameter("RootSystem.Uptake.Vmax", s.dumux_str(3.26e-6 * 24 * 3600 * 1.e4))  # g /cm^2 / s - > g / m^2 / d
#s.setParameter("RootSystem.Uptake.Km", s.dumux_str(5.8e-3 * 1.e4))  # g / cm3 todo setter

s.setVGParameters([loam])
s.initializeProblem()
s.setCriticalPressure(-15000)  # cm pressure head

if rank == 0:
    print(s)  # TODO make super nice

fig, (ax1, ax2) = plt.subplots(1, 2)

times = [0., 10., 20.]  # days
s.ddt = 1.e-5

molarMassWat = 18. # [g/mol]
densityWat = 1. #[g/cm3]
# [mol/cm3] = [g/cm3] /  [g/mol] 
molarDensityWat =  densityWat / molarMassWat # [mol/cm3] 

points = s.getDofCoordinates()
x = np.array(s.getSolutionHead())
write_file_array("pressureHead",x.flatten())
write_file_array("coord",points.flatten())

write_file_array("solute_conc1", np.array(s.getSolution_(1)).flatten()*molarDensityWat) 
                                                     

for i, dt in enumerate(np.diff(times)):

    if rank == 0:
        print("*****", "external time step", dt, " d, simulation time", s.simTime, "d, internal time step", s.ddt, "d")

    s.solve(dt)

    points = s.getDofCoordinates()
    x = np.array(s.getSolutionHead())
    
    write_file_array("pressureHead",x.flatten())
    write_file_array("coord",points.flatten())

    write_file_array("solute_conc1", np.array(s.getSolution_(1)).flatten()*molarDensityWat) 
