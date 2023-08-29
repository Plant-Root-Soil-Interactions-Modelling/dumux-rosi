import sys; sys.path.append("../modules/"); sys.path.append("../../../CPlantBox/");  sys.path.append("../../build-cmake/cpp/python_binding/")
sys.path.append("../../../CPlantBox/src/python_modules")


from rosi_richards10c_cyl import Richards10CCylFoam  # C++ part (Dumux binding)

from richards import RichardsWrapper  # Python part

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()

def write_file_array(name, data, space =","):
    name2 = './results/'+ name+ '.txt'
    with open(name2, 'a') as log:
        log.write(space.join([num for num in map(str, data)])  +'\n')
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

s = RichardsWrapper(Richards10CCylFoam())

s.initialize()

# theta_r, theta_s, alpha, n, Ks
loam = [0.045, 0.43, 0.04, 1.6, 50]

nCells = 500
s.createGrid([0.02], [0.6], [nCells])  # [cm]
s.setParameter( "Soil.Grid.Cells", str(nCells))

s.setHomogeneousIC(-100.)  # cm pressure head
s.setOuterBC("constantFluxCyl", 0)  #  [cm/day]
s.setInnerBC("constantFluxCyl", -0.26)  #  [cm/day]

#s.setICZ_solute(0.)  # [kg/m2] 

molarMassWat = 18. # [g/mol]
densityWat = 1. #[g/cm3]
# [mol/cm3] = [g/cm3] /  [g/mol] 
molarDensityWat =  densityWat / molarMassWat # [mol/cm3] 

solidDensity = 2700 # [kg/m^3 solid]
solidMolarMass = 60.08e-3 # [kg/mol] 
# [mol / m3 solid] =[kg/m^3 solid] / [kg/mol] 
solidMolDensity = solidDensity/solidMolarMass
# [mol / m3 scv] = [mol / m3 solid] * [m3 solid /m3 space]
bulkDensity_m3 = solidMolDensity*(1.-0.43)

MolarMass = 1.8e-2 #[kg/mol] 0.02003 #[kg/mol]
exud = 1. # [mol/cm2/d]#1.0* MolarMass *1000# [mol/cm2/d] * [kg/mol] * [g/kg] =  [g/cm2/d]
Ds = 1e-8 # m^2/s
Dl = 1e-9

numComp = 8
numFluidComp = 2


s.setParameter( "Soil.MolarMass", str(solidMolarMass))
s.setParameter( "Soil.solidDensity", str(solidDensity))

s.setParameter( "Soil.BC.Bot.C1Type", str(3))
s.setParameter( "Soil.BC.Top.C1Type", str(3))
s.setParameter( "Soil.BC.Bot.C1Value", str(exud)) 
s.setParameter( "Soil.BC.Top.C1Value", str(0 )) 

s.setParameter("Component.LiquidDiffusionCoefficient", str(Ds)) #m^2/s

    
gradient = False
if gradient:
    #C_S = np.array([0.1, 0.3, 0.4, 0.5, 9])  #mol/cm3 wat
    #s.setParameter("Soil.IC.C1Z", "0.0001 0.0003 0.0004 0.0005 0.009" )  #mol/cm3 / mol/cm3 = mol/mol 
    #s.setParameter("Soil.IC.C1", str(C_S/molarDensityWat)[1:(len(str(C_S))-1)])   #mol/cm3 / mol/cm3 = mol/mol
    C_S = np.array([6, 0.2])  #mol/cm3 wat
    s.setParameter("Soil.IC.C1Z", "0.0002 0.006" )  #mol/cm3 / mol/cm3 = mol/mol 
    s.setParameter("Soil.IC.C1", str(C_S/molarDensityWat)[1:(len(str(C_S/molarDensityWat))-1)])   #mol/cm3 / mol/cm3 = mol/mol 
else:
    C_S = 0.1  #mol/cm3 wat
    s.setParameter("Soil.IC.C1", str(C_S/ molarDensityWat) )  #mol/cm3 / mol/cm3 = mol/mol 

s.setVGParameters([loam])



s.setParameter("Problem.EnableGravity", "false")
s.setParameter("Problem.verbose", "1")
s.setParameter("Flux.UpwindWeight", "0.5")
s.setParameter("Newton.EnableAbsoluteResidualCriterion", "true")
s.setParameter("Newton.MaxAbsoluteResidual", "1.e-10")
s.setParameter("Newton.EnableChop", "true")
s.setParameter("Newton.EnableResidualCriterion", "true")
s.setParameter("Newton.EnableShiftCriterion", "true")
s.setParameter("Newton.MaxAbsoluteResidual", "1e-10")

s.setParameter("Newton.MaxRelativeShift", "1e-10")

s.setParameter("Newton.MaxSteps", "30")
s.setParameter("Newton.ResidualReduction", "1e-10")
s.setParameter("Newton.SatisfyResidualAndShiftCriterion", "true")
s.setParameter("Newton.TargetSteps", "10")
s.setParameter("Newton.UseLineSearch", "false")
s.setParameter("Newton.EnablePartialReassembly", "true")
s.setParameter("Grid.Overlap", "0")  #no effec5

s.initializeProblem()
s.setParameter("Flux.UpwindWeight", "1")
s.setCriticalPressure(-15000)  # cm pressure head


fig, (ax1, ax2) = plt.subplots(1, 2)

times = [0., 5./24., 10./24.]  # days
s.ddt = 1.e-5

points = s.getDofCoordinates().flatten()

x = np.array(s.getSolutionHead())
write_file_array("pressureHead",x.flatten())
write_file_array("coord",points)
theta = np.array(s.getWaterContent()).flatten()
write_file_array("theta",theta)
write_file_array("getSaturation",np.array(s.getSaturation()).flatten())
write_file_array("krs",np.array(s.getKrw()).flatten())
# [g/cm3] * [mol/kg] * [kg/g] = [mol/cm3]

for i in range(numFluidComp):
    write_file_array("solute_conc"+str(i+1), np.array(s.getSolution_(i+1)).flatten()*molarDensityWat ) 
for i in range(numFluidComp, numComp):
    write_file_array("solute_conc"+str(i+1), np.array(s.getSolution_(i+1)).flatten()* bulkDensity_m3 /1e6 ) 


for i, dt in enumerate(np.diff(times)):

    if rank == 0:
        print("*****", "external time step", dt, " d, simulation time", s.simTime, "d, internal time step", s.ddt, "d")

    s.solve(dt)#, maxDt = 2500/(24*3600))

    x = np.array(s.getSolutionHead())
    write_file_array("pressureHead",x.flatten())
    #print(x.flatten())
    write_file_array("coord",points)
    write_file_array("theta",np.array(s.getWaterContent()).flatten())
    write_file_array("getSaturation",np.array(s.getSaturation()).flatten())
    write_file_array("krs",np.array(s.getKrw()).flatten())
    print("min konz",min(np.array(s.getSolution_(i+1)).flatten()*molarDensityWat))
    for i in range(numFluidComp):
        write_file_array("solute_conc"+str(i+1), np.array(s.getSolution_(i+1)).flatten()*molarDensityWat) 
    for i in range(numFluidComp, numComp):
        write_file_array("solute_conc"+str(i+1), np.array(s.getSolution_(i+1)).flatten()* bulkDensity_m3 /1e6 ) 




