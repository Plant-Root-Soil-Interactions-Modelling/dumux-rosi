import sys; sys.path.append("../modules/"); sys.path.append("../../../CPlantBox/");  sys.path.append("../../build-cmake/cpp/python_binding/")
sys.path.append("../../../CPlantBox/src/python_modules")


useC3 = True
if useC3:
    from rosi_richards3c_cyl import Richards3CCylFoam  # C++ part (Dumux binding)
else:
    from rosi_richardsnc_cyl import RichardsNCCylFoam  # C++ part (Dumux binding)

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

if useC3:
    s = RichardsWrapper(Richards3CCylFoam())
else:
    s = RichardsWrapper(RichardsNCCylFoam())
    
s.initialize()

# theta_r, theta_s, alpha, n, Ks
loam = [0.045, 0.43, 0.04, 1.6, 50]

s.createGrid([0.02], [0.6], [500])  # [cm]

s.setHomogeneousIC(-100.)  # cm pressure head
s.setOuterBC("constantFluxCyl", 0)  #  [cm/day]
s.setInnerBC("constantFluxCyl", 0)  #  [cm/day]

#s.setICZ_solute(0.)  # [kg/m2] 

MolarMass_kg = 5e-2 #[kg/mol]
MolarMass_g = MolarMass_kg * 1000 #[kg/mol]
exud_mol = 0.1 # [mol/cm2/d] 
exud_g = exud_mol * MolarMass_g #[g/cm2/d]

if useC3:
    s.setParameter( "Soil.BC.Bot.C1Type", str(3))
    s.setParameter( "Soil.BC.Top.C1Type", str(3))
    s.setParameter( "Soil.BC.Bot.C1Value", str(exud_g )) #g/cm^2/d
    s.setParameter( "Soil.BC.Top.C1Value", str(0 )) 
    s.setParameter("1.Component.MolarMass", str(MolarMass_kg)) 
    s.setParameter("1.Component.LiquidDiffusionCoefficient", "1.e-9") #m^2/s

    s.setParameter( "Soil.BC.Top.C2Type", str(3))
    s.setParameter( "Soil.BC.Bot.C2Type", str(3))
    s.setParameter( "Soil.BC.Bot.C2Value", str(0 )) 
    s.setParameter( "Soil.BC.Top.C2Value", str(0 )) 
    s.setParameter("2.Component.MolarMass", str(MolarMass_kg)) 
    s.setParameter("2.Component.LiquidDiffusionCoefficient", "1.e-7") #m^2/s
    
    #s.setParameter("Soil.IC.C1Z", "0.0002 0.0003 0.0004 0.0005 0.006") 
    
    #s.setParameter("Soil.IC.C1", "0.2 0.3 0.4 0.5 6") 
    #s.setParameter("Soil.IC.C1", "0.02 0.03 0.04 0.05 0.6") 
    #s.setParameter("Soil.IC.C1", "0.6 0.05 0.04 0.03 0.02") 
else:
    s.setParameter( "Soil.BC.Bot.SType", str(3))
    s.setParameter( "Soil.BC.Top.SType", str(3))
    s.setParameter( "Soil.BC.Bot.CValue", str(exud_g )) 
    s.setParameter( "Soil.BC.Top.CValue", str(0 )) 
    s.setParameter("Component.MolarMass", str(MolarMass_kg)) 
    s.setParameter("Component.LiquidDiffusionCoefficient", "1.e-9") #m^2/s
    s.setParameter("Soil.IC.CZ", "0.0002 0.0003 0.0004 0.0005 0.006") 
    s.setParameter("Soil.IC.C", "0.02 0.03 0.04 0.05 0.6") 
    
    
 



ci = 1

s.setParameter("Component.BufferPower", "0")  
           

s.setVGParameters([loam])
s.setParameter("Problem.EnableGravity", "false")
s.setParameter("Flux.UpwindWeight", "1")
s.setParameter("Newton.EnableAbsoluteResidualCriterion", "true")
s.setParameter("Newton.EnableChop", "false")
s.setParameter("Newton.EnableResidualCriterion", "true")
s.setParameter("Newton.EnableShiftCriterion", "true")
s.setParameter("Newton.MaxAbsoluteResidual", "1e-8")
s.setParameter("Newton.MaxRelativeShift", "1e-8")
s.setParameter("Newton.MaxSteps", "30")
s.setParameter("Newton.ResidualReduction", "1e-12")
s.setParameter("Newton.SatisfyResidualAndShiftCriterion", "true")
s.setParameter("Newton.TargetSteps", "10")
s.setParameter("Newton.UseLineSearch", "true")
s.setParameter("Newton.EnablePartialReassembly", "true")
s.setParameter("Grid.Overlap", "0")  #no effec5
s.initializeProblem()
s.setCriticalPressure(-15000)  # cm pressure head


if rank == 0:
    print(s)  # TODO make super nice


times = [0., 1./24, 2./24.]  # days
s.ddt = 1.e-5

points = s.getDofCoordinates()

x = np.array(s.getSolutionHead())
write_file_array("pressureHead",x.flatten())
write_file_array("coord",points.flatten())
write_file_array("theta",np.array(s.getWaterContent()).flatten())
write_file_array("getSaturation",np.array(s.getSaturation()).flatten())
write_file_array("krs",np.array(s.getKrw()).flatten())
# [g/cm3] * [mol/kg] * [kg/g] = [mol/cm3]
solute_mC = s.getSolution_(1).flatten() # [g/cm3 solution]
solute_mol = solute_mC / MolarMass_g # [g/cm3] * [mol/g] = [mol/cm3]
write_file_array("solute_conc_mol", solute_mol ) #[mol/cm3 solution]
write_file_array("solute_conc_g", solute_mC ) #[g/cm3 solution] == [g/g water]
#raise Exception
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
    
    solute_mC = s.getSolution_(1).flatten() # [g/g solution]
    avgRho = rhowat
    solute_mol = solute_mC / MolarMass_g # [g/g solution] * [mol/g] = [mol/g solution]
    print(solute_mC[0], solute_mol[0], MolarMass_g)
    write_file_array("solute_conc_mol", solute_mol ) #[mol/g solution]
    write_file_array("solute_conc_g", solute_mC ) #[mol/g solution]

