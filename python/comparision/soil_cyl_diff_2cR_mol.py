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
s.setOuterBC("constantFluxCyl",  0)  #  [cm/day]
s.setInnerBC("constantFluxCyl", -1)  #  [cm/day]

#s.setICZ_solute(0.)  # [kg/m2] 

MolarMass = 1.2e-2 #[kg/mol] 
MolarMass2 = MolarMass #[kg/mol]
MolarMass_g = MolarMass * 1000 #[g/mol]
MolarMass_g2 = MolarMass2 * 1000 #[g/mol]
exud = 1 # [mol/cm2/d] 



if useC3:
    s.setParameter( "Soil.v_maxL", str(5e5))# [1/d]
    s.setParameter( "Soil.C_oa", str(0.0002338 ))#mol/cm3
    s.setParameter( "Soil.K_L", str(10e-2 ))#[mol/cm^3]
    s.setParameter( "Soil.BC.Bot.C1Type", str(3))
    s.setParameter( "Soil.BC.Top.C1Type", str(6))
    s.setParameter( "Soil.BC.Bot.C1Value", str(exud )) 
    s.setParameter( "Soil.BC.Top.C1Value", str(0 )) 
    s.setParameter("1.Component.MolarMass", str(MolarMass)) 
    s.setParameter("1.Component.LiquidDiffusionCoefficient", "1.e-8") #m^2/s

    s.setParameter( "Soil.BC.Bot.C2Type", str(3))
    s.setParameter( "Soil.BC.Top.C2Type", str(6))
    s.setParameter( "Soil.BC.Bot.C2Value", str(exud/3)) 
    s.setParameter( "Soil.BC.Top.C2Value", str(0 )) 
    s.setParameter("2.Component.MolarMass", str(MolarMass)) 
    s.setParameter("2.Component.LiquidDiffusionCoefficient", "1.e-9") #m^2/s
    
    # s.setParameter("Soil.IC.C1", "1")
    # s.setParameter("Soil.IC.C2", "1")
    # s.setParameter("Soil.IC.C1Z", "0.0002 0.0003 0.0004 0.0005 0.006") 
    # s.setParameter("Soil.IC.C1", "0.02 0.03 0.04 0.05 0.6") 
    # s.setParameter("Soil.IC.C2Z", "0.0002 0.0003 0.0004 0.0005 0.006") 
    # s.setParameter("Soil.IC.C2", "0.02 0.03 0.04 0.05 0.6") 
    # s.setParameter("Soil.IC.C1", "0.6 0.05 0.04 0.03 0.02") 
else:
    s.setParameter( "Soil.BC.Bot.SType", str(3))
    s.setParameter( "Soil.BC.Top.SType", str(6))
    s.setParameter( "Soil.BC.Bot.CValue", str(exud )) 
    s.setParameter( "Soil.BC.Top.CValue", str(0 )) 
    s.setParameter("Component.MolarMass", str(MolarMass)) 
    s.setParameter("Component.LiquidDiffusionCoefficient", "1.e-9") #m^2/s
    
    
 

s.setParameter("Component.BufferPower", "0")  
s.setVGParameters([loam])
s.setParameter("Problem.EnableGravity", "false")
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
s.setCriticalPressure(-15000)  # cm pressure head


if rank == 0:
    print(s)  # TODO make super nice


times = [0.,1.,2.]  # days
s.ddt = 1.e-5

points = s.getDofCoordinates()

x = np.array(s.getSolutionHead())
write_file_array("pressureHead",x.flatten())
write_file_array("coord",points.flatten())

solute_mol = s.getSolution_(1).flatten() # [mol/mol solution]
solute_mol2 = s.getSolution_(2).flatten() # [mol/mol solution]

molarMassWat = 18 # [g/mol]
densityWat = 1 #[g/cm3]
# [mol/cm3] = [g/cm3] /  [g/mol] 
molarDensityWat =  densityWat / molarMassWat # [mol/cm3]
write_file_array("solute_molFr", solute_mol) # [mol solute/mol solution]
write_file_array("solute_molKonz", solute_mol * molarDensityWat) # [mol solute / cm3 solution] 
write_file_array("solute2_molFr", solute_mol2) # [mol solute/mol solution]
write_file_array("solute2_molKonz", solute_mol2 * molarDensityWat) # [mol solute / cm3 solution] 



for i, dt in enumerate(np.diff(times)):

    if rank == 0:
        print("*****", "external time step", dt, " d, simulation time", s.simTime, "d, internal time step", s.ddt, "d")

    s.solve(dt)
    
    x = np.array(s.getSolutionHead())
    write_file_array("pressureHead",x.flatten())
    write_file_array("coord",points.flatten())
    
    solute_mol = s.getSolution_(1).flatten() # [mol/mol solution]
    solute_mol2 = s.getSolution_(2).flatten() # [mol/mol solution]

    molarMassWat = 18 # [g/mol]
    densityWat = 1 #[g/cm3]
    # [mol/cm3] = [g/cm3] /  [g/mol] 
    molarDensityWat =  densityWat / molarMassWat # [mol/cm3]
    write_file_array("solute_molFr", solute_mol) # [mol solute/mol solution]
    write_file_array("solute_molKonz", solute_mol * molarDensityWat) # [mol solute / cm3 solution] 
    write_file_array("solute2_molFr", solute_mol2) # [mol solute/mol solution]
    write_file_array("solute2_molKonz", solute_mol2 * molarDensityWat) # [mol solute / cm3 solution] 
