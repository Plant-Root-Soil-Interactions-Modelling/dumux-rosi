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

MolarMass = 1.2e-2 #[kg/mol] 
MolarMass2 = MolarMass #[kg/mol]
MolarMass_g = MolarMass * 1000 #[g/mol]
MolarMass_g2 = MolarMass2 * 1000 #[g/mol]
exud = 1. # [mol/cm2/d] 

if useC3:
    s.setParameter( "Soil.v_maxL", str(0))
    s.setParameter( "Soil.BC.Bot.C1Type", str(3))
    s.setParameter( "Soil.BC.Top.C1Type", str(3))
    s.setParameter( "Soil.BC.Bot.C1Value", str(exud )) 
    s.setParameter( "Soil.BC.Top.C1Value", str(0 )) 
    s.setParameter("1.Component.MolarMass", str(MolarMass)) 
    s.setParameter("1.Component.LiquidDiffusionCoefficient", "1.e-8") #m^2/s

    s.setParameter( "Soil.BC.Top.C2Type", str(3))
    s.setParameter( "Soil.BC.Bot.C2Type", str(3))
    s.setParameter( "Soil.BC.Bot.C2Value", str(exud/3 )) 
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
    
    
 



ci = 1
# s.setParameter("Soil.IC.C1Z", "0.0002 0.0003 0.0004 0.0005 0.006") 
# s.setParameter("Soil.IC.C1", "0.02 0.03 0.04 0.05 0.6") 
# s.setParameter("Soil.IC.C1", "0.6 0.05 0.04 0.03 0.02") 
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

fig, (ax1, ax2) = plt.subplots(1, 2)

times = [0., 5./24, 10./24.]  # days
s.ddt = 1.e-5

col = ["r*", "b*", "g*", "c*", "m*", "y*", ]
points = s.getDofCoordinates()

x = np.array(s.getSolutionHead())
write_file_array("pressureHead",x.flatten())
write_file_array("coord",points.flatten())
#write_file_array("theta",np.array(s.getWaterContent()).flatten())
#write_file_array("getSaturation",np.array(s.getSaturation()).flatten())
#write_file_array("krs",np.array(s.getKrw()).flatten())
# [g/cm3] * [mol/kg] * [kg/g] = [mol/cm3]


rho_avg = 1#np.array(s.getAvgDensity()).flatten()/1000 #g/cm3
#write_file_array("density", rho_avg) 

solute_mol = s.getSolution_(1).flatten() # [mol/mol solution] == [mol/cm3 solution]
# [g/g solution] * [g solution / cm3 solution] * [mol/g] = [mol/cm3 solution]
# solute_mol = solute_mC * rho_avg / MolarMass_g 
# write_file_array("solute_massFr", solute_mC)  # [g/g solution]
# write_file_array("solute_massKonz", solute_mC * rho_avg) # [g solution / cm3 solution] 
write_file_array("solute_molKonz", solute_mol) # [g solution / cm3 solution] 


solute_mC = s.getSolution_(2).flatten() # [g/g solution]
# [g/g solution] * [g solution / cm3 solution] * [mol/g] = [mol/cm3 solution]
solute_mol = solute_mC * rho_avg / MolarMass_g2 
write_file_array("solute2_massFr", solute_mC)  # [g/g solution]
write_file_array("solute2_massKonz", solute_mC * rho_avg) # [g solution / cm3 solution] 
write_file_array("solute2_molKonz", solute_mol) # [g solution / cm3 solution] 



for i, dt in enumerate(np.diff(times)):

    if rank == 0:
        print("*****", "external time step", dt, " d, simulation time", s.simTime, "d, internal time step", s.ddt, "d")

    s.solve(dt)
    rho_avg = 1#np.array(s.getAvgDensity()).flatten()/1000 #g/cm3
    #write_file_array("density", rho_avg) 
    
    x = np.array(s.getSolutionHead())
    write_file_array("pressureHead",x.flatten())
    write_file_array("coord",points.flatten())
    #write_file_array("theta",np.array(s.getWaterContent()).flatten())
    #write_file_array("getSaturation",np.array(s.getSaturation()).flatten())
    #write_file_array("krs",np.array(s.getKrw()).flatten())
    #write_file_array("solute_conc", np.array(s.getSolution_(1)).flatten()) 
    #write_file_array("solute_conc2", np.array(s.getSolution_(2)).flatten()) 
    write_file_array("density", np.array(s.getAvgDensity()).flatten()) #kg/m3
    
    solute_mC = s.getSolution_(1).flatten() # [g/g solution]
    # [g/g solution] * [g solution / cm3 solution] * [mol/g] = [mol/cm3 solution]
    solute_mol = solute_mC * rho_avg / MolarMass_g 
    write_file_array("solute_massFr", solute_mC)  # [g/g solution]
    write_file_array("solute_massKonz", solute_mC * rho_avg) # [g solution / cm3 solution] 
    write_file_array("solute_molKonz", solute_mol) # [g solution / cm3 solution] 


    solute_mC = s.getSolution_(2).flatten() # [g/g solution]
    # [g/g solution] * [g solution / cm3 solution] * [mol/g] = [mol/cm3 solution]
    solute_mol = solute_mC * rho_avg / MolarMass_g2 
    write_file_array("solute2_massFr", solute_mC)  # [g/g solution]
    write_file_array("solute2_massKonz", solute_mC * rho_avg) # [g solution / cm3 solution] 
    write_file_array("solute2_molKonz", solute_mol) # [g solution / cm3 solution] 
