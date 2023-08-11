import sys; sys.path.append("../modules/"); sys.path.append("../../../CPlantBox/");  sys.path.append("../../build-cmake/cpp/python_binding/")
sys.path.append("../../../CPlantBox/src/python_modules")


useC3 = True
if useC3:
    from rosi_richards5c_cyl import Richards5CCylFoam  # C++ part (Dumux binding)
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
    s = RichardsWrapper(Richards5CCylFoam())
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

molarMassWat = 18. # [g/mol]
densityWat = 1. #[g/cm3]
# [mol/cm3] = [g/cm3] /  [g/mol] 
molarDensityWat =  densityWat / molarMassWat # [mol/cm3]  
m3_2_cm3 = 1e6 # cm3/m3

MolarMass = 1.8e-2 #[kg/mol] 0.02003 #[kg/mol]
exud = 0. # [mol/cm2/d]#1.0* MolarMass *1000# [mol/cm2/d] * [kg/mol] * [g/kg] =  [g/cm2/d]
Dl = 1.e-8

numComp = 6
numFluidComp = 2

if useC3:
    s.setParameter( "Soil.MolarMass", str(60.08e-3))
    s.setParameter( "Soil.solidDensity", str(2700))
    
    s.setParameter( "Soil.BC.Bot.C1Type", str(3))
    s.setParameter( "Soil.BC.Top.C1Type", str(3))
    s.setParameter( "Soil.BC.Bot.C1Value", str(exud)) 
    s.setParameter( "Soil.BC.Top.C1Value", str(0 )) 
    
    s.setParameter("1.Component.LiquidDiffusionCoefficient", str(Dl)) #m^2/s

    s.setParameter( "Soil.BC.Bot.C2Type", str(3))
    s.setParameter( "Soil.BC.Top.C2Type", str(3))
    s.setParameter( "Soil.BC.Bot.C2Value", str(exud/3.)) 
    s.setParameter( "Soil.BC.Top.C2Value", str(0 )) 
    s.setParameter("2.Component.LiquidDiffusionCoefficient", str(Dl/2.)) #m^2/s
    
    for i in range(numFluidComp + 1, numComp+1):
        print("Soil.BC.Bot.C"+str(i)+"Type")
        s.setParameter( "Soil.BC.Bot.C"+str(i)+"Type", str(3))
        s.setParameter( "Soil.BC.Top.C"+str(i)+"Type", str(3))
        s.setParameter( "Soil.BC.Bot.C"+str(i)+"Value", str(0)) 
        s.setParameter( "Soil.BC.Top.C"+str(i)+"Value", str(0 )) 
        
else:
    s.setParameter( "Soil.BC.Bot.SType", str(3))
    s.setParameter( "Soil.BC.Top.SType", str(6))
    s.setParameter( "Soil.BC.Bot.CValue", str(exud )) 
    s.setParameter( "Soil.BC.Top.CValue", str(0 )) 
    s.setParameter("Component.MolarMass", str(MolarMass)) 
    s.setParameter("Component.LiquidDiffusionCoefficient", "1.e-9") #m^2/s
    
    

s.setParameter("Soil.m_maxOBis", str(0.01 *0))
s.setParameter("Soil.k_DO", str(1 *0 )) 
s.setParameter("Soil.k_RO", str(1 *0 ))
s.setParameter("Soil.k_SOBis", str(0)) #m^3/mol/d
s.setParameter("Soil.k_SO", str(3.5e-05 )) #m^3/mol/d
s.setParameter("Soil.k_growthBis", str(0)) #m^3/mol/d


s.setParameter("Soil.m_maxO", str(0.01))
s.setParameter("Soil.micro_maxO", str(0.02))
s.setParameter("Soil.k_decay2", str(0.5 ))
s.setParameter("Soil.k_decay", str(0.75 ))
s.setParameter("Soil.betaO", str(0.1 ))
s.setParameter("Soil.k_growthO", str(0.8 ))
s.setParameter("Soil.k_phi", str(0.1 ))
s.setParameter("Soil.C_S_W_thresO", str(0.001/(12*1000) )) #mol/cm3

C_S = 5e5 # mol/m^3
C_L = 1e5 # mol/m^3

s.setParameter("Soil.IC.C1", str(C_S / m3_2_cm3 / molarDensityWat)) # [mol/m3] / [cm3/m3] / [mol/cm3] = [mol C/mol water]
s.setParameter("Soil.IC.C2", str(C_L / m3_2_cm3 / molarDensityWat)) #[mol C/mol water]

COa = 0.0002338 * 1e6 #mol C / m3 space
s.setParameter("Soil.IC.C3", str(COa/ ((2700./60.08e-3)*(1.-0.43)) )) #mol C / mol Soil 
#s.setParameter("Soil.IC.C3", str(0.009127163)) #[mol/mol soil] == 233.8 [mol/m3 bulk soil]
s.setParameter("Soil.IC.C4", str(COa/((2700./60.08e-3)*(1.-0.43)) *0/ 5))
s.setParameter("Soil.IC.C5", "0.000") 
s.setParameter("Soil.IC.C6", "0.000") 
#s.setParameter("Soil.IC.C3Z", "0.0002 0.0003 0.0004 0.0005 0.006") 
#s.setParameter("Soil.IC.C3", "0.02 0.03 0.04 0.05 0.6") 
#s.setParameter("Soil.IC.C1", "0.6 0.05 0.04 0.03 0.02") 
s.setParameter("Component.BufferPower", "0")
s.setParameter("Soil.v_maxL", str(0))#[d-1]
s.setParameter("Soil.K_L", str(10e-2))#[mol/cm3]
         

s.setVGParameters([loam])
# # s.setParameter("Assembly.NumericDifferenceMethod", "1")
# # s.setParameter("TimeLoop.DtInitial", "1")
# s.setParameter("TimeLoop.MaxTimeStepSize", "36000")
# s.setParameter("Problem.EnableGravity", "false")
# s.setParameter("Problem.verbose", "1")
# s.setParameter("Flux.UpwindWeight", "0.5")
# s.setParameter("Newton.EnableAbsoluteResidualCriterion", "true")
# s.setParameter("Newton.EnableChop", "false")
# s.setParameter("Newton.EnableResidualCriterion", "true")
# s.setParameter("Newton.EnableShiftCriterion", "true")
# s.setParameter("Newton.MaxAbsoluteResidual", "1e-10")
# s.setParameter("Newton.MaxRelativeShift", "1e-10")
# s.setParameter("Newton.MaxSteps", "30")
# s.setParameter("Newton.ResidualReduction", "1e-10")
# s.setParameter("Newton.SatisfyResidualAndShiftCriterion", "true")
# s.setParameter("Newton.TargetSteps", "10")
# s.setParameter("Newton.UseLineSearch", "false")
# s.setParameter("Newton.EnablePartialReassembly", "true")


# s.setParameter("LinearSolver.MaxIterations", "250")
# s.setParameter("LinearSolver.PreconditionerIterations", "1")
# s.setParameter("LinearSolver.PreconditionerRelaxation", "1.0")
# s.setParameter("LinearSolver.ResidualReduction", "1e-13")
# s.setParameter("LinearSolver.Verbosity", "0")


# s.setParameter("Grid.Overlap", "0")  #no effec5


s.setParameter("Problem.EnableGravity", "false")
s.setParameter("Problem.verbose", "0")
s.setParameter("Flux.UpwindWeight", "0.5")
s.setParameter("Newton.EnableAbsoluteResidualCriterion", "true")
s.setParameter("Newton.MaxAbsoluteResidual", "1.e-10")
s.setParameter("Newton.EnableChop", "true")
s.setParameter("Newton.EnableResidualCriterion", "true")
s.setParameter("Newton.EnableShiftCriterion", "true")
s.setParameter("Newton.MaxAbsoluteResidual", "1e-10")

s.setParameter("Newton.MaxRelativeShift", "1e-15")

s.setParameter("Newton.MaxSteps", "30")
s.setParameter("Newton.ResidualReduction", "1e-15")
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

times = [0., 5./24., 10./24.]  # days
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

for i in range(numFluidComp):
    write_file_array("solute_conc"+str(i+1), np.array(s.getSolution_(i+1)).flatten()*molarDensityWat) 
for i in range(numFluidComp, numComp):
    write_file_array("solute_conc"+str(i+1), np.array(s.getSolution_(i+1)).flatten()* ((2700./60.08e-3)*(1.-0.43)) ) 
#write_file_array("solute_conc2", np.array(s.getSolution_(2)).flatten()) 
write_file_array("solute_conc1bis", np.array(s.getSolution_(1)).flatten() * molarDensityWat) 
write_file_array("solute_conc2bis", np.array(s.getSolution_(2)).flatten() * molarDensityWat) 
                                                                         



for i, dt in enumerate(np.diff(times)):

    if rank == 0:
        print("*****", "external time step", dt, " d, simulation time", s.simTime, "d, internal time step", s.ddt, "d")

    s.solve(dt)

    x = np.array(s.getSolutionHead())
    write_file_array("pressureHead",x.flatten())
    #print(x.flatten())
    write_file_array("coord",points.flatten())
    write_file_array("theta",np.array(s.getWaterContent()).flatten())
    write_file_array("getSaturation",np.array(s.getSaturation()).flatten())
    write_file_array("krs",np.array(s.getKrw()).flatten())
    for i in range(numFluidComp):
        write_file_array("solute_conc"+str(i+1), np.array(s.getSolution_(i+1)).flatten()*molarDensityWat) 
    for i in range(numFluidComp, numComp):
        write_file_array("solute_conc"+str(i+1), np.array(s.getSolution_(i+1)).flatten()* ((2700./60.08e-3)*(1.-0.43)) ) 
    #write_file_array("solute_conc1", np.array(s.getSolution_(1)).flatten()) 
    #write_file_array("solute_conc2", np.array(s.getSolution_(2)).flatten()) 
    write_file_array("solute_conc1bis", np.array(s.getSolution_(1)).flatten() * molarDensityWat) 
                                          
                                                                                                  
    write_file_array("solute_conc2bis", np.array(s.getSolution_(2)).flatten() * molarDensityWat) 
                                                                             



