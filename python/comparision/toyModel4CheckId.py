import sys; sys.path.append("../modules_fpit/"); 
sys.path.append("../../../CPlantBox/");  
sys.path.append("../../build-cmake/cpp/python_binding/")
sys.path.append("../../../CPlantBox/src/python_modules")


from rosi_richards10c import Richards10CSP  # C++ part (Dumux binding)

from richards import RichardsWrapper  # Python part

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()

def write_file_array(name, data, space =","):
    name2 = './results3D/'+ name+ '.txt'
    with open(name2, 'a') as log:
        log.write(space.join([num for num in map(str, data)])  +'\n')
""" 
Cylindrical 1D model, diffusion only (DuMux), Michaelis Menten

everything scripted, no input file needed, also works parallel with mpiexec
"""

results_dir="./results3D/"
if not os.path.exists(results_dir):
    os.makedirs(results_dir)
else:
    test = os.listdir(results_dir)
    for item in test:
        try:
            os.remove(results_dir+item)
        except:
            pass
usemoles = True
s = RichardsWrapper(Richards10CSP(), usemoles)

s.initialize()

# theta_r, theta_s, alpha, n, Ks
loam = [0.045, 0.43, 0.04, 1.6, 50]

min_b = [0., 0., 0.] 
max_b = [400., 200, 100.]
cell_number = [2,2,2]
# s.createGrid([0.02], [0.6], [nCells])  # [cm]
s.createGrid(min_b, max_b, cell_number, False)
cell_number = np.array(cell_number)
s.setParameter( "Soil.Grid.Cells", str(cell_number)[1:(len(str(cell_number))-1)])

s.setHomogeneousIC(-100., equilibrium = True)  # cm pressure head
    
s.setOuterBC("constantFlux", 0)  #  [cm/day]
s.setInnerBC("constantFlux",0* -10)  #  [cm/day]
s.setParameter("Problem.doSoluteFlow", "0")

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
exud = 1.*0# [mol/cm2/d]#1.0* MolarMass *1000# [mol/cm2/d] * [kg/mol] * [g/kg] =  [g/cm2/d]
Ds = 1e-8 # m^2/s
Dl = 1e-9

numComp = 8
numFluidComp = 2


s.setParameter( "Soil.MolarMass", str(solidMolarMass))
s.setParameter( "Soil.solidDensity", str(solidDensity))

s.setParameter( "Soil.BC.Bot.C1Type", str(2))
s.setParameter( "Soil.BC.Top.C1Type", str(2))
s.setParameter( "Soil.BC.Bot.C1Value", str(exud)) 
s.setParameter( "Soil.BC.Top.C1Value", str(0 )) 

s.setParameter("1.Component.LiquidDiffusionCoefficient", str(Ds)) #m^2/s

s.setParameter( "Soil.BC.Bot.C2Type", str(2))
s.setParameter( "Soil.BC.Top.C2Type", str(2))
s.setParameter( "Soil.BC.Bot.C2Value", str(exud)) 
s.setParameter( "Soil.BC.Top.C2Value", str(0 )) 
s.setParameter("2.Component.LiquidDiffusionCoefficient", str(Dl)) #m^2/s

for i in range(numFluidComp + 1, numComp+1):
    print("Soil.BC.Bot.C"+str(i)+"Type")
    s.setParameter( "Soil.BC.Bot.C"+str(i)+"Type", str(2))
    s.setParameter( "Soil.BC.Top.C"+str(i)+"Type", str(2))
    s.setParameter( "Soil.BC.Bot.C"+str(i)+"Value", str(0)) 
    s.setParameter( "Soil.BC.Top.C"+str(i)+"Value", str(0 )) 
        
    
    
s.setParameter("Soil.betaC", str(0.001 ))
s.setParameter("Soil.betaO", str(0.1 ))
s.setParameter("Soil.C_S_W_thresC", str(0.1 /1e6)) #mol/cm3
s.setParameter("Soil.C_S_W_thresO", str(0.1 /1e6 )) #mol/cm3
s.setParameter("Soil.k_decay", str(0.2 ))
s.setParameter("Soil.k_decay2", str(0.6 ))
#s.setParameter("Soil.k_decay3", str(1 ))
s.setParameter("Soil.k_DC", str(1  )) # 1/d
s.setParameter("Soil.k_DO", str(1  )) # 1/d
#s.setParameter("Soil.k_growthBis", str(1 )) #bool
s.setParameter("Soil.k_growthC", str(0.2 ))
s.setParameter("Soil.k_growthO", str(0.2 ))
s.setParameter("Soil.K_L", str(8.3))#[mol/cm3]
s.setParameter("Soil.k_phi", str(0.1 ))
s.setParameter("Soil.k_RC", str(0.1 ))
s.setParameter("Soil.k_RO", str(0.1 ))

k_SC = 1
k_SO = 10
s.setParameter("Soil.k_SC", str(k_SC )) #cm^3/mol/d
#s.setParameter("Soil.k_SCBis", str(k_SC )) #cm^3/mol/d
s.setParameter("Soil.k_SO", str(k_SO )) #cm^3/mol/d
#s.setParameter("Soil.k_SOBis", str(k_SO )) #cm^3/mol/d

m_maxC = 0.1 
m_maxO = 0.02 
s.setParameter("Soil.m_maxC", str(m_maxC  ))# 1/d
s.setParameter("Soil.m_maxO", str(m_maxO  ))# 1/d
s.setParameter("Soil.micro_maxC", str(2 ))# 1/d
s.setParameter("Soil.micro_maxO", str(0.01 ))# 1/d
s.setParameter("Soil.v_maxL", str(1.5 ))#[d-1]

k_sorp = 0.2*100
s.setParameter("Soil.k_sorp", str(k_sorp)) # mol / cm3
f_sorp = 0.9
s.setParameter("Soil.f_sorp", str(f_sorp)) #[-]
CSSmax = 1e-4*10000*0
s.setParameter("Soil.CSSmax", str(CSSmax)) #[mol/cm3 scv]
s.setParameter("Soil.alpha", str(0.1*0)) #[1/d]

C_S = 6.0 *0#mol/cm3 wat
s.setParameter("Soil.IC.C1", str(C_S/ molarDensityWat) )  #mol/cm3 / mol/cm3 = mol/mol 


C_L = 10  *0#mol/cm3 wat
s.setParameter("Soil.IC.C2", str(C_L/ molarDensityWat) )  #mol/cm3 / mol/cm3 = mol/mol
COa = 0.011 * 1e6*0#mol C / m3 space
s.setParameter("Soil.IC.C3",str(COa/ bulkDensity_m3)) #mol C / mol Soil 
#s.setParameter("Soil.IC.C3", str(0.009127163)) #[mol/mol soil] == 233.8 [mol/m3 bulk soil]
COd = 0.05 * 1e6 *0#mol C / m3 space
s.setParameter("Soil.IC.C4", str(COd/bulkDensity_m3 ))
CCa = 0.011 * 1e6*0#mol C / m3 space
s.setParameter("Soil.IC.C5", str(CCa/ bulkDensity_m3)) 
CCd = 0.05 * 1e6*0  #mol C / m3 space
s.setParameter("Soil.IC.C6", str(CCd/bulkDensity_m3 ))
CSS2_init = CSSmax*1e6 * (C_S/(C_S+ k_sorp)) * (1 - f_sorp) #mol C/ m3 scv

s.setParameter("Soil.IC.C7", str(0))#CSS2_init/bulkDensity_m3))#[1:(len(str(CSS2_init/bulkDensity_m3))-1)] ) #mol C / mol scv
s.setParameter("Soil.IC.C8", str(0 ))


s.setVGParameters([loam])



s.setParameter("Problem.reactionExclusive", "1")
s.setParameter("Problem.EnableGravity", "true")
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

times = [0., 0.4]#5./24., 10./24.]  # days
s.ddt = 1.e-5

points = s.getCellCenters()#.flatten()
print('getCellCenters',points)



cidx = np.array(s.base.getCellIndices())
cidx_sorted = np.sort(cidx)
if (cidx != cidx_sorted).any():
    print('too many threads for  the number of cells: ,',cidx,cidx_sorted)
    raise Exception
print('cidx',cidx)
print('cidx_sorted',cidx_sorted)


