import sys; sys.path.append("../modules_fpit/"); sys.path.append("../../../CPlantBox/");  
sys.path.append("../../build-cmake/cpp/python_binding/")
sys.path.append("../../../CPlantBox/src/python_modules")

sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../modules/");
sys.path.append("/data");
sys.path.append("../../../CPlantBox/src/python_modules");
sys.path.append("../../../CPlantBox/src/functional/");
sys.path.append("../../../CPlantBox/src/rsml/");
sys.path.append("../../../CPlantBox/src/visualisation/")
sys.path.append("../../../CPlantBox/src/structural/")
sys.path.append("../../../CPlantBox/src/external/")
sys.path.append("../../../CPlantBox/src/")
sys.path.append("../../../CPlantBox/");

from rosi_richards10c_cyl import Richards10CCylFoam  # C++ part (Dumux binding)
from rosi_richards10c import Richards10CSP  # C++ part (Dumux binding), macroscopic soil model

# from richards import RichardsWrapper  # Python part
from richards_no_mpi import RichardsNoMPIWrapper as RichardsWrapper  # Python part of cylindrcial model (a single cylindrical model is not allowed to run in parallel)

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import functional.van_genuchten as vg
from functional.xylem_flux import *
from functional.phloem_flux import PhloemFluxPython  # root system Python hybrid solver
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()

def write_file_array(name, data, space =",", directory_ ="./results/", fileType = '.txt', allranks = False ):
    np.array(data).reshape(-1)
    try:
        if (rank == 0) or allranks:
            name2 = directory_+ name+ fileType
            #print('write_file_array',name)
            with open(name2, 'a') as log:
                log.write(space.join([num for num in map(str, data)])  +'\n')
    except:
        print(name, data,data.shape)
        raise Exception
    
""" 
Cylindrical 1D model, diffusion only (DuMux), Michaelis Menten

everything scripted, no input file needed, also works parallel with mpiexec
"""

results_dir="./results/"
dirResults = results_dir
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
s = RichardsWrapper(Richards10CSP(), usemoles)  # water and N solute          

min_b = [-5, -5, -10.] 
max_b = [5, 5, 0.] 
cell_number = [5,5,20]
p_mean_ = -100
do1D = False
s.soil_ = [0.045, 0.43, 0.04, 1.6, 50]
s.vg_soil = vg.Parameters(s.soil_) 
#@see dumux-rosi\cpp\python_binding\solverbase.hh
s.betaC = 0.001 
s.betaO = 0.1 
C_S = 0. # in mol/m3 water
s.C_S_W_thresC = 0.1/1e6#C_S/1e6 # in mol/cm3 water
s.C_S_W_thresO = 0.05/1e6
s.k_decay = 0.2
s.k_decay2 = 0.6 
s.k_DC = 1. 
s.k_DO = 1. 
s.k_growthC = 0.2
s.k_growthO = 0.2 
s.K_L = 8.3
s.k_phi = 0.1 
s.k_RC = 0.1 
s.k_RO = 0.1 

s.k_sorp = 0.4*1e6
s.f_sorp = 0.5

demoType = "dumux_w"
if demoType == "dumux_10c":
    s.CSSmax =C_S/10/1e6 # 1e-4*10000*0.
    s.alpha =0.1# 0.
    unitConversion = 1e3
    doBio = 1.
    CSS2_init = s.CSSmax*1e6 * (C_S/(C_S+ s.k_sorp*1e6)) * (1 - s.f_sorp)#mol C/ m3 scv
elif demoType == "dumux_3c":
    s.CSSmax =0. # 1e-4*10000*0.
    s.alpha =0.1# 0.
    unitConversion = 1e3
    doBio = 0.
    CSS2_init = s.CSSmax*1e6 * (C_S/(C_S+ s.k_sorp*1e6)) * (1 - s.f_sorp)#mol C/ m3 scv
elif demoType == "dumux_w":
    s.CSSmax = 0.
    s.alpha = 0.
    doBio = 0.
    unitConversion = 0.
    CSS2_init = 0.
else:
    print('demoType',demoType)
    raise Exception


s.C_aOLim=1.e-10*0.
s.C_aCLim=1.e-10*0.
s.setParameter("Soil.C_aOLim", str(s.C_aOLim)) #[molC/cm3 scv]
s.setParameter("Soil.C_aCLim", str(s.C_aCLim)) #[molC/cm3 scv]

s.ICcc = np.array([C_S *unitConversion,
                   0 *unitConversion,
                    (C_S/10+s.C_aOLim)* unitConversion *doBio,
                    C_S/10* unitConversion *doBio,
                    (C_S/10+s.C_aCLim)* unitConversion *doBio,
                    C_S/10* unitConversion *doBio,
                    CSS2_init, 0.])# in mol/m3 water or mol/m3 scv


s.initialize()

s.createGrid(min_b, max_b, cell_number, False)  # [cm] 

#cell_number = str(cell_number)
cell_number_ = cell_number
cell_number= s.dumux_str(cell_number)#.replace("[", "");cell_number=cell_number.replace("]", "");cell_number=cell_number.replace(",", "");
s.setParameter( "Soil.Grid.Cells", cell_number)    
s.setParameter("Problem.reactionExclusive", s.dumux_str(int(not do1D)))

# BC
s.setTopBC("noFlux")
s.setBotBC("noFlux") #in acc. with Jorda et al. (2022), however, they assume inflow if h>0
s.solidDensity = 2700 # [kg/m^3 solid]
s.solidMolarMass = 60.08e-3 # [kg/mol] 
# [mol / m3 solid] =[kg/m^3 solid] / [kg/mol] 
s.solidMolDensity = s.solidDensity/s.solidMolarMass
# [mol / m3 scv] = [mol / m3 solid] * [m3 solid /m3 space]
s.bulkDensity_m3 = s.solidMolDensity*(1.- s.vg_soil.theta_S)

s.setParameter( "Soil.MolarMass", str(s.solidMolarMass))
s.setParameter( "Soil.solidDensity", str(s.solidDensity))

s.Ds = 1e-9 # m^2/s
s.Dl = 3e-12
s.numComp = 8
s.numFluidComp = 2
s.setParameter("1.Component.LiquidDiffusionCoefficient", str(s.Ds)) #m^2/s
s.setParameter("2.Component.LiquidDiffusionCoefficient", str(s.Dl)) #m^2/s

s.decay = 0. #1.e-5

for i in range(1, s.numComp+1):
    s.setParameter( "Soil.BC.Bot.C"+str(i)+"Type", str(2))
    s.setParameter( "Soil.BC.Top.C"+str(i)+"Type", str(2))
    s.setParameter( "Soil.BC.Bot.C"+str(i)+"Value", str(0)) 
    s.setParameter( "Soil.BC.Top.C"+str(i)+"Value", str(0 )) 
for i in range(s.numComp):
    molarC = s.ICcc[i] / s.phaseDensity(isDissolved = (i < s.numFluidComp)) #mol/m3 to mol/mol
    s.setParameter( "Soil.IC.C"+str(i+1), str(molarC ))

s.setParameter("Soil.betaC", str(s.betaC))
s.setParameter("Soil.betaO", str(s.betaO ))
s.setParameter("Soil.C_S_W_thresC", str(s.C_S_W_thresC )) #mol/cm3
s.setParameter("Soil.C_S_W_thresO", str(s.C_S_W_thresO )) #mol/cm3
s.setParameter("Soil.k_decay", str(s.k_decay ))
s.setParameter("Soil.k_decay2", str(s.k_decay2))
#s.setParameter("Soil.k_decay3", str(1 ))
s.setParameter("Soil.k_DC", str(s.k_DC )) # 1/d
s.setParameter("Soil.k_DO", str(s.k_DO )) # 1/d
#s.setParameter("Soil.k_growthBis", str(1 )) #bool
s.setParameter("Soil.k_growthC", str(s.k_growthC ))
s.setParameter("Soil.k_growthO", str(s.k_growthO ))
s.setParameter("Soil.K_L", str(s.K_L))#[mol/cm3]
s.setParameter("Soil.k_phi", str(s.k_phi ))
s.setParameter("Soil.k_RC", str(s.k_RC))
s.setParameter("Soil.k_RO", str(s.k_RO))

s.k_SC = 1
s.k_SO = 10
s.setParameter("Soil.k_SC", str(s.k_SC )) #cm^3/mol/d
#s.setParameter("Soil.k_SCBis", str(k_SC )) #cm^3/mol/d
s.setParameter("Soil.k_SO", str(s.k_SO )) #cm^3/mol/d
#s.setParameter("Soil.k_SOBis", str(k_SO )) #cm^3/mol/d

s.m_maxC = 0.1 
s.m_maxO = 0.02 
s.setParameter("Soil.m_maxC", str(s.m_maxC  ))# 1/d
s.setParameter("Soil.m_maxO", str(s.m_maxO  ))# 1/d
s.micro_maxC = 2
s.micro_maxO = 0.01
s.setParameter("Soil.micro_maxC", str(s.micro_maxC ))# 1/d
s.setParameter("Soil.micro_maxO", str(s.micro_maxO ))# 1/d
s.v_maxL = 1.5
s.setParameter("Soil.v_maxL", str(s.v_maxL))#[d-1]

s.setParameter("Soil.k_sorp", str(s.k_sorp)) # mol / cm3
s.setParameter("Soil.f_sorp", str(s.f_sorp)) #[-]
s.setParameter("Soil.CSSmax", str(s.CSSmax)) #[mol/cm3 scv]
s.setParameter("Soil.alpha", str(s.alpha)) #[1/d]


# Paramters
#dumux-rosi\python\modules\richards.py
s.setVGParameters([s.soil_])
#@see dumux-rosi\cpp\python_binding\solverbase.hh
#s.setParameter("Newton.EnableAbsoluteResidualCriterion", "true")
s.MaxRelativeShift = 1e-8
s.setParameter("Newton.MaxRelativeShift", str(s.MaxRelativeShift))
s.setParameter("Problem.verbose", "0")

# IC
if isinstance(p_mean_,(int,float)):
    #print('set pmean float',p_mean_)
    s.setHomogeneousIC(p_mean_, equilibrium = not do1D)  # cm pressure head
elif isinstance(p_mean_,type(np.array([]))):
    pass
else:
    print(type(p_mean_))
    raise Exception

s.initializeProblem()

if isinstance(p_mean_,(int,float)):
    pass
elif isinstance(p_mean_,type(np.array([]))):
    s.setInitialConditionHead(p_mean_)
else:
    print(type(p_mean_))
    raise Exception

s.wilting_point = -15000
s.setCriticalPressure(s.wilting_point)  # for boundary conditions constantFlow, constantFlowCyl, and atmospheric
s.ddt = 1.e-5  # [day] initial Dumux time step

write_file_array('getWaterContent',s.getWaterContent(), directory_ =dirResults, fileType = '.csv')
write_file_array('getSolutionHead',s.getSolutionHead(), directory_ =dirResults, fileType = '.csv')

solute_conc = np.array(s.getSolution(1))
if rank == 0:
    try:
        assert min(solute_conc) >=0
    except:
        print("soil_sol_fluxes", solute_conc)
        print("min(solute_conc)",min(solute_conc))
        raise Exception

cidx = np.array(s.base.getCellIndices())
cidx_sorted = np.sort(cidx)
if (cidx != cidx_sorted).any():
    print('too many threads for  the number of cells: ,',cidx,cidx_sorted)
    raise Exception
print('s.solve')
s.solve( 0.010416666666666666)
print('s.solve_end')
raise Exception

print("finished")

