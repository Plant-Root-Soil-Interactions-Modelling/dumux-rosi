import sys; sys.path.append("../modules_fpit/"); sys.path.append("../../../CPlantBox/");  
sys.path.append("../../build-cmake/cpp/python_binding/")
sys.path.append("../../../CPlantBox/src/python_modules")



# from richards import RichardsWrapper  # Python part
from richards import RichardsWrapper  # Python part of cylindrcial model (a single cylindrical model is not allowed to run in parallel)

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()
from rosi_richards10c import Richards10CSP  # C++ part (Dumux binding), macroscopic soil model

def write_file_array(name, data, space =",", rank_ = rank):
    name2 = './results/'+ str(rank_)+"_"+name+ '.txt'
    if rank_ == 0:
        with open(name2, 'a') as log:
            log.write(space.join([num for num in map(str, data)])  +'\n')
    else:
        with open(name2, 'a') as log:
            log.write(str(data))
    
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
usemoles = True
s = RichardsWrapper(Richards10CSP(), usemoles)

s.initialize()

    
# theta_r, theta_s, alpha, n, Ks
loam = [0.045, 0.43, 0.04, 1.6, 50]


min_b = [-5, -5, -10.] 
max_b = [5, 5, 0.] 
cell_number = [5,5,20]

maxDt_temp = 250/(24*3600)


s.createGrid(min_b, max_b, cell_number, False)  # [cm]
cell_number_ = cell_number
cell_number= s.dumux_str(cell_number)
s.setParameter( "Soil.Grid.Cells", cell_number)

initPHead = -100
s.setHomogeneousIC(initPHead)  # cm pressure head

s.setOuterBC("constantFlux", 0)  #  [cm/day]
s.setInnerBC("constantFlux",0* -10)  #  [cm/day]

s.setParameter("Problem.reactionExclusive", "1")   
#s.setICZ_solute(0.)  # [kg/m2] 

molarMassWat = 18. # [g/mol]
densityWat = 1. #[g/cm3]
# [mol/cm3] = [g/cm3] /  [g/mol] 
molarDensityWat =  densityWat / molarMassWat # [mol/cm3] 
s.molarDensityWat = molarDensityWat
solidDensity = 2700 # [kg/m^3 solid]
solidMolarMass = 60.08e-3 # [kg/mol] 
# [mol / m3 solid] =[kg/m^3 solid] / [kg/mol] 
solidMolDensity = solidDensity/solidMolarMass
# [mol / m3 scv] = [mol / m3 solid] * [m3 solid /m3 space]
bulkDensity_m3 = solidMolDensity*(1.-0.43)

MolarMass = 1.8e-2 #[kg/mol] 0.02003 #[kg/mol]
exud = 1.*0 # [mol/cm2/d]#1.0* MolarMass *1000# [mol/cm2/d] * [kg/mol] * [g/kg] =  [g/cm2/d]
exuds_in = exud*0
exudl_in = exud * 0
Qexud =0.# (exuds_in+ exudl_in)* (2 * np.pi * r_in * l)

Ds = 1e-8 # m^2/s
Dl = 1e-9

numComp = 8
numFluidComp = 2
s.numFluidComp = numFluidComp
s.numComp = numComp

s.solidDensity = solidDensity
s.solidMolarMass=solidMolarMass
s.solidMolDensity=solidMolDensity
s.bulkDensity_m3 =bulkDensity_m3 #mol / m3 bulk soil

s.setParameter( "Soil.MolarMass", str(solidMolarMass))
s.setParameter( "Soil.solidDensity", str(solidDensity))

s.setParameter( "Soil.BC.Bot.C1Type", str(2))
s.setParameter( "Soil.BC.Top.C1Type", str(2))
s.setParameter( "Soil.BC.Bot.C1Value", str(0)) 
s.setParameter( "Soil.BC.Top.C1Value", str(0)) 

s.setParameter("1.Component.LiquidDiffusionCoefficient", str(Ds)) #m^2/s

s.setParameter( "Soil.BC.Bot.C2Type", str(2))
s.setParameter( "Soil.BC.Top.C2Type", str(2))
s.setParameter( "Soil.BC.Bot.C2Value", str(0)) 
s.setParameter( "Soil.BC.Top.C2Value", str(0)) 
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

s.k_sorp = 0.2*10000
s.setParameter("Soil.k_sorp", str(s.k_sorp)) # mol / cm3
s.f_sorp = 0.1
s.setParameter("Soil.f_sorp", str(s.f_sorp)) #[-]
s.CSSmax = 0.
s.setParameter("Soil.CSSmax", str(s.CSSmax)) #[mol/cm3 scv]
s.setParameter("Soil.alpha", str(0.1*0)) #[1/d]
C_S = 0.
s.setParameter("Soil.IC.C1", str(C_S/ molarDensityWat) )  #mol/cm3 / mol/cm3 = mol/mol 

C_L = 10*0  #mol/cm3 wat
s.setParameter("Soil.IC.C2", str(C_L/ molarDensityWat) )  #mol/cm3 / mol/cm3 = mol/mol
COa = 0.011 * 1e6*0#mol C / m3 space
s.setParameter("Soil.IC.C3",str(COa/ bulkDensity_m3)) #mol C / mol Soil 
#s.setParameter("Soil.IC.C3", str(0.009127163)) #[mol/mol soil] == 233.8 [mol/m3 bulk soil]
COd = 0.05 * 1e6*0#mol C / m3 space
s.setParameter("Soil.IC.C4", str(COd/bulkDensity_m3 ))
CCa = 0.011 * 1e6*0#mol C / m3 space
s.setParameter("Soil.IC.C5", str(CCa/ bulkDensity_m3)) 
CCd = 0.05 * 1e6*0  #mol C / m3 space
s.setParameter("Soil.IC.C6", str(CCd/bulkDensity_m3 ))
CSS2_init = s.CSSmax*1e6 * (C_S/(C_S+ s.k_sorp*1e6)) * (1 - s.f_sorp) *0#mol C/ m3 scv

s.setParameter("Soil.IC.C7", str(0))#CSS2_init/bulkDensity_m3))#[1:(len(str(CSS2_init/bulkDensity_m3))-1)] ) #mol C / mol scv
s.setParameter("Soil.IC.C8", str(0 ))


s.setVGParameters([loam])

s.setParameter("Problem.verbose", "0")

s.initializeProblem()

s.setCriticalPressure(-15000)  # cm pressure head


times = [0.,0.010416666666666666]  # days
s.ddt = 1.e-5

points = s.getDofCoordinates()

x = np.array(s.getSolutionHead())

QWexudIn  = 0. # cm3/day
QWexudOut = 0.0069749956986250895 # cm3/day
QWexud    = 0.#QWexudIn + QWexudOut


vols = 1.  #cm3 scv
theta = np.array(s.getWaterContent())
buWBefore = theta*vols
  
for i, dt in enumerate(np.diff(times)):

    if rank == 0:
        print("*****", "external time step", dt, " d, simulation time", s.simTime, "d, internal time step", s.ddt, "d")

    s.solve(dt, maxDt = maxDt_temp)
    
    x = np.array(s.getSolutionHead())
    theta = np.array(s.getWaterContent())
    buWAfter = theta*vols
    
    print('tot wat vol', buWAfter, qOut)   
    print("content",sum(buWBefore), sum(buWAfter), QWexud*dt)
    print(sum(buWBefore) + QWexud*dt - sum(buWAfter) , sum(buWBefore) - sum(buWAfter) )
    buWBefore = buWAfter
    s.reset()

print("finished")

