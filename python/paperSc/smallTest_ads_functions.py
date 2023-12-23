import sys; sys.path.append("../modules_fpit/"); sys.path.append("../../../CPlantBox/");  sys.path.append("../../build-cmake/cpp/python_binding/")
sys.path.append("../../../CPlantBox/src/python_modules")
sys.path.append("../../../CPlantBox/src/functional/");
sys.path.append("../../../CPlantBox/src/rsml/");
sys.path.append("../../../CPlantBox/src/visualisation/")
sys.path.append("../../../CPlantBox/src/structural/")
sys.path.append("../../../CPlantBox/src/external/")
sys.path.append("../../../CPlantBox/");


import matplotlib; matplotlib.use('agg')

from rosi_richards10c_cyl import Richards10CCylFoam  # C++ part (Dumux binding)

from richards_no_mpi import RichardsNoMPIWrapper  # Python part

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import van_genuchten as vg
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()

def write_file_array(name, data, space =","):
    name2 = './results/smallTest_ads/'+ name+ '.txt'
    with open(name2, 'a') as log:
        log.write(space.join([num for num in map(str, data)])  +'\n')
""" 
Cylindrical 1D model, diffusion only (DuMux), Michaelis Menten

everything scripted, no input file needed, also works parallel with mpiexec
"""

def setShape(s,r_in, r_out,doLogarithmic=True):
    nCells = 10
    logbase = 0.5
    s.l = 1.0 #length in cm
    doLogarithmic = True
    s.r_in = r_in
    s.r_out = r_out
    if doLogarithmic:
        s.points = np.logspace(np.log(r_in) / np.log(logbase), np.log(r_out) / np.log(logbase), nCells, base = logbase)
        #s.points = np.array([0.05,  0.06747797, 0.09106553, 0.12289834, 0.1658586,  0.22383603,
        #                     0.30208002, 0.40767492 ,0.55018151, 0.74250263])
        s.createGrid1d(s.points)
        nCells -= 1
    else:
        s.createGrid([r_in], [r_out], [nCells])  # [cm]
    s.setParameter( "Soil.Grid.Cells", str(nCells))
    return s

def setShape1D(s,r_in, r_out,doLogarithmic=True):
    nCells = 10
    logbase = 0.5
    s.l = 1.0000000000000002 #length in cm
    doLogarithmic = True
    s.r_in = r_in
    s.r_out = r_out
    if doLogarithmic:
        s.points = np.logspace(np.log(r_in) / np.log(logbase), np.log(r_out) / np.log(logbase), nCells, base = logbase)
        #s.points = np.array([0.05,  0.06747797, 0.09106553, 0.12289834, 0.1658586,  0.22383603,
        #                     0.30208002, 0.40767492 ,0.55018151, 0.74250263])
        s.createGrid1d(s.points)
        nCells -= 1
    else:
        s.createGrid([r_in], [r_out], [nCells])  # [cm]
    s.setParameter( "Soil.Grid.Cells", str(nCells))
    return s

    
def setDefault(s,css1Function_=0 ):
    s.setParameter("Problem.reactionExclusive", "0")
    molarMassWat = 18. # [g/mol]
    densityWat = 1. #[g/cm3]
    # [mol/cm3] = [g/cm3] /  [g/mol] 
    molarDensityWat =  densityWat / molarMassWat # [mol/cm3] 
    s.molarDensityWat = molarDensityWat

    s.css1Function = css1Function_
    s.setParameter( "Soil.css1Function", str(s.css1Function))
    s.MaxRelativeShift = 1e-8
    s.setParameter("Newton.MaxRelativeShift", str(s.MaxRelativeShift))
    s.setParameter("Problem.verbose", "0")
    s.setParameter("Newton.EnableChop", "true")# force solute mole fraction > 0 and water in possible pressure ranges
    # s.setParameter("Problem.verbose_local_residual", "true")# for debug
    s.setParameter("Flux.UpwindWeight", "1")#very important because we get high solute gradient.
    return s

def setSoilParam(s,paramIdx):
    paramSet = pd.read_csv('./TraiRhizoparam_ordered.csv').loc[paramIdx]
    s.bulkDensity =  paramSet['ro_B']*1000 #g/cm3 => kg/m3
    s.solidDensity = 2650 # [kg/m^3 solid] #taken from google docs TraiRhizo
    s.solidMolarMass = 60.08e-3 # [kg/mol] 

    # theta_r, theta_s, alpha, n, Ks
    s.soil = [0.045, np.nan, 0.04, 1.6, 50]

    s.soil[1] = 1- s.bulkDensity/s.solidDensity #== size of air volume
    s.vg_soil = vg.Parameters(s.soil) 
    # [mol / m3 solid] =[kg/m^3 solid] / [kg/mol] 
    s.solidMolDensity = s.solidDensity/s.solidMolarMass
    # [mol / m3 scv] = [mol / m3 solid] * [m3 solid /m3 space]
    s.bulkDensity_m3 = s.solidMolDensity*(1.- s.vg_soil.theta_S)

    s.setParameter( "Soil.MolarMass", str(s.solidMolarMass))
    s.setParameter( "Soil.solidDensity", str(s.solidDensity))
    s.setVGParameters([s.soil])
    return s


def getBiochemParam(s,paramIdx, noAds):
    s.numFluidComp = 2
    paramSet = pd.read_csv('./TraiRhizoparam_ordered.csv').iloc[paramIdx].to_dict()
    s.mg_per_molC = 12000
    s.betaC = paramSet['beta_C'] # -
    s.betaO = paramSet['beta_O'] # -

    ## ATT! C_S_W_thres has to be > 0
    s.C_S_W_thresC =  paramSet['C_thres,C']/s.mg_per_molC# mgC/cm3 water => mol/cm3 water
    s.C_S_W_thresO =  paramSet['C_thres,O']/s.mg_per_molC# mgC/cm3 water => mol/cm3 water

    Ds = paramSet['DS_W'] #cm^2/d
    Dl = 0.003456 #cm^2/d

    s.k_SC =  paramSet['k_C,S'] * s.mg_per_molC/paramSet['theta'] #[cm^3 water /mgC/d] => [cm^3 bulk soil/mol/d ]
    s.k_DC =  paramSet['k_d,C'] #1/d
    s.k_DO =  paramSet['k_d,O']#1/d

    s.K_L =  paramSet['K_L'] /s.mg_per_molC#[mol cm-3 soil]
    s.k_SO =  paramSet['k_O,S'] * s.mg_per_molC/paramSet['theta']

    s.k_RC =  paramSet['k_r,C']#1/d
    s.k_RO =  paramSet['k_r,O']#1/d

    s.m_maxC =  paramSet['m_max,C']#1/d
    s.m_maxO =  paramSet['m_max,O']#1/d
    s.k_decay2 =  paramSet['p_L'] # -

    s.micro_maxC =  paramSet['u_max,C']#1/d
    s.micro_maxO =  paramSet['u_max,O']#1/d

    s.v_maxL = paramSet['v_max,L'] #1/d
    s.k_decay = paramSet['Y'] # -
    s.k_growthC = paramSet['Y_C'] # -    
    s.k_growthO = paramSet['Y_O'] # -    
    if noAds:
        s.CSSmax = 0.
    else:
        s.CSSmax = paramSet['CSS_max']/s.mg_per_molC # mol C/cm3 bulk soil
    # can cause issue
    s.k_sorp = paramSet['k_sorp'] /s.mg_per_molC# mol C/cm3 bulk soil
    s.alpha = 0.1# -
    s.f_sorp = 0.5
    s.k_phi = 0.1
    

    s.C_aOLim=1.e-10*0.
    s.C_aCLim=1.e-10*0.
    s.Ds = Ds /(24*3600) /10000 # m^2/s
    s.Dl = Dl /(24*3600) /10000# m^2/s
    return s

def setBiochemParam(s):

    s.setParameter("Soil.C_aOLim", str(s.C_aOLim)) #[molC/cm3 scv]
    s.setParameter("Soil.C_aCLim", str(s.C_aCLim)) #[molC/cm3 scv]

    s.setParameter("1.Component.LiquidDiffusionCoefficient", str(s.Ds)) #m^2/s
    s.setParameter("2.Component.LiquidDiffusionCoefficient", str(s.Dl)) #m^2/s


    s.setParameter("Soil.betaC", str(s.betaC)) # -
    s.setParameter("Soil.betaO", str(s.betaO )) # -
    s.setParameter("Soil.C_S_W_thresC", str(s.C_S_W_thresC )) #mol/cm3
    s.setParameter("Soil.C_S_W_thresO", str(s.C_S_W_thresO )) #mol/cm3
    s.setParameter("Soil.k_decay", str(s.k_decay ))
    s.setParameter("Soil.k_decay2", str(s.k_decay2))
    s.setParameter("Soil.k_DC", str(s.k_DC )) # 1/d
    s.setParameter("Soil.k_DO", str(s.k_DO )) # 1/d
    s.setParameter("Soil.k_growthC", str(s.k_growthC ))
    s.setParameter("Soil.k_growthO", str(s.k_growthO ))
    s.setParameter("Soil.K_L", str(s.K_L))#[mol/cm3]
    s.setParameter("Soil.k_phi", str(s.k_phi ))
    s.setParameter("Soil.k_RC", str(s.k_RC))
    s.setParameter("Soil.k_RO", str(s.k_RO))

    s.setParameter("Soil.k_SC", str(s.k_SC )) #cm^3/mol/d
    s.setParameter("Soil.k_SO", str(s.k_SO )) #cm^3/mol/d

    s.setParameter("Soil.m_maxC", str(s.m_maxC  ))# 1/d
    s.setParameter("Soil.m_maxO", str(s.m_maxO  ))# 1/d
    s.setParameter("Soil.micro_maxC", str(s.micro_maxC ))# 1/d
    s.setParameter("Soil.micro_maxO", str(s.micro_maxO ))# 1/d
    s.setParameter("Soil.v_maxL", str(s.v_maxL))#[d-1]

    s.setParameter("Soil.k_sorp", str(s.k_sorp)) # mol / cm3
    s.setParameter("Soil.f_sorp", str(s.f_sorp)) #[-]
    s.setParameter("Soil.CSSmax", str(s.CSSmax)) #[mol/cm3 scv]
    s.setParameter("Soil.alpha", str(s.alpha)) #[1/d]
    return s

def setIC(s,icw, paramIdx):
    s.setHomogeneousIC(icw) #-97.5)  # cm pressure head

    paramSet = pd.read_csv('./TraiRhizoparam_ordered.csv').loc[paramIdx]
    C_S = paramSet['CS_init'] /s.mg_per_molC## in mol/cm3 water
    C_L = paramSet['CL_init'] /s.mg_per_molC## in mol/cm3 water
    
    CSS2_init = s.CSSmax * (C_S/(C_S+ s.k_sorp)) * (1 - s.f_sorp)#mol C/ cm3 scv
    unitConversion = 1.0e6 # mol/cm3  => mol/m3 
    s.ICcc = np.array([C_S *unitConversion,
                       C_L*unitConversion,
                        9.16666666666667e-07* unitConversion,
                        8.33333333333333e-06* unitConversion,
                        8.33333333333333e-07* unitConversion,
                        8.33333333333333e-06* unitConversion,
                        CSS2_init*unitConversion,
                       0.])# in mol/m3 water or mol/m3 scv
    
    for i in range(s.numComp):
        #mol/m3 to mol/mol
        molarC = s.ICcc[i] / s.phaseDensity(isDissolved = (i < s.numFluidComp)) 
        print('C',i+1,s.ICcc[i], s.phaseDensity(isDissolved = (i < s.numFluidComp)),molarC )
        s.setParameter( "Soil.IC.C"+str(i+1), str(molarC ))
    return s

def getBC(s):    
    s.QCflowOut =np.array([-3.66483704e-06, -1.15332047e-11]  )
    wout = 0.1 
    s.QflowOut =wout * (2 * np.pi * s.r_out * s.l) 
    s.win = -0.057400356841653545
    s.WatBC =  s.win* (2 * np.pi * s.r_in * s.l) + s.QflowOut#wout * (2 * np.pi * s.r_out * l) 
    exud = 0.
    s.exuds_in = exud
    s.exudl_in = exud

    s.Qexud = (s.exuds_in+ s.exudl_in)* (2 * np.pi * s.r_in * s.l) +sum(s.QCflowOut)
    return s

def setBC(s):
    
    s.setOuterBC("constantFluxCyl", 0.)  #  [cm/day]
    s.setInnerBC("constantFluxCyl", s.win)  #  [cm/day]
    
    s.setParameter( "Soil.BC.Bot.C1Type", str(3))
    s.setParameter( "Soil.BC.Top.C1Type", str(3))
    s.setParameter( "Soil.BC.Bot.C1Value", str(s.exuds_in)) 
    s.setParameter( "Soil.BC.Top.C1Value", str(0.)) 


    s.setParameter( "Soil.BC.Bot.C2Type", str(3))
    s.setParameter( "Soil.BC.Top.C2Type", str(3))
    s.setParameter( "Soil.BC.Bot.C2Value", str(s.exudl_in)) 
    s.setParameter( "Soil.BC.Top.C2Value", str(0.)) 

    for i in range(s.numFluidComp + 1, s.numComp+1):
        s.setParameter( "Soil.BC.Bot.C"+str(i)+"Type", str(3))
        s.setParameter( "Soil.BC.Top.C"+str(i)+"Type", str(3))
        s.setParameter( "Soil.BC.Bot.C"+str(i)+"Value", str(0)) 
        s.setParameter( "Soil.BC.Top.C"+str(i)+"Value", str(0 )) 

def doInit(s):
    s.initializeProblem()
    s.wilting_point = -15000
    s.setCriticalPressure(s.wilting_point)  # for boundary conditions constantFlow, 
    return s

def setSink(s):
    qOut = s.distributeSource(s.QflowOut, 0, s.l, s.numFluidComp)
    valueTopBC = s.distributeSources(s.QCflowOut,
            np.array([nc+1 for nc in range(s.numFluidComp+1)]),
                                 s.l, s.numFluidComp)
    return s



def forget():
    aaa = np.array(["Soil.MolarMass", s.solidMolarMass, "Soil.solidDensity", s.solidDensity,
    "Soil.betaC", str(s.betaC ),"Soil.betaO", str(s.betaO),"Soil.C_S_W_thresC", str(s.C_S_W_thresC ),"Soil.C_S_W_thresO", str(s.C_S_W_thresO ),
    "Soil.k_decay", str(s.k_decay),"Soil.k_decay2", str(s.k_decay2 ),"Soil.k_DC", str(s.k_DC  ),
    "Soil.k_DO", str(s.k_DO  ),"Soil.k_growthC", str(s.k_growthC),"Soil.k_growthO", str(s.k_growthO),"Soil.K_L", str(s.K_L),"Soil.k_phi", str(s.k_phi ),"Soil.k_RC", str(s.k_RC),"Soil.k_RO", str(s.k_RO ),"Soil.k_SC", str(s.k_SC ),"Soil.k_SO", str(s.k_SO ),"Soil.m_maxC", str(s.m_maxC  ),"Soil.m_maxO", str(s.m_maxO  ),"Soil.micro_maxC", str(s.micro_maxC ),"Soil.micro_maxO", str(s.micro_maxO ),"Soil.v_maxL", str(s.v_maxL),"Soil.k_sorp", str(s.k_sorp),"Soil.f_sorp", str(s.f_sorp),"Soil.CSSmax", str(s.CSSmax),"Soil.alpha", str(s.alpha),"Soil.C_aOLim", str(s.C_aOLim),"Soil.C_aCLim", str(s.C_aCLim),"1.Component.LiquidDiffusionCoefficient", str(s.Ds),"2.Component.LiquidDiffusionCoefficient", str(s.Dl),"Newton.MaxRelativeShift",str(s.MaxRelativeShift),'wiltingPoint',s.wilting_point, 'bulkDensity',s.bulkDensity,'bulkDensity_m3',s.bulkDensity_m3,'solidDensity',s.solidDensity, 'solidMolarMass',s.solidMolarMass,' solidMolDensity', s.solidMolDensity, "Soil.css1Function", str(s.css1Function),'s.soil',s.soil,'s.points',s.points],dtype=object) 

    name2 = './results/smallTest_ads/'+ "params"+ '.csv'
    space ="\n"
    with open(name2, 'a') as log:
        log.write(space.join([num for num in map(str, aaa)])  +'\n')