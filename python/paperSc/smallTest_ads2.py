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

from richards import RichardsWrapper  # Python part
from richards_no_mpi import RichardsNoMPIWrapper  # Python part of cylindrcial 

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import van_genuchten as vg     
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()

def write_file_array(name, data, space =","):
    name2 = './results/smallTest_ads2/'+ name+ '.txt'
    with open(name2, 'a') as log:
        log.write(space.join([num for num in map(str, data)])  +'\n')
""" 
Cylindrical 1D model, diffusion only (DuMux), Michaelis Menten

everything scripted, no input file needed, also works parallel with mpiexec
"""

results_dir="./results/smallTest_ads2/"
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
#s = RichardsWrapper(Richards10CCylFoam(), usemoles)
s=RichardsNoMPIWrapper(Richards10CCylFoam(), usemoles)

s.numComp = 8
s.numFluidComp = 2
s.initialize()

def getTotCContent( kjg, cyl, l, init=False):
    vols = cyl.getCellSurfacesCyl()  * l  #cm3 scv
    totC = 0
    for i in range(cyl.numComp):
        isDissolved = (i < 2)
        totC += cyl.getContentCyl(i+1, isDissolved,l)#mol
        print(i+1, sum(cyl.getContentCyl(i+1, isDissolved,l)))
    # mol/cm3
    C_S_W = np.array(cyl.getSolution_(1)).flatten()*cyl.molarDensityWat
    if init:
        css1 = cyl.CSSmax * (C_S_W/(C_S_W+cyl.k_sorp)) * cyl.f_sorp
    else:
        css1 = np.array(s.base.getCSS1_out()).flatten()/1e6 #mol C / cm3 scv
        assert css1[-1] == 0.
        css1 = css1[:(len(css1)-1)]
        
        
    totC += css1*vols
    #print("CSW", np.array(cyl.getSolution_(1)).flatten()*cyl.molarDensityWat)
    #print("css1Test",cyl.CSSmax * (C_S_W/(C_S_W+cyl.k_sorp)) * cyl.f_sorp*vols,  (C_S_W/(C_S_W+cyl.k_sorp)) )
    print("css1", sum(css1*vols))
    assert len(np.array(totC).shape)==1
    return totC
    

nCells = 10
logbase = 0.5
r_in = 0.02
r_out = 0.2500184914854141
l = 1.0000000000000002 #length in cm
doLogarithmic = True
if doLogarithmic:
    s.points = np.logspace(np.log(r_in) / np.log(logbase), np.log(r_out) / np.log(logbase), nCells, base = logbase)
    #s.points = np.array([0.05,  0.06747797, 0.09106553, 0.12289834, 0.1658586,  0.22383603,
    #                     0.30208002, 0.40767492 ,0.55018151, 0.74250263])
    s.createGrid1d(s.points)
    nCells -= 1
else:
    s.createGrid([r_in], [r_out], [nCells])  # [cm]
s.setParameter( "Soil.Grid.Cells", str(nCells))

s.setParameter("Problem.reactionExclusive", "0")
s.setHomogeneousIC(-97.5)  # cm pressure head
#QflowOut =  0.004092286784778323*0
QCflowOut =np.array([-3.66483704e-06, -1.15332047e-11]  )
#wout = 0.1*0 
#win = -0.008902766607560108
wout = 0.1 
QflowOut =wout * (2 * np.pi * r_out * l) 
win = -0.057400356841653545
s.setOuterBC("constantFluxCyl", 0.)  #  [cm/day]
s.setInnerBC("constantFluxCyl", win)  #  [cm/day]
WatBC =  win* (2 * np.pi * r_in * l) + QflowOut#wout * (2 * np.pi * r_out * l) 
#s.setICZ_solute(0.)  # [kg/m2] 



exud = 0.
exuds_in = exud
exudl_in = exud

Qexud = (exuds_in+ exudl_in)* (2 * np.pi * r_in * l) +sum(QCflowOut)
molarMassWat = 18. # [g/mol]
densityWat = 1. #[g/cm3]
# [mol/cm3] = [g/cm3] /  [g/mol] 
molarDensityWat =  densityWat / molarMassWat # [mol/cm3] 
s.molarDensityWat = molarDensityWat

paramSet = pd.read_csv('./TraiRhizoparam_ordered.csv').loc[1646]
mg_per_molC = 12000
s.betaC = paramSet['beta_C'] # -
s.betaO = paramSet['beta_O'] # -

## ATT! C_S_W_thres has to be > 0
s.C_S_W_thresC =  paramSet['C_thres,C']/mg_per_molC# mgC/cm3 water => mol/cm3 water
s.C_S_W_thresO =  paramSet['C_thres,O']/mg_per_molC# mgC/cm3 water => mol/cm3 water

s.setParameter( "Soil.MolarMass", str(s.solidMolarMass))
s.setParameter( "Soil.solidDensity", str(s.solidDensity))
Ds = paramSet['DS_W'] #cm^2/d
Dl = 0.003456 #cm^2/d

s.k_SC =  paramSet['k_C,S'] * mg_per_molC/paramSet['theta'] #[cm^3 water /mgC/d] => [cm^3 bulk soil/mol/d ]
s.k_DC =  paramSet['k_d,C'] #1/d
s.k_DO =  paramSet['k_d,O']#1/d

s.K_L =  paramSet['K_L'] /mg_per_molC#[mol cm-3 soil]
s.k_SO =  paramSet['k_O,S'] * mg_per_molC/paramSet['theta']

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
s.k_growthO = paramSet['Y_O'] #- 
s.CSSmax = paramSet['CSS_max']/mg_per_molC*0. # mol C/cm3 bulk soil
# can cause issue
s.k_sorp = paramSet['k_sorp'] /mg_per_molC# mol C/cm3 bulk soil
s.alpha = 0.1# -
s.f_sorp = 0.5
s.k_phi = 0.1


C_S = paramSet['CS_init'] /mg_per_molC## in mol/cm3 water
C_L = paramSet['CL_init'] /mg_per_molC## in mol/cm3 water
unitConversion = 1.0e6 # mol/cm3  => mol/m3 
s.bulkDensity =  paramSet['ro_B']*1000 #g/cm3 => kg/m3
s.solidDensity = 2650 # [kg/m^3 solid] #taken from google docs TraiRhizo
s.solidMolarMass = 60.08e-3 # [kg/mol] 
# [mol / m3 solid] =[kg/m^3 solid] / [kg/mol] 
s.solidMolDensity = s.solidDensity/s.solidMolarMass

# theta_r, theta_s, alpha, n, Ks
s.soil = [0.045, np.nan, 0.04, 1.6, 50]

s.soil[1] = 1- s.bulkDensity/s.solidDensity #== size of air volume
s.vg_soil = vg.Parameters(s.soil) 
# [mol / m3 scv] = [mol / m3 solid] * [m3 solid /m3 space]
# BC
# [mol / m3 scv] = [mol / m3 solid] * [m3 solid /m3 space]
s.bulkDensity_m3 = s.solidMolDensity*(1.- s.vg_soil.theta_S)

CSS2_init = s.CSSmax*1e6 * (C_S/(C_S+ s.k_sorp*1e6)) * (1 - s.f_sorp)#mol C/ m3 scv



s.css1Function = 0
s.C_aOLim=1.e-10*0.
s.C_aCLim=1.e-10*0.
s.setParameter("Soil.C_aOLim", str(s.C_aOLim)) #[molC/cm3 scv]
s.setParameter("Soil.C_aCLim", str(s.C_aCLim)) #[molC/cm3 scv]

s.setParameter( "Soil.css1Function", str(s.css1Function))
s.ICcc = np.array([C_S *unitConversion,
                   C_L*unitConversion,
                    9.16666666666667e-07* unitConversion,
                    8.33333333333333e-06* unitConversion,
                    8.33333333333333e-07* unitConversion,
                    8.33333333333333e-06* unitConversion,
                    CSS2_init, 0.])# in mol/m3 water or mol/m3 scv




numComp = 8
numFluidComp = 2
s.numFluidComp = numFluidComp
s.numComp = numComp


for i in range(s.numComp):
    molarC = s.ICcc[i] / s.phaseDensity(isDissolved = (i < s.numFluidComp)) #mol/m3 to mol/mol
    print('C',i+1,s.ICcc[i], s.phaseDensity(isDissolved = (i < s.numFluidComp)),molarC )
    s.setParameter( "Soil.IC.C"+str(i+1), str(molarC ))

s.setParameter( "Soil.BC.Bot.C1Type", str(3))
s.setParameter( "Soil.BC.Top.C1Type", str(3))
s.setParameter( "Soil.BC.Bot.C1Value", str(exuds_in)) 
s.setParameter( "Soil.BC.Top.C1Value", str(0.)) 
                                               
                                                        
                                               

s.setParameter("1.Component.LiquidDiffusionCoefficient", str(Ds)) #m^2/s

s.setParameter( "Soil.BC.Bot.C2Type", str(3))
s.setParameter( "Soil.BC.Top.C2Type", str(3))
s.setParameter( "Soil.BC.Bot.C2Value", str(exudl_in)) 
s.setParameter( "Soil.BC.Top.C2Value", str(0.)) 
s.setParameter("2.Component.LiquidDiffusionCoefficient", str(Dl)) #m^2/s

for i in range(s.numFluidComp + 1, s.numComp+1):
    print("Soil.BC.Bot.C"+str(i)+"Type")
    s.setParameter( "Soil.BC.Bot.C"+str(i)+"Type", str(3))
    s.setParameter( "Soil.BC.Top.C"+str(i)+"Type", str(3))
    s.setParameter( "Soil.BC.Bot.C"+str(i)+"Value", str(0)) 
    s.setParameter( "Soil.BC.Top.C"+str(i)+"Value", str(0 )) 
        
    

    
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


s.setVGParameters([s.soil])



s.setParameter("Problem.verbose", "0")

s.setParameter("Newton.MaxRelativeShift", "1e-8")


s.initializeProblem()

qOut = s.distributeSource(wout * (2 * np.pi * r_out * l) , 0, l, s.numFluidComp)
s.setCriticalPressure(-15000)  # cm pressure head


valueTopBC = s.distributeSources(QCflowOut,
    np.array([nc+1 for nc in range(s.numFluidComp+1)]),
                                 l, s.numFluidComp)


times = np.array([0.,0.001, 1., 2.] )/24. # days
s.ddt = 1.e-5


buTotCBefore = getTotCContent(0, s,l, init=True)
buWBefore_ =  s.getWaterVolumesCyl(l)

doReset = True
for i, dt in enumerate(np.diff(times)):

    assert (sum(buWBefore_) +WatBC*dt  ) > 0
    assert (sum(buTotCBefore) +Qexud*dt ) > 0
    if rank == 0:
        print("*****", "external time step", dt, " d, simulation time", s.simTime, 
              "d, internal time step", s.ddt, "d")

    s.solve(dt, maxDt = 1)
    
    buTotCAfter = getTotCContent(0, s,l)
    buWAfter_ =  s.getWaterVolumesCyl(l)
    
    print("content, before",sum(buTotCBefore),'after', sum(buTotCAfter), 'added',Qexud*dt)
    print(sum(buTotCBefore) +Qexud*dt - sum(buTotCAfter) )
    print("% error sollutes ",(sum(buTotCBefore) +Qexud*dt - sum(buTotCAfter))/(sum(buTotCBefore) +Qexud*dt)*100 )
    
    print("water",sum(buWBefore_), sum(buWAfter_), WatBC*dt)
    print(sum(buWBefore_) +WatBC*dt - sum(buWAfter_) )
    print("% error water ",(sum(buWBefore_) +WatBC*dt - sum(buWAfter_))/sum(buWAfter_)*100 )
    #print("theta", np.array(s.getWaterContent()).flatten() )
    #print("ccs1",  np.array(s.base.getCSS1_out()).flatten())
    
    if doReset:
        s.reset()
    else:
        buTotCBefore = buTotCAfter
        buWBefore_  = buWAfter_



