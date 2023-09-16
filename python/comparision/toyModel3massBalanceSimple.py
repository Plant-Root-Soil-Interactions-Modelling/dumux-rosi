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
usemoles = True
s = RichardsWrapper(Richards10CCylFoam(), usemoles)

s.initialize()

def getTotCContent( kjg, cyl, l, init=False):
    vols = cyl.getCellSurfacesCyl()  * l  #cm3 scv
    totC = 0
    for i in range(cyl.numComp):
        isDissolved = (i < 2)
        totC += cyl.getContentCyl(i+1, isDissolved,l)#mol
        print(i, sum(cyl.getContentCyl(i+1, isDissolved,l)))
    # mol/cm3
    C_S_W = np.array(cyl.getSolution_(1)).flatten()*cyl.molarDensityWat
    if init:
        css1 = cyl.CSSmax * (C_S_W/(C_S_W+cyl.k_sorp)) * cyl.f_sorp
    else:
        css1 = np.array(s.base.getCSS1_out()).flatten()/1e6 #mol C / cm3 scv
        assert css1[-1] == 0.
        css1 = css1[:(len(css1)-1)]
        
        
    totC += css1*vols
    print("css1", sum(css1*vols))
    assert len(np.array(totC).shape)==1
    return totC
    
# theta_r, theta_s, alpha, n, Ks
loam = [0.045, 0.43, 0.04, 1.6, 50]

nCells = 50
r_in = 0.02
r_out = 0.6
l = 1 #length in cm
s.createGrid([r_in], [r_out], [nCells])  # [cm]
s.setParameter( "Soil.Grid.Cells", str(nCells))

s.setHomogeneousIC(-100.)  # cm pressure head
s.setOuterBC("constantFluxCyl", 0.1*0)  #  [cm/day]
s.setInnerBC("constantFluxCyl", -1)  #  [cm/day]

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
exud = 1. # [mol/cm2/d]#1.0* MolarMass *1000# [mol/cm2/d] * [kg/mol] * [g/kg] =  [g/cm2/d]
exuds_in = exud
exudl_in = exud 
exuds_out = - exuds_in / 24
exudl_out = - exudl_in / 3 
Qexud = (exuds_in+ exudl_in)* (2 * np.pi * r_in * l) + (exuds_out+ exudl_out)* (2 * np.pi * r_out * l) 

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

s.setParameter( "Soil.BC.Bot.C1Type", str(3))
s.setParameter( "Soil.BC.Top.C1Type", str(3))
s.setParameter( "Soil.BC.Bot.C1Value", str(exuds_in)) 
s.setParameter( "Soil.BC.Top.C1Value", str(exuds_out)) 

s.setParameter("1.Component.LiquidDiffusionCoefficient", str(Ds)) #m^2/s

s.setParameter( "Soil.BC.Bot.C2Type", str(3))
s.setParameter( "Soil.BC.Top.C2Type", str(3))
s.setParameter( "Soil.BC.Bot.C2Value", str(exudl_in)) 
s.setParameter( "Soil.BC.Top.C2Value", str(exudl_out)) 
s.setParameter("2.Component.LiquidDiffusionCoefficient", str(Dl)) #m^2/s

for i in range(numFluidComp + 1, numComp+1):
    print("Soil.BC.Bot.C"+str(i)+"Type")
    s.setParameter( "Soil.BC.Bot.C"+str(i)+"Type", str(3))
    s.setParameter( "Soil.BC.Top.C"+str(i)+"Type", str(3))
    s.setParameter( "Soil.BC.Bot.C"+str(i)+"Value", str(0)) 
    s.setParameter( "Soil.BC.Top.C"+str(i)+"Value", str(0 )) 
        
    
    
s.setParameter("Soil.betaC", str(0.001 ))
s.setParameter("Soil.betaO", str(0.1 ))
s.setParameter("Soil.C_S_W_thresC", str(0.1 )) #mol/cm3
s.setParameter("Soil.C_S_W_thresO", str(0.05 )) #mol/cm3
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

s.k_sorp = 0.3
s.setParameter("Soil.k_sorp", str(s.k_sorp)) # mol / cm3
s.f_sorp = 0.3
s.setParameter("Soil.f_sorp", str(s.f_sorp)) #[-]
s.CSSmax = 3*0#1e-4*1e5
s.setParameter("Soil.CSSmax", str(s.CSSmax)) #[mol/cm3 scv]
s.setParameter("Soil.alpha", str(0.1*0)) #[1/d]

gradient = False
if gradient:
    #C_S = np.array([0.1, 0.3, 0.4, 0.5, 9])  #mol/cm3 wat
    #s.setParameter("Soil.IC.C1Z", "0.0001 0.0003 0.0004 0.0005 0.009" )  #mol/cm3 / mol/cm3 = mol/mol 
    #s.setParameter("Soil.IC.C1", str(C_S/molarDensityWat)[1:(len(str(C_S))-1)])   #mol/cm3 / mol/cm3 = mol/mol
    minV = r_in
    maxV = r_out
    rez = 2
    C_S = np.array([ minV + (maxV - minV)* (kk/rez) for kk in range(rez + 1)])  # np.array([0.02,6])  #mol/cm3 wat
    print(C_S)
    assert min(C_S) == minV
    assert max(C_S)  == maxV
    ZC_S = C_S/100
    s.setParameter("Soil.IC.C1Z",s.dumux_str(ZC_S))   # "0.0002 0.006" )  #mol/cm3 / mol/cm3 = mol/mol 
    s.setParameter("Soil.IC.C1", s.dumux_str(C_S/molarDensityWat))   #mol/cm3 / mol/cm3 = mol/mol 
else:
    C_S = 1  #mol/cm3 wat
    s.setParameter("Soil.IC.C1", str(C_S/ molarDensityWat) )  #mol/cm3 / mol/cm3 = mol/mol 

C_L = 10  #mol/cm3 wat
s.setParameter("Soil.IC.C2", str(C_L/ molarDensityWat) )  #mol/cm3 / mol/cm3 = mol/mol
COa = 0.011 * 1e6#mol C / m3 space
s.setParameter("Soil.IC.C3",str(COa/ bulkDensity_m3)) #mol C / mol Soil 
#s.setParameter("Soil.IC.C3", str(0.009127163)) #[mol/mol soil] == 233.8 [mol/m3 bulk soil]
COd = 0.05 * 1e6#mol C / m3 space
s.setParameter("Soil.IC.C4", str(COd/bulkDensity_m3 ))
CCa = 0.011 * 1e6#mol C / m3 space
s.setParameter("Soil.IC.C5", str(CCa/ bulkDensity_m3)) 
CCd = 0.05 * 1e6  #mol C / m3 space
s.setParameter("Soil.IC.C6", str(CCd/bulkDensity_m3 ))
CSS2_init = s.CSSmax*1e6 * (C_S/(C_S+ s.k_sorp*1e6)) * (1 - s.f_sorp) *0#mol C/ m3 scv

s.setParameter("Soil.IC.C7", str(0))#CSS2_init/bulkDensity_m3))#[1:(len(str(CSS2_init/bulkDensity_m3))-1)] ) #mol C / mol scv
s.setParameter("Soil.IC.C8", str(0 ))


s.setVGParameters([loam])



s.setParameter("Problem.EnableGravity", "false")
s.setParameter("Problem.verbose", "10")
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

times = [0., 1., 2.]  # days
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
print("len output", np.array(s.getSolution_(1)).flatten().shape, points.shape)
for i in range(numFluidComp):
    write_file_array("solute_conc"+str(i+1), np.array(s.getSolution_(i+1)).flatten()*molarDensityWat ) 
for i in range(numFluidComp, numComp):
    write_file_array("solute_conc"+str(i+1), np.array(s.getSolution_(i+1)).flatten()* bulkDensity_m3 /1e6 ) 

C_S_W = np.array(s.getSolution_(1)).flatten()*molarDensityWat*1e6
cstot_cls = (np.array(s.getWaterContent()).flatten() 
            * C_S_W) 
cstot_css1 = s.CSSmax*1e6 * (C_S_W/(C_S_W+s.k_sorp*1e6)) * s.f_sorp  

cstot_css2 = np.array(s.getSolution_(7)).flatten()*bulkDensity_m3 
cstot = cstot_cls + cstot_css1 + cstot_css2 
write_file_array("css1",cstot_css1/1e6)
write_file_array("css2",cstot_css2/1e6)
write_file_array("cl",cstot_cls/1e6)
write_file_array("cstot",cstot/1e6)


write_file_array("totcs",np.array([np.round(cstot_cls.sum()),
                np.round(cstot_css1.sum()), 
                np.round(cstot_css2.sum()), 
                np.round(cstot.sum())]), "\t")

contents0 = np.array([(cstot_cls * points).sum(),
                (cstot_css1 * points).sum(), 
                (cstot_css2 * points).sum(), 
                (cstot * points).sum()])
write_file_array("contents", contents0/1e6)
#print(contents0)

RF = 1 + s.f_sorp*(1/theta)*s.CSSmax*1e6*((s.k_sorp*1e6)/((s.k_sorp*1e6+C_S_W)**2))
write_file_array("RF",RF)
buTotCBefore = getTotCContent(0, s,l, init=True)

vols = s.getCellSurfacesCyl()  * l  #cm3 scv
print(sum(buTotCBefore))
#raise Exception
for i, dt in enumerate(np.diff(times)):

    if rank == 0:
        print("*****", "external time step", dt, " d, simulation time", s.simTime, "d, internal time step", s.ddt, "d")

    s.solve(dt, maxDt = 2500/(24*3600))
    
    buTotCAfter = getTotCContent(0, s,l)
    
    print("content",sum(buTotCBefore), sum(buTotCAfter), Qexud*dt)
    print(sum(buTotCBefore) +Qexud*dt - sum(buTotCAfter) )
    print("% error",(sum(buTotCBefore) +Qexud*dt - sum(buTotCAfter))/sum(buTotCAfter)*100 )
    buTotCBefore = buTotCAfter
    x = np.array(s.getSolutionHead())
    write_file_array("pressureHead",x.flatten())
    #print(x.flatten())
    write_file_array("coord",points)
    write_file_array("theta",np.array(s.getWaterContent()).flatten())
    write_file_array("getSaturation",np.array(s.getSaturation()).flatten())
    write_file_array("krs",np.array(s.getKrw()).flatten())
    for i in range(numFluidComp):
        write_file_array("solute_conc"+str(i+1), np.array(s.getSolution_(i+1)).flatten()*molarDensityWat) 
    for i in range(numFluidComp, numComp):
        write_file_array("solute_conc"+str(i+1), np.array(s.getSolution_(i+1)).flatten()* bulkDensity_m3 /1e6 ) 
        
        
    C_S_W = np.array(s.getSolution_(1)).flatten()*molarDensityWat*1e6
    cstot_cls = (np.array(s.getWaterContent()).flatten() 
                * C_S_W) #* points
    cstot_css1 = np.array(s.base.getCSS1_out()).flatten()
    cstot_css1 = cstot_css1[:(len(cstot_css1)-1)]
    #print("cstot_css1",cstot_css1,cstot_css1* (vols/1e6), sum(cstot_css1* (vols/1e6)))
    
    RF = np.array(s.base.getRF_out()).flatten()
    RF = RF[:(len(RF)-1)]
    
    # cstot_css1Bis =  CSSmax*1e6 * (C_S_W/(C_S_W+k_sorp*1e6)) * f_sorp  #* points
    cstot_css2 = np.array(s.getSolution_(7)).flatten()*bulkDensity_m3#* points
    cstot = cstot_cls + cstot_css1 + cstot_css2

    write_file_array("totcs",np.array([np.round(cstot_cls.sum()),
                    np.round(cstot_css1.sum()), 
                    np.round(cstot_css2.sum()), 
                    np.round(cstot.sum())]), "\t")
                    
    write_file_array("RF",RF)
    write_file_array("css1",cstot_css1/1e6)
    # write_file_array("cstot_css1Bis",cstot_css1Bis/1e6)
    write_file_array("css2",cstot_css2/1e6)
    write_file_array("cl",cstot_cls/1e6)
    write_file_array("cstot",cstot/1e6)
        
    contents = np.array([(cstot_cls * points).sum(),
                    (cstot_css1 * points).sum(), 
                    (cstot_css2 * points).sum(), 
                    (cstot * points).sum()])
    #print(contents)
    #print(min(cstot_cls))#contents[3]/contents0[3]*100-100)
    write_file_array("contents", contents/1e6)



