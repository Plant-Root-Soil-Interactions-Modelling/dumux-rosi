import sys; sys.path.append("../modules_fpit/"); sys.path.append("../../../CPlantBox/");  sys.path.append("../../build-cmake/cpp/python_binding/")
sys.path.append("../../../CPlantBox/src/python_modules")

import matplotlib; matplotlib.use('agg')

from rosi_richards10c_cyl import Richards10CCylFoam  # C++ part (Dumux binding)

from richards import RichardsWrapper  # Python part
from richards_no_mpi import RichardsNoMPIWrapper  # Python part of cylindrcial 

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()
import smallTest_ads_functions as stf

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
#s = RichardsWrapper(Richards10CCylFoam(), usemoles)
s=RichardsNoMPIWrapper(Richards10CCylFoam(), usemoles)

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
    #print("CSW", np.array(cyl.getSolution_(1)).flatten()*cyl.molarDensityWat)
    #print("css1Test",cyl.CSSmax * (C_S_W/(C_S_W+cyl.k_sorp)) * cyl.f_sorp*vols,  (C_S_W/(C_S_W+cyl.k_sorp)) )
    print("css1", sum(css1*vols))
    assert len(np.array(totC).shape)==1
    return totC
    
r_in = 0.02
r_out = 0.7981351842686684
paramIdx = 1646
stf.setShape(s,r_in, r_out)
l = s.l
stf.setDefault(s)
stf.setSoilParam(s,paramIdx)
stf.getBiochemParam(s,paramIdx)
stf.setBiochemParam(s)
stf.setIC(s, -108.5, paramIdx)
stf.getBC(s)
stf.setBC(s)
    
s.setParameter("Newton.MaxRelativeShift", "1e-5")

s.setParameter("Problem.EnableGravity", "false")
s.setParameter("Problem.verbose", "0")
s.setParameter("Newton.MaxSteps", "30")
s.setParameter("Newton.ResidualReduction", "1e-10")
s.setParameter("Newton.SatisfyResidualAndShiftCriterion", "true")
s.setParameter("Newton.TargetSteps", "10")
s.setParameter("Newton.UseLineSearch", "false")
s.setParameter("Newton.EnablePartialReassembly", "true")
s.setParameter("Newton.EnableAbsoluteResidualCriterion", "true")
s.setParameter("Newton.MaxAbsoluteResidual", "1.e-10")
s.setParameter("Newton.EnableChop", "true")
s.setParameter("Newton.EnableResidualCriterion", "true")
s.setParameter("Newton.EnableShiftCriterion", "true")
s.setParameter("Newton.MaxAbsoluteResidual", "1e-10")

s.initializeProblem()

qOut = s.distributeSource(s.QflowOut , 0, s.l, s.numFluidComp)
s.setCriticalPressure(-15000)  # cm pressure head

times = np.array([0.,0.013888888888888888,0.013888888888888888*2.] ) # days

s.ddt = 1.e-5
buTotCBefore = getTotCContent(0, s,s.l, init=True)
buWBefore_ =  s.getWaterVolumesCyl(l)
doReset = True 
for i, dt in enumerate(np.diff(times)):
    s.ddt =min( 1.e-5,s.ddt)

    assert (sum(buWBefore_) +s.WatBC*dt  ) > 0
    assert (sum(buTotCBefore) +s.Qexud*dt ) > 0
    if rank == 0:
        print("*****", "external time step", dt, 
              " d, simulation time", s.simTime, 
              "d, internal time step", s.ddt, "d")

    s.solve(dt, maxDt = 1)
    
    buTotCAfter = getTotCContent(0, s,s.l)
    buWAfter_ =  s.getWaterVolumesCyl(s.l)
    
    print("content, before",sum(buTotCBefore),'after', sum(buTotCAfter), 'added',s.Qexud*dt)
    print(sum(buTotCBefore) +s.Qexud*dt - sum(buTotCAfter) )
    print("% error sollutes ",(sum(buTotCBefore) +s.Qexud*dt - sum(buTotCAfter))/(sum(buTotCBefore) +s.Qexud*dt)*100 )
    
    print("water",sum(buWBefore_), sum(buWAfter_), s.WatBC*dt)
    print(sum(buWBefore_) +s.WatBC*dt - sum(buWAfter_) )
    print("% error water ",(sum(buWBefore_) +s.WatBC*dt - sum(buWAfter_))/sum(buWAfter_)*100 )
    
    if doReset:
        s.reset()
    else:
        buTotCBefore = buTotCAfter
        buWBefore_  = buWAfter_



