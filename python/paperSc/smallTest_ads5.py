import sys; sys.path.append("../modules_fpit/"); sys.path.append("../../../CPlantBox/");  sys.path.append("../../build-cmake/cpp/python_binding/")
sys.path.append("../../../CPlantBox/src/python_modules")
sys.path.append("../../../CPlantBox/");
sys.path.append("../../../CPlantBox/src")

import matplotlib; matplotlib.use('agg')

from rosi_richards10c_cyl import Richards10CCylFoam  # C++ part (Dumux binding)

from richards import RichardsWrapper  # Python part
from richards_no_mpi import RichardsNoMPIWrapper  # Python part of cylindrcial 

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()
import scenario_setup as stf #
import smallTest_ads_functions as stf2
import scenario_setup
from scenario_setup import write_file_array, write_file_float

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

    
r_in = 0.02
r_out = 0.2104870824716315
paramIdx = 1640
stf2.setShape(s,r_in, r_out)
l = s.l
stf.setDefault(s)
stf.setSoilParam(s,paramIdx)
stf.getBiochemParam(s,paramIdx, noAds = False)
stf.setBiochemParam(s)
#ICZ = np.array([8.27842215e-05,0.01063999,4.77502731e-05 ])
stf.setIC3D(s,  paramIdx)

s.QCflowOut =np.array([0.,0.]) #np.array([ -5.136688687529069e-06 ,-2.910181366106239e-12]  )
s.QflowOut = 0.
s.win =  0.# -0.07004566384961891
s.exuds_in = 0.1
s.exudl_in = 0.
s.WatBC =  s.win* (2 * np.pi * s.r_in * s.l) + s.QflowOut
s.Qexud = (s.exuds_in+ s.exudl_in)* (2 * np.pi * s.r_in * s.l) +sum(s.QCflowOut)
#stf.getBC(s)
stf2.setBC(s)
css1Function=5
s.css1Function = css1Function
s.setParameter( "Soil.css1Function", str(s.css1Function))
s.setParameter( "Problem.segLength", str(s.l))

s.initializeProblem()

qOut = s.distributeSource(s.QflowOut , 0, s.l, s.numFluidComp)

valueTopBC = s.distributeSources(s.QCflowOut,
    np.array([nc+1 for nc in range(s.numFluidComp+1)]),
                                 l, s.numFluidComp)
s.setCriticalPressure(-15000)  # cm pressure head

times = np.array([0.,0.013888888888888888]) #,0.013888888888888888*2.,0.013888888888888888*3.,0.013888888888888888*4.] ) # days

s.ddt = 1.e-5
buTotCBefore = s.getTotCContent_each(s.l)
buWBefore_ =  s.getWaterVolumesCyl(l)
doReset = True 
print('buTotCBefore',buTotCBefore.sum(axis = 1))

for i, dt in enumerate(np.diff(times)):
    s.ddt =min( 1.e-5,s.ddt)

    #assert (sum(buWBefore_) +s.WatBC*dt  ) > 0
    #assert (sum(buTotCBefore) +s.Qexud*dt ) > 0
    if rank == 0:
        print("*****", "external time step", dt, 
              " d, simulation time", s.simTime, 
              "d, internal time step", s.ddt, "d")

    s.solve(dt, maxDt = 250./(24.*3600.))
    
    buTotCAfter = s.getTotCContent_each(s.l)
    buWAfter_ =  s.getWaterVolumesCyl(s.l)
    print('buTotCAfter',buTotCAfter.sum(axis = 1))
    print('buTotCBefore',sum(buTotCBefore.sum(axis = 1)),
            'buTotCAfter',sum(buTotCAfter.sum(axis = 1)),
            's.Qexud*dt',s.Qexud*dt,
            'diff',sum(buTotCAfter.sum(axis = 1)) - sum(buTotCBefore.sum(axis = 1)) - s.Qexud*dt)
    
    if False:
        print("content, before",sum(buTotCBefore),'after', sum(buTotCAfter), 'added',s.Qexud*dt)
        print(sum(buTotCBefore) +s.Qexud*dt - sum(buTotCAfter) )
        print("% error sollutes ",
              (sum(buTotCBefore) +s.Qexud*dt - sum(buTotCAfter))/(sum(buTotCBefore) +s.Qexud*dt)*100 )
        
        print("water",sum(buWBefore_), sum(buWAfter_), s.WatBC*dt)
        print(sum(buWBefore_) +s.WatBC*dt - sum(buWAfter_) )
        print("% error water ",(sum(buWBefore_) +s.WatBC*dt - sum(buWAfter_))/sum(buWAfter_)*100 )
        
    if doReset:
        s.reset()
    else:
        buTotCBefore = buTotCAfter
        buWBefore_  = buWAfter_


