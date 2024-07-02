""" 
Jan's new scenario

with the classical sink

complicated MPI support (a non-mpi version of richards_cyl is needed, see script dumux3_nompi.sh)
"""
import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src")

from rosi_richards import RichardsSP  # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part

import matplotlib.pyplot as plt
import numpy as np
import time
import functional.van_genuchten as vg
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()
if False:
    import sys; 
    sys.path.append("../../python/modules/"); 
    sys.path.append("../../../CPlantBox/");  
    sys.path.append("../../../CPlantBox/src")
    #sys.path.append("../../../../../../DUMUXexudDune27/DUMUX/dumux-rosi/build-cmake/cpp/python_binding/");
    sys.path.append("../../../../../DUMUX/DUMUX/dumux-rosi/build-cmake/cpp/python_binding/");

    from rosi_richards import RichardsSP  # C++ part (Dumux binding)
    from richards import RichardsWrapper  # Python part
    


raise Exception
s = RichardsWrapper(RichardsSP())
raise Exception
s.initialize()
print('did init')
s.createGrid([-5., -5., -200.], [5., 5., 0.], [5,5,5])  # [cm]
p_mean_ = -15000
s.soil = [0.049, 0.352, 0.019, 4.887, 421.67]
s.vg_soil = vg.Parameters(s.soil) 
s.setHomogeneousIC(p_mean_, equilibrium = False)  # cm pressure head
s.setTopBC("noFlux")
s.setBotBC("noFlux") 
s.setVGParameters([s.soil])
s.initializeProblem()

pressureinit = s.getSolutionHead()
print('expected pressure',p_mean_)
try:
    print('expected theta', vg.water_content( p_mean_, s.vg_soil))
except:
    pass
print('obtained pressure', min(pressureinit))
#sat = s.getSaturation()
#print('obtained saturation', min(sat), max(sat),'obtained theta (from saturation)', 
#            min(sat)*s.vg_soil.theta_S, max(sat)*s.vg_soil.theta_S)
thetainit = s.getWaterContent_()
print('obtained theta ', min(thetainit), 'difference',min(thetainit)- vg.water_content( p_mean_, s.vg_soil) )

print('Phead from obtained saturation and theta')
#print(vg.pressure_head( min(sat)*s.vg_soil.theta_S, s.vg_soil))
print(vg.pressure_head(  min(thetainit), s.vg_soil))