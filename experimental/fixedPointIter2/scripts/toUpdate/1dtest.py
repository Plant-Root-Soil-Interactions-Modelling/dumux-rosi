import sys;
import os

#os.chdir('experimental/fixedPointIter2/scripts')
sys.path.append("../modules/");
sys.path.append("../inputDataTraiRhizo/");
sys.path.append("../../../../CPlantBox/");
sys.path.append("../../../../CPlantBox/src")
sys.path.append("../../../build-cmake/cpp/python_binding/");

import matplotlib; matplotlib.use('agg')

from rosi_richards10c_cyl import RichardsNCCylFoam  # C++ part (Dumux binding)

from richards import RichardsWrapper  # Python part
from richards_no_mpi import RichardsNoMPIWrapper  # Python part of cylindrcial

import matplotlib.pyplot as plt
import numpy as np
from numpy import array
import pandas as pd
import os
import scenario_setup

""" 
Cylindrical 1D model, diffusion only (DuMux), Michaelis Menten

everything scripted, no input file needed, also works parallel with mpiexec
"""

# directory where the results will be printed
results_dir="./results/1dtest/"
if not os.path.exists(results_dir):
    os.makedirs(results_dir)
else:
    test = os.listdir(results_dir)
    for item in test:
        try:
            os.remove(results_dir+item)
        except:
            pass


                          

def initialize_dumux_nc_(  gId=0, a_in=0.02,
                a_out=0.2 ,seg_length=1.0 ,
                x=[-100] ,   # cm
                cAll = [0,0,0,0,0,0,0,0 ],             # mol/mol scv
                                         Cells = [],NC = 10,
                                         logbase = 0.5):                                   # cm
    verbose = False
    print('tocyl')
    lId =gId
    
    if a_in < a_out:
    
        cyl = RichardsNoMPIWrapper(RichardsNCCylFoam(), True)  # only works for RichardsCylFoam compiled without MPI

        cyl.initialize(verbose = False) # No parameter file found. Continuing without parameter file.
        soilVG = [0.08, 0.43, 0.04, 1.6, 50]
        cyl.setVGParameters([soilVG])
        lb =  logbase
        
        points = np.logspace(np.log(a_in) / np.log(lb), np.log(a_out) / np.log(lb), 
                              NC, base = lb)
        
        cyl.createGrid1d(points)# cm
            
        cyl.setParameter("SpatialParams.Temperature","293.15") # todo: redefine at each time step
        cyl.setParameter("Soil.BC.dzScaling", "1")
        cyl.setParameter( "Soil.css1Function", str(9))
        cyl.setParameter("Problem.verbose", "0")
        cyl.setParameter("Problem.reactionExclusive", "0")    
        cyl.setParameter("Soil.CriticalPressure", str(-15000))
        cyl.seg_length = seg_length
        cyl.setParameter("Problem.segLength", str(seg_length))   # cm
        cyl.l = seg_length   
        cyl.setParameter( "Soil.Grid.Cells",str( NC-1)) # -1 to go from vertices to cell (dof)
        if verbose:
            print("Soil.IC.P", cyl.dumux_str(x))
        cyl.setParameter("Soil.IC.P", cyl.dumux_str(x))# cm
        
        #default: no flux
        cyl.setInnerBC("fluxCyl", 0.)  # [cm/day] #Y 0 pressure?
        cyl.setOuterBC("fluxCyl", 0.)
        
        
        for j in range( 1, cyl.numComp):
            cyl.setParameter( "Soil.BC.Bot.C"+str(j)+"Type", str(3))
            cyl.setParameter( "Soil.BC.Top.C"+str(j)+"Type", str(3))
            cyl.setParameter( "Soil.BC.Bot.C"+str(j)+"Value", str(0)) 
            cyl.setParameter( "Soil.BC.Top.C"+str(j)+"Value", str(0 ))

        # mol/cm2/day
        for j in range( 1, cyl.numComp):       
            cyl.setParameter("Soil.IC.C"+str(j),"0" ) 


        cyl.maxDt = 250/(3600*24) # soilModel.maxDt_1DS
        cyl.setParameter("Problem.verbose", "0")
        cyl.setParameter("Newton.Verbosity", "0")
        cyl.initializeProblem(maxDt=cyl.maxDt ) # ewton solver configured with the following options and parameters:


        cyl.setCriticalPressure(-15000)  # cm pressure head
        
        
        
        return cyl
    else:
        print("RhizoMappedSegments.initialize_dumux_: Warning, segment {:g} might not be in domain, radii [{:g}, {:g}] cm".format(i, a_in, a_out))
        return []
        


cyl = initialize_dumux_nc_()
print(cyl.numSoluteComp)
cyl.solve(5)
print(cyl.getSolutionHead())