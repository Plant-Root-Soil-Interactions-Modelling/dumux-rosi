import sys;
import os
#os.chdir('experimental/fixedPointIter2/scripts')
sys.path.append("../modules/");
sys.path.append("../inputDataPuptake/");
sys.path.append("../../../../CPlantBox/");
sys.path.append("../../../../CPlantBox/src")
sys.path.append("../../../build-cmake/cpp/python_binding/");

import matplotlib; matplotlib.use('agg')

import numpy as np
from numpy import array
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()
# import smallTest_ads_functions as stf
import scenario_setup
import weatherFunctions
import helpfull
import ctypes
import tempfile

from rosi_richards10c import RichardsNCSPILU as RichardsNCSP  # C++ part (Dumux binding), macroscopic soil model
from richards import RichardsWrapper  # Python part, macroscopic soil model

# directory where the results will be printed
results_dir="./results/3dtest/"
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
s = RichardsWrapper(RichardsNCSP(), usemoles)  # water and N solute

s.initialize(np.array(["10c3d.input"])) 


boundsMin = np.array([-0.015 , -0.06 , -0.4 ])
boundsMax = np.array([0.015 , 0.06 , 0])
numberOfCells = np.array([3,12,40])
s.base.createGrid(boundsMin,boundsMax, numberOfCells, True) 
s.base.Parameters_init("10c3d.input" )
maxDt = 250#/(24*3600)
s.base.initializeProblem(maxDt);
s.base.solve(100)#/(24*3600)