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
import smallTest_ads_functions as stf
import scenario_setup
from scenario_setup import write_file_array, write_file_float

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
year = 2013
min_b = np.array([-5, -5, -10.] )
max_b =np.array( [5, 5, 0.] )
cell_number = np.array([2,2,2])
soil_type = "loam"
genotype = "WT"
comp = "phenolics"
usemoles = True
soil_ = scenario_setup.vg_SPP(0)
mode = "dumux_10c"
p_mean = -100
css1Function_ = 3
paramIdx = 5
dt = 1./3./24.
times = np.array([0.,dt])
s, soil = scenario_setup.create_soil_model(soil_type, year, soil_,#comp, 
                min_b, max_b, cell_number, demoType = mode, times = None, net_inf = None,
                usemoles = usemoles, dirResults = results_dir, p_mean_ = p_mean, 
                                         css1Function = css1Function_,
                                        paramIndx=paramIdx,
                                        noAds = False)



cell_volumes = comm.bcast(s.getCellVolumes() , root = 0) #cm3
s.numFluidComp = 2
water_content = comm.bcast( np.array(s.getWaterContent()),root= 0)  # theta per cell [1]
soil_water = np.multiply(water_content, cell_volumes)  # water per cell [cm3]
soil_solute_conz = comm.bcast(np.array([np.array(
    s.getSolution(i+1)) for i in range(s.numComp)]),
                                    root = 0) # mol/mol

s.solve(30/(24*60))
s.base.printParams()