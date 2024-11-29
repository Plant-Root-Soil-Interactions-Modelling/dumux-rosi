import sys;
#os.chdir('experimental/fixedPointIter2/scripts')
sys.path.append("../modules/");
sys.path.append("../inputDataPuptake/");
sys.path.append("../../../../CPlantBox/");
sys.path.append("../../../../CPlantBox/src")
sys.path.append("../../../build-cmake/cpp/python_binding/");

import matplotlib; matplotlib.use('agg')

import numpy as np
from numpy import array
import os
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()
# import smallTest_ads_functions as stf
import scenario_setup
import weatherFunctions
import helpfull
import ctypes
import tempfile

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

# outer time step (outside of fixed-point iteration loop)
dt = 20/60/24
dt_inner_init = 1/60/24 # dt
# min, max, objective number of iteration for the fixed-point iteration
minIter = 4 # empirical minimum number of loop to reduce error
k_iter = 30
targetIter= 5
# which functional modules to implement
doSoluteFlow = True # only water (False) or with solutes (True)
noAds = False # stop adsorption?
doPhloemFlow = False
doPhotosynthesis = False # photosynthesis-Transpiration (True) or just xylem flow (False)?
# when there is no transpiration, we use the plant wat. pot.
# at the beginning of the time step. Otherwise does not converge
beforeAtNight = True  
# static or growing organism
static_plant = False
# print debug messages in cyl3plant faile
mpiVerbose = False
# print debug messages in rhizo_mappedPlant
mpiVerboseInner = False
# how many files are printed. use 'False' in debug mode
# ATT: for short ismulations only
doMinimumPrint =  True
# use moles (mol) and not mass (g) in dumux
usemoles = True

paramIndx_ = 44
spellData = {'scenario':"none",'spellStart':0,'spellEnd':11, 'condition':"wet"}
weatherInit = weatherFunctions.weather(1.,dt, spellData)
s = scenario_setup.create_soil_model(usemoles = usemoles, 
                                           results_dir = results_dir, 
                                           p_mean_ = -weatherInit['p_mean'], 
                                    paramIndx=paramIndx_,
                                    noAds = noAds, doSoluteFlow = doSoluteFlow)


print('s.numSoluteComp',s.numSoluteComp)

#s.solve(dt_inner_init)
print('solution head after solve',s.getSolutionHead())

