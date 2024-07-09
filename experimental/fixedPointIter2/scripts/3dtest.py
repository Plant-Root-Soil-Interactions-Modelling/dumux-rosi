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


def sendSource2dumux(s, SSL, idComp):
    # convert array to dictionnary
    test_values = list(SSL)
    test_keys = np.array([i for i in range(len(test_values))])
    res = {}
    for key in test_keys:
        for value in test_values:
            res[key] = value
            test_values.remove(value)
            break
    # send to dumux
    print('to set source', idComp, res)
    s.setSource(res.copy(), eq_idx=idComp)  # [mol/day], in modules/richards.py

SSLs={  0 : [-1.63454452e+01 ,-4.86715484e+00, -1.83306654e+01, -1.08568187e-02,
  7.35655294e-04, -8.81541736e-04, -8.01383865e-01, -1.04270112e+00,
 -3.12108695e-01, -9.62598348e-02, -1.38178350e-02, -8.60528561e-02,
 -3.35928275e-01, -2.33806565e-03, -7.50076951e-02, -3.31687877e+00,
 -2.81513572e+00, -1.00766871e+00, -1.72151391e+00, -1.33918174e+00,
 -1.80016800e+00, -6.46870046e+00, -6.93422203e+00, -3.59652243e+00,
 -4.88997450e+00, -8.05473910e-01, -1.07733677e-01],
        1 : [-1.77723125e-07, -1.29168650e-07, -1.08133850e-07, -4.82988423e-08,
 -5.14958944e-08, -1.38697326e-07, -1.80575566e-07, -2.26362494e-07,
 -2.39602582e-07, -1.13959276e-07, -2.28595100e-08, -9.21221533e-08,
 -6.42909176e-08, -5.43552068e-09, -4.19748465e-08, -1.59564690e-07,
 -1.24683058e-07, -1.02375676e-07, -6.14917403e-08, -2.35603801e-08,
 -6.18693313e-08, -7.11324113e-08, -6.25536627e-08 ,-7.65804191e-08,
 -1.12330749e-07, -8.98053837e-09, -1.27902545e-09]
        }
for idComp in SSLs:
    print(idComp, SSLs[idComp])
    sendSource2dumux(s, SSLs[idComp], idComp)

s.solve(dt_inner_init)

