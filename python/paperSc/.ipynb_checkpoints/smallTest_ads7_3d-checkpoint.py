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
times = np.array([0.,dt]) #,0.013888888888888888*2.,0.013888888888888888*3.,0.013888888888888888*4.] ) # days
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
# * 
soil_solute_conz[:2] *= s.molarDensityWat_cm3
soil_solute_conz[2:] *= s.bulkDensity_m3/1e6

CSS1 = soil_solute_conz[-1] 
CSW = soil_solute_conz[0] 
CSS2 = soil_solute_conz[6] 
CSWinit = 4.707833666666667e-06 #mol/cm3
print('CSS1, CSW, CSS2',CSS1, CSW, CSS2, soil_solute_conz)
print('CSS1[0], CSW[0], CSS2[0]',CSS1[0], CSW[0], CSS2[0])
print(s.computeDtCSS2(CSS1[0]* s.cm3_per_m3, CSW[0]* s.cm3_per_m3, CSS2[0]* s.cm3_per_m3), "mol/cm3/d")
print(s.computeDtCSS2(CSS1[0]* s.cm3_per_m3, CSW[0]* s.cm3_per_m3, CSS2[0]* s.cm3_per_m3)*1e6/(24*3600), "mol/m3/s")
print('initcss2',s.computeInitCSS2(CSS1[0]* s.cm3_per_m3, CSW[0]* s.cm3_per_m3))
print(s.computeDtCSS2(CSS1[0]* s.cm3_per_m3, CSW[0]* s.cm3_per_m3,
                      s.computeInitCSS2(CSS1[0]* s.cm3_per_m3, CSW[0]* s.cm3_per_m3)* s.cm3_per_m3))
#print(s.DtCSS2(CSS1[0]* s.cm3_per_m3, CSW[0]* s.cm3_per_m3, CSS2[0]* s.cm3_per_m3),"mol/m3/s")

#print(s.DtCSS2(CSS1[0]* s.cm3_per_m3, CSW[0]* s.cm3_per_m3, CSS2[0]* s.cm3_per_m3),"mol/m3/s")
#print(s.DtCSS2(CSS1[0]* s.cm3_per_m3, CSWinit* s.cm3_per_m3, CSS2[0]* s.cm3_per_m3),"mol/m3/s")

print(s.k_sorp)
print(s.CSSmax * s.cm3_per_m3 , 
      ( CSW[0]* s.cm3_per_m3/( CSW[0]* s.cm3_per_m3+s.k_sorp* s.cm3_per_m3)) , 
      ( CSWinit/( CSWinit+s.k_sorp)), 
      CSS2[0]* s.cm3_per_m3)
print('CSWinit',CSWinit,'CSW',CSW[0])
print(( CSW[0]* s.cm3_per_m3/( CSW[0]* s.cm3_per_m3+s.k_sorp* s.cm3_per_m3))-0.5306318020538945)
print(( CSWinit* s.cm3_per_m3/(CSWinit* s.cm3_per_m3+s.k_sorp* s.cm3_per_m3))-0.5306318020538945)
#css2init = 
print(s.CSS2_init-CSS2[0])
# s.getCSS2Init =  lambda CSW: s.CSSmax * (CSW/(CSW+ s.k_sorp))#mol C/ cm3 scv
# DtCSS2 =  lambda CSS1, CSW, CSS2: s.alpha /(24.*60.*60.) * (s.CSSmax * cm3_per_m3 * (CSW/(CSW+s.k_sorp * cm3_per_m3))- CSS2)
# (s.CSSmax * cm3_per_m3 * (CSW/(CSW+s.k_sorp * cm3_per_m3)) - s.CSSmax * (CSW/(CSW+ s.k_sorp)) * cm3_per_m3)

s.solve(30/(24*60))
print(s.base.getReac_CSS2())

s.base.printParams()