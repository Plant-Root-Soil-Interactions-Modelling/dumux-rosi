import sys; sys.path.append("../modules_fpit/"); sys.path.append("../../../CPlantBox/");  sys.path.append("../../build-cmake/cpp/python_binding/")
sys.path.append("../../../CPlantBox/src/")

import matplotlib; matplotlib.use('agg')

from rosi_richards10c_cyl import Richards10CCylFoam  # C++ part (Dumux binding)

from richards import RichardsWrapper  # Python part
from richards_no_mpi import RichardsNoMPIWrapper  # Python part of cylindrcial 

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()
import scenario_setup as scenario
import functional.van_genuchten as vg
import scenario_setup as stf

""" 
Cylindrical 1D model, diffusion only (DuMux), Michaelis Menten

everything scripted, no input file needed, also works parallel with mpiexec
"""

usemoles = True
#s = RichardsWrapper(Richards10CCylFoam(), usemoles)
s=RichardsNoMPIWrapper(Richards10CCylFoam(), usemoles)

s.initialize()


# DUMUX_NUM_THREADS=1 mpirun -n 1 python3 XcGrowth.py 9.5 dumux_10c 10.5 $1 noAds
initsim = 9.5
mode = "dumux_10c"
simMax = 10.5
extraName =  'noAds'
paramIdx = 98
    
dt_ = 1/3/24
k_iter = 20
l_ks =  "dx_2"#"root", "dx", "dx_2"
organism = "plant"# "RS"#
weightBefore = False
SRIBefore = False
beforeAtNight = True
adaptRSI_  = False
static_plant = False
useOuterFluxCyl_w = False
useOuterFluxCyl_sol = False
css1Function_ = 0
lightType =""#+- "nolight" # or ""
mpiVerbose = False
noAds = (extraName == 'noAds')

results_dir="./results/errorSoil/"
data_dir = "./results/paramIndxnoAds14760dx_2dumux_10c_9.5to10.5_20mn_0s_1_100_saved/"
comm.barrier()
print('results_dir','DUMUXexudDune27/DUMUX/dumux-rosi/python/paperSc/',results_dir, flush = True)
comm.barrier()
if rank == 0:
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)
    else:
        test = os.listdir(results_dir)
        for item in test:
            try:
                os.remove(results_dir+item)
            except:
                pass
"""scenario"""
year = 2013
soil_type = "loam"
genotype = "WT"
comp = "phenolics"
usemoles = True
""" parameters   """
soil_ = scenario.vg_SPP(0)

min_b = [-5, -5, -10.] 
max_b = [5, 5, 0.] 
cell_number = [5,5,20]


""" rhizosphere model parameters """
recreateComsol = False

periodic = False
nc = 10

logbase = 0.5  # according to Mai et al. (2019)

""" initialize """

##
#reset wat and sol content according to the last values in paramIndxnoAds14760dx_2dumux_10c_9.5to10.5_20mn_0s_1_100
# also get source

p_mean = -111.5

##
##### 1D
cell2rhizoId =0
demoType = mode
dirResults = results_dir
p_mean_ = p_mean
css1Function = css1Function_
paramIndx=paramIdx; 

s1d=RichardsNoMPIWrapper(Richards10CCylFoam(), usemoles)      
s1d.dirResults = dirResults

#volCel = cell_volumes[cell2rhizoId]
r_in = 0.02
length_ = 1.
r_out = 0.111388368129964#np.sqrt(volCel/(np.pi*length_)+r_in*r_in)#0.2
# (r_out*r_out-r_in*r_in)*np.pi*length_
s1d.initialize()
#
#print(s1d.dimWorld,s.dimWorld)
stf.setShape1D(s1d,r_in, r_out,length = length_,nCells = 10,doLogarithmic=True)

stf.setDefault(s1d)
s1d.setParameter("Problem.reactionExclusive", "0")  
#s1d.setParameter("Newton.MaxRelativeShift","1e-10")
stf.setSoilParam(s1d,paramIndx)
stf.getBiochemParam(s1d,paramIndx,noAds)
stf.setBiochemParam(s1d)
s1d.win = 0.
stf.setIC(s1d,paramIndx )
    
s1d, __ = stf.setupOther(s1d, css1Function, p_mean_)
#reset wat and sol content according to the last values in paramIndxnoAds14760dx_2dumux_10c_9.5to10.5_20mn_0s_1_100
# also get source


pheadinit_cm = -75

pheadinit_Pa = s.to_pa(pheadinit_cm)
print('pheadinit_cm_all',s1d.to_head(pheadinit_Pa), pheadinit_cm, s1d.getSolutionHead())
pheadinit_PaCyl = np.full(s1d.getSolutionHead().shape,pheadinit_Pa)
s1d.base.setSolution(pheadinit_PaCyl,0 )
print( s1d.getSolutionHead())
   
cyl = s1d
    
pheadOld = cyl.getSolutionHead()
cyl_cell_volumes = cyl.getCellVolumesCyl()  #cm3 scv                       
oldTheta = cyl.getWaterContent()

nc_content = np.array([cyl.getContentCyl(nc+1, nc < 2)  for nc in range(s1d.numFluidComp)])# mol
cyl.base.setSolution(pheadinit_PaCyl,0 )

newTheta = cyl.getWaterContent()
newWatMol = (cyl.getWaterContent() * cyl_cell_volumes) * (1/1e6) * s1d.molarDensityWat_m3 
nc_molFr =np.array( [nc_c/newWatMol for nc_c in nc_content])
for nc in range(s1d.numFluidComp):
    cyl.base.setSolution(nc_molFr[nc],nc+1 )
nc_content_new = np.array([cyl.getContentCyl(nc+1, nc < 2)  for nc in range(s1d.numFluidComp)])# mol

try:
    assert (abs(nc_content.reshape(-1) - nc_content_new.reshape(-1)) <= np.minimum(abs(nc_content.reshape(-1)),abs( nc_content_new.reshape(-1)))*1e-6  ).all()
    print('sucess')
    print('nc_content',nc_content)
    print('nc_content_new',nc_content_new)
    print(nc_content.reshape(-1) - nc_content_new.reshape(-1))
except:
    print('the solute content changed')
    print('nc_content_old\t',nc_content)
    print('nc_content_new\t',nc_content_new)
    print(nc_content.reshape(-1) - nc_content_new.reshape(-1))
    print('vols',cyl_cell_volumes)
    print('points',cyl.getPoints())
    print(np.minimum(abs(nc_content.reshape(-1)),abs( nc_content_new.reshape(-1)))*1e-6)
    print(abs(nc_content.reshape(-1) == nc_content_new.reshape(-1)) <= np.minimum(abs(nc_content.reshape(-1)),abs( nc_content_new.reshape(-1)))*1e-6 )
    #print('newWatMol',newWatMol,'nc_molFr',nc_molFr)
    
