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
paramIndx_ = 1476
    
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

p_mean = -100.#pd.read_csv(data_dir + '/TraiRhizoparam_ordered.csv').iloc[-1].to_numpy()
#ICcc = np.array([min(pd.read_csv().iloc[-1].to_numpy()) if (ncomp < 3) else max(pd.read_csv().iloc[-1].to_numpy()) for ncomp in range(1,9)]
##
s, soil = scenario.create_soil_model(soil_type, year, soil_,#comp, 
                min_b, max_b, cell_number, demoType = mode, times = None, net_inf = None,
                usemoles = usemoles, dirResults = results_dir, p_mean_ = p_mean, 
                                         css1Function = css1Function_,
                                        paramIndx=paramIndx_,
                                        noAds = noAds)
##
#reset wat and sol content according to the last values in paramIndxnoAds14760dx_2dumux_10c_9.5to10.5_20mn_0s_1_100
# also get source
thetaInit = pd.read_csv(data_dir +'theta.txt').iloc[-1].to_numpy()

pheadinit_cm =  np.array([vg.pressure_head( tlo, s.vg_soil) for tlo in thetaInit])
g_ = 9.81 #9.80665; // cm / s^2 (for type conversions)
rho_ = 1.e3 #kg / m^3 (for type conversions)
pRef_ = 1.e5 #Pa 101300; // Pa
pheadinit_hPa = pheadinit_cm/100. *rho_*g_ +pRef_
s.base.setSolution(pheadinit_hPa,0 )
for eqIdx in range(1,9):
    solInit = pd.read_csv(data_dir +'Soil_solute_conc'+str(eqIdx)+'.txt').iloc[-1].to_numpy()
    s.base.setSolution(solInit,eqIdx )
    
# source
soil_solute_content=comm.bcast(np.array([
    np.array( s.getContent(i+1, isDissolved = (i < s.numFluidComp))) for i in range(s.numComp)]),root = 0)
water_content = comm.bcast( np.array(s.getWaterContent()),root= 0)  # theta per cell [1]

soil_contents = np.concatenate((np.array([water_content]),soil_solute_content )) 
soil_source = np.full((s.numComp+1,5*5*20),0. )
if True:
    for idComp in range(s.numComp+1):
        if idComp != 7:
            try:
                soil_source[idComp] = pd.read_csv(data_dir +'setsourceLim2_'+str(idComp)+'.csv').iloc[-1].to_numpy()

                SSL = soil_source[idComp].copy()
                toAdd= np.maximum(0., -(soil_contents[idComp]/dt_ + SSL))
                print('toAdd',  min(toAdd), max(toAdd), np.where(toAdd>0.), toAdd[np.where(toAdd>0.)], 
                     -(soil_contents[idComp]/dt_ + SSL)[np.where(toAdd>0.)], soil_contents[idComp][np.where(toAdd>0.)])
                print('source: ',idComp,' min(soil_contents[idComp] ', min(soil_contents[idComp] ),
                      min(soil_contents[idComp] + SSL*dt_),min(soil_contents[idComp] +  soil_source[idComp]*dt_), )
                print('test1',np.where(SSL*dt_ < -soil_contents[idComp]), SSL[np.where(SSL*dt_ < -soil_contents[idComp])] ,
                      -soil_contents[idComp][np.where(SSL*dt_ < -soil_contents[idComp])]/dt_)
                SSL[np.where(toAdd>0.)] += toAdd[np.where(toAdd>0.)] #+ 1e-20
                print('test2',np.where(SSL*dt_ < -soil_contents[idComp]), SSL[np.where(SSL*dt_ < -soil_contents[idComp])] ,
                       -(soil_contents[idComp]/dt_ + SSL)[np.where(toAdd>0.)])
                epsilon = 1e-25
                while not (SSL*dt_ >= -soil_contents[idComp]).all():
                    SSL[np.where(SSL*dt_ < -soil_contents[idComp])] += epsilon
                    epsilon *= 10
                print('test3',np.where(SSL*dt_ < -soil_contents[idComp]), SSL[np.where(SSL*dt_ < -soil_contents[idComp])] ,
                       -(soil_contents[idComp]/dt_ + SSL)[np.where(toAdd>0.)])
                #SSL[np.where(SSL*dt_ < -soil_contents[idComp])] += np.maximumsoil_contents[idComp][np.where(SSL*dt_ < -soil_contents[idComp])]/dt_ + 
                #print('test1',np.where(SSL*dt_ < -soil_contents[idComp]), SSL[np.where(SSL*dt_ < -soil_contents[idComp])] ,
                #      -soil_contents[idComp][np.where(SSL*dt_ < -soil_contents[idComp])]/dt_)
                assert (SSL*dt_ >= -soil_contents[idComp]).all()

                test_values = list(SSL.copy())
                test_keys = np.array([i for i in range(len(test_values))])
                res = {}
                for key in test_keys:
                    for value in test_values:
                        res[key] = value
                        test_values.remove(value)
                        break                        
                scenario.write_file_float("setsource_"+str(idComp), res, directory_ =results_dir) 
                scenario.write_file_array("setsourceLim2_"+str(idComp), SSL, directory_ =results_dir, fileType =".csv") 
                s.setSource(res.copy(), eq_idx = idComp)  # [mol/day], in modules/richards.py
            except  Exception as err:
                print(idComp, f"{err=}, {type(err)=}", 'no source set for element')
                raise Exception

s.Qexud = sum(soil_source[1:][:].reshape(-1)) #sum changes of C
s.WatBC = sum(soil_source[0])#sum changes of W
##

s.ddt = 1.e-5
buTotCBefore = s.getTotCContent()
buWBefore_ =  s.getWaterVolumes()
doReset = True 

s.setParameter("Flux.UpwindWeight", "1")
phead = np.array(s.getSolutionHead())
scenario.write_file_array("solScompStart", phead, directory_ =results_dir, fileType =".csv")
print('before solve: min, max vals',min(phead), max(phead))
for ncomp in range(s.numComp):
    solScomp =  np.array(s.getSolution(ncomp + 1)).flatten()
    print('before solve: s.getSolution(ncomp + 1)',ncomp + 1,
         min(solScomp), max(solScomp))
    scenario.write_file_array("solScompStart", solScomp, directory_ =results_dir, fileType =".csv") 

times = np.array([ step_* dt_ for step_ in range(2)] ) # days
for i, dt in enumerate(np.diff(times)):
    s.ddt =min( 1.e-5,s.ddt)

    assert (sum(buWBefore_) +s.WatBC*dt  ) > 0
    assert (sum(buTotCBefore) +s.Qexud*dt ) > 0
    if rank == 0:
        print("*****", "external time step", dt, 
              " d, simulation time", s.simTime, 
              "d, internal time step", s.ddt, "d")

    s.solve(dt, maxDt = 1)
    
    phead = np.array(s.getSolutionHead())
    print('after solve: min, max vals',min(phead), max(phead))
    scenario.write_file_array("solScompEnd", phead, directory_ =results_dir, fileType =".csv") 
    for ncomp in range(s.numComp):
        solScomp =  np.array(s.getSolution(ncomp + 1)).flatten()
        print('after solve: s.getSolution(ncomp + 1)',ncomp + 1,
             min(solScomp), max(solScomp))
        scenario.write_file_array("solScompEnd", solScomp, directory_ =results_dir, fileType =".csv") 
        
    buTotCAfter = s.getTotCContent()
    buWAfter_ =  s.getWaterVolumes()
    
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




