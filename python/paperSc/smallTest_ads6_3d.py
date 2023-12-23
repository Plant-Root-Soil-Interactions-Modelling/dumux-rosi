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
min_b = [-5, -5, -10.] 
max_b = [5, 5, 0.] 
cell_number = [1,1,1]
soil_type = "loam"
genotype = "WT"
comp = "phenolics"
usemoles = True
soil_ = scenario_setup.vg_SPP(0)
mode = "dumux_10c"
p_mean = -100
css1Function_ = 0
paramIdx = 1640
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
soil_solute_content = comm.bcast(np.array([np.array(
    s.getContent(i+1, isDissolved = (i < s.numFluidComp))) for i in range(s.numComp+1)]),
                                    root = 0) # mol
soil_source_sol=np.full((s.numComp+1,s.numberOfCellsTot),0.)        
soil_source_sol[0] =  0.00205364944098248   #c_s                        
# soil_source_sol[6] =  0.0105992908371624    #C_ss2                      
soil_source_sol[8] =  0.0173749033084651      #C_SS1               
soil_fluxes_limited=np.full(s.numberOfCellsTot,0.)

soil_sources_limited = np.concatenate((np.array([soil_fluxes_limited]),soil_source_sol ))
soil_contents = np.concatenate((np.array([water_content]),soil_solute_content ))
buTotCBeforeEach = comm.bcast(s.getTotCContent_each(), root = 0) 
buTotCBeforeAll = buTotCBeforeEach.sum(axis = 0)
buTotCBefore = buTotCBeforeAll
print('buTotCBeforeEach',buTotCBeforeEach)

for idComp in range(s.numComp+1):#mol/day

    SSL = soil_sources_limited[idComp].copy()
    soil_contents_temp = soil_contents[idComp].copy()
    if idComp == 1:# add source of css1
        SSL += soil_sources_limited[s.numComp+1].copy()
        soil_contents_temp += soil_contents[s.numComp+1].copy()
        
    if (max(abs(SSL)) != 0.):
        SSL = np.maximum(SSL, -soil_contents_temp/dt)
        
        toAdd= np.maximum(0., -(soil_contents_temp/dt + SSL))
        SSL[np.where(toAdd>0.)] += toAdd[np.where(toAdd>0.)] #+ 1e-20
        
        k_limit_source3d = 0
        epsilon_source3d = 1e-25
        while (not (SSL*dt >= -soil_contents_temp).all()) and (k_limit_source3d <= 10):
            SSL[np.where(SSL*dt < -soil_contents_temp)] += epsilon_source3d
            epsilon_source3d *= 10
            k_limit_source3d += 1
        
        try:
            assert min(soil_contents_temp + SSL*dt) >=0.
        except:
            print(soil_sources_limited[idComp], SSL,dt,  min(soil_contents_temp + SSL*dt) )
            write_file_array("setsourceLim2_"+str(idComp), SSL, directory_ =results_dir) 
            raise Exception
        
        test_values = list(SSL)
        test_keys = np.array([i for i in range(len(test_values))])
        res = {}
        for key in test_keys:
            for value in test_values:
                res[key] = value
                test_values.remove(value)
                break                        
        write_file_float("setsource_"+str(idComp), res, directory_ =results_dir) 
        write_file_array("setsourceLim1_"+str(idComp),  soil_sources_limited[idComp], 
                         directory_ =results_dir) 
        write_file_array("setsourceLim2_"+str(idComp), SSL, directory_ =results_dir) 
        print('res',soil_sources_limited[idComp],res)
        s.setSource(res.copy(), eq_idx = idComp)  # [mol/day], in modules/richards.py
write_file_array("setsourceLim1_"+str(s.numComp+1),  soil_sources_limited[s.numComp+1], 
                         directory_ =results_dir) 
        
#outer_R_bc_sol =  sources_sol*0.          

s.ddt = 1.e-5

buWBefore_ =  (soil_water)
doReset = True 

s.Qexud = sum(soil_source_sol)
s.WatBC = sum(soil_fluxes_limited)
for i, dt in enumerate(np.diff(times)):
    s.ddt =min( 1.e-5,s.ddt)

    if rank == 0:
        print("*****", "external time step", dt, 
              " d, simulation time", s.simTime, 
              "d, internal time step", s.ddt, "d")

    s.solve(dt, maxDt = 250./(24.*3600.))

    outer_R_bc = np.transpose(s.getFlux_10c())
    bulkSoil_sources = np.transpose(s.getSource_10c()) 

    outer_R_bc = np.transpose(s.getFlux_10c())
    bulkSoil_sources = np.transpose(s.getSource_10c()) #buTotCBeforeEach[-1], buTotCAfterEach[-1]))
    outer_R_bc_wat = outer_R_bc[0]# [cm3] #new_soil_water - soil_water - soil_fluxes_limited *dt # change in water per cell [cm3]
    sources_wat =  bulkSoil_sources[0]    
    
    outer_R_bc_sol = outer_R_bc[1:]#np.array([soil_solute_content_new[i] - soil_solute_content[i] - soil_source_sol[i]*dt for i in range(r.numComp)])# mol
    sources_sol =  bulkSoil_sources[1:]
    
    buTotCAfterEach = comm.bcast(s.getTotCContent_each(), root = 0) 
    print('buTotCAfterEach',buTotCAfterEach)
    print('buTotCBefore',sum(buTotCBeforeEach),
            'buTotCAfter',sum(buTotCAfterEach),
            's.Qexud*dt',s.Qexud*dt,
            'diff',sum(buTotCAfterEach) - sum(buTotCBeforeEach) - s.Qexud*dt)
    
    buTotCAfterAll = buTotCAfterEach.sum(axis = 0)
    buTotCAfter = buTotCAfterAll
    water_content =comm.bcast( np.array(s.getWaterContent()), root = 0) 
    new_soil_water = np.multiply(water_content, cell_volumes)  # calculate net flux
            
            
    buWAfter_ = (new_soil_water)
    
    print("content, before",sum(buTotCBefore),'after', sum(buTotCAfter), 'added',s.Qexud*dt)
    print(sum(buTotCBefore) +s.Qexud*dt - sum(buTotCAfter) )
    print("% error sollutes ",
          (sum(buTotCBefore) +s.Qexud*dt - sum(buTotCAfter))/(sum(buTotCBefore) +s.Qexud*dt)*100 )
    
    print("water",sum(buWBefore_), sum(buWAfter_), s.WatBC*dt)
    print(sum(buWBefore_) +s.WatBC*dt - sum(buWAfter_) )
    print("% error water ",(sum(buWBefore_) +s.WatBC*dt - sum(buWAfter_))/sum(buWAfter_)*100 )
    

    s.bulkMassErrorWaterAll_real = new_soil_water - (soil_water + sources_wat + outer_R_bc_wat)
    s.bulkMassErrorWaterAll_abs = abs(s.bulkMassErrorWaterAll_real)
    s.bulkMassErrorWaterAll_rel = abs(s.bulkMassErrorWaterAll_abs /new_soil_water )*100
    
    s.bulkMassCError1dsAll_real = buTotCAfterEach - ( buTotCBeforeEach + sources_sol + outer_R_bc_sol)
    print('s.bulkMassCError1dsAll_real_A',s.bulkMassCError1dsAll_real)
    s.bulkMassCError1dsAll_real[0] += s.bulkMassCError1dsAll_real[-1]#put together CS and CSS1
    s.bulkMassCError1dsAll_real[-1][:] = 0.
    print('s.bulkMassCError1dsAll_real_B',s.bulkMassCError1dsAll_real)
    
    if doReset:
        s.reset()
    else:
        buTotCBefore = buTotCAfter
        buWBefore_  = buWAfter_
        buTotCBeforeEach = buTotCAfterEach




print('after solve: s.getSolutionHead()',np.array(s.getSolutionHead()))
for ncomp in range(s.numComp):
    print('after solve: s.getSolution(ncomp + 1)',ncomp + 1,
          np.array(s.getSolution(ncomp + 1)).flatten())
