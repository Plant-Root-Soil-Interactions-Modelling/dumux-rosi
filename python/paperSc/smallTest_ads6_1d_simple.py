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
import scenario_setup as stf #
import smallTest_ads_functions as stf2
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
cell_number = [2,1,1]
soil_type = "loam"
genotype = "WT"
comp = "phenolics"
usemoles = True
soil_ = scenario_setup.vg_SPP(0)
mode = "dumux_10c"
p_mean = -100
css1Function_ = 8
paramIndx = 1640
dt = 1./3./24.
times = np.array([0.,dt,dt*2])  # days
#s3d, soil = scenario_setup.create_soil_model(soil_type, year, soil_,#comp, 
#                min_b, max_b, cell_number, 
demoType = mode
dirResults = results_dir
p_mean_ = p_mean
css1Function = css1Function_
noAds = True

s=RichardsNoMPIWrapper(Richards10CCylFoam(), usemoles)

   
r_in = 0.02
r_out = 0.2
paramIdx = 1640
s.initialize()

stf.setShape1D(s,r_in, r_out,length = 1.,nCells = 10,doLogarithmic=False)

stf.setDefault(s)
s.setParameter("Newton.MaxSteps", "200")
s.setParameter("Newton.MaxTimeStepDivisions", "100")
#s.setParameter("Newton.EnableResidualCriterion", "true") # sometimes helps, sometimes makes things worse
#s.setParameter("Newton.EnableAbsoluteResidualCriterion", "true")
#s.setParameter("Newton.SatisfyResidualAndShiftCriterion", "true")
s.setParameter("Problem.reactionExclusive", "0")  
s.setParameter("Newton.MaxRelativeShift","1e-10")
stf.setSoilParam(s,paramIndx)
stf.getBiochemParam(s,paramIndx,noAds)
stf.setBiochemParam(s)
stf.setIC(s,paramIndx)
s.win = 0.; s.exudl_in = 0.
s.exuds_in =1e-6#02;#( 0.00205364944098248+0.0173749033084651)
stf.setBC1D(s)

s.initializeProblem()
print('get length')
l = s.segLength
s.Qexud = s.exuds_in * (2 * np.pi * r_in * s.segLength)
print('set default', l)
s.wilting_point = -15000
s.setCriticalPressure(s.wilting_point)  # for boundary conditions

cell_volumes = comm.bcast(s.getCellVolumes() , root = 0) #cm3
s.numFluidComp = 2

buTotCBeforeEach = comm.bcast(s.getTotCContent_each(), root = 0) 

print('buTotCBeforeEach',buTotCBeforeEach[[0,-1]])
#outer_R_bc_sol =  sources_sol*0.          

s.ddt = 1.e-5

#s.base.reportParams()
#print(s.Qexud/ 24. / 3600.)

for i, dt in enumerate(np.diff(times)):
    s.ddt =min( 1.e-5,s.ddt)

    if rank == 0:
        print("*****", "external time step", dt, 
              " d, simulation time", s.simTime, 
              "d, internal time step", s.ddt, "d")

    s.solve(dt, maxDt = 250./(24.*3600.))

    outer_R_bc = -np.transpose(s.getFlux_10c())# mol
    bulkSoil_sources = np.transpose(s.getSource_10c()) 


    outer_R_bc_sol = outer_R_bc[1:]# mol
    sources_sol =  bulkSoil_sources[1:] #mol
    #  mol C / cm3 scv
    
    print('outer_R_bc_sol',#outer_R_bc_sol,
            outer_R_bc_sol[0])
    buTotCAfterEach = comm.bcast(s.getTotCContent_each(), root = 0) 
    print('buTotCAfterEach',buTotCAfterEach[[0,-1]])
    print('buTotCBefore',sum(buTotCBeforeEach.flatten()),
            'buTotCAfter',sum(buTotCAfterEach.flatten()),
            's.Qexud*dt',s.Qexud*dt,
            'diff',sum(buTotCAfterEach.flatten()) - sum(buTotCBeforeEach.flatten()) - s.Qexud*dt)
    
    s.bulkMassCError1dsAll_real = buTotCAfterEach - buTotCBeforeEach 
    print('shapes', s.bulkMassCError1dsAll_real.shape,
        sources_sol.shape,outer_R_bc_sol.shape)
    s.bulkMassCError1dsAll_real[:-1] -= sources_sol + outer_R_bc_sol *0
    print('s.bulkMassCError1dsAll_real_A',#s.bulkMassCError1dsAll_real, 
            'cs',s.bulkMassCError1dsAll_real[0],
           'css1', s.bulkMassCError1dsAll_real[-1], 
           'cstot',s.bulkMassCError1dsAll_real[0] + s.bulkMassCError1dsAll_real[-1])#.flatten())
    s.bulkMassCError1dsAll_real[0] += s.bulkMassCError1dsAll_real[-1]#put together CS and CSS1
    s.bulkMassCError1dsAll_real[-1][:] = 0.
    print('s.bulkMassCError1dsAll_real_B',s.bulkMassCError1dsAll_real[0],
            sum(s.bulkMassCError1dsAll_real.flatten()))
    #print(sum(np.array(s.base.computeCSS1s())*cell_volumes/1e6)*s.f_sorp)
    buTotCBeforeEach = buTotCAfterEach
    
    
print('cell_volumes',cell_volumes)
print('points',s.getPoints())
flfl =s.getFlux_10c_()
flfl_cs = np.array([ff[1] for ff in flfl])
print('getFlux_10c_',flfl_cs)
print('getFace2CellIds_', s.getFace2CellIds_())


