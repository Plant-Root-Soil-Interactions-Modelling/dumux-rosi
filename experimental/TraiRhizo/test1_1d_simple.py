import sys; sys.path.append("./modules_fpit");  
sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src")
# dumux38Exud/dumux/dumux-rosi/experimental/TraiRhizo
import matplotlib; matplotlib.use('agg')

#from rosi_richards10c_cyl import RichardsNCCylFoam  # C++ part (Dumux binding)
#from rosi_richards_cyl import RichardsCylFoam as RichardsNCCylFoam  # C++ part (Dumux binding)
from rosi_richardsnc_cyl import RichardsNCCylFoam  # C++ part (Dumux binding)

usemoles = False

from richards import RichardsWrapper  # Python part
from richards_no_mpi import RichardsNoMPIWrapper  # Python part of cylindrcial 

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()
import scenario_setup as stf #
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

#s = RichardsWrapper(Richards10CCylFoam(), usemoles)
year = 2013
soil_type = "loam"
genotype = "WT"
comp = "phenolics"
soil_ = scenario_setup.vg_SPP(0)
mode = "dumux_10c"
min_b = [-4., -4., -15.]  # cm
max_b = [4., 4., 0.]  # cm
p_mean = -659.8 + (max_b[2] - min_b[2]) / 2 
css1Function_ = 9
paramIndx = 5
dt = 30 / (24 * 3600) 
times = np.array([0.,dt])#,dt*2])  # days
#s3d, soil = scenario_setup.create_soil_model(soil_type, year, soil_,#comp, 
#                min_b, max_b, cell_number, 
demoType = mode
dirResults = results_dir
p_mean_ = p_mean
noAds = True

s=RichardsNoMPIWrapper(RichardsNCCylFoam(), usemoles)

s.dirResults = dirResults

s.css1Function = css1Function_

def getCSS2Init(CSW):
    return 0., 0.
s.getCSS2Init = getCSS2Init
s.noAds=noAds
r_in = 0.02
r_out = 0.2

s.initialize()
length_ = 2.
stf.setShape1D(s,r_in, r_out,length = length_,nCells = 10,doLogarithmic=True)

stf.setDefault(s)

s.setParameter("Component.MolarMass", "3.1e-2")  #  kg/mol
s.setParameter("Newton.MaxSteps", "200")
s.setParameter("Newton.MaxTimeStepDivisions", "100")
s.setParameter("Problem.EnableGravity", "False")
#s.setParameter("Newton.EnableResidualCriterion", "true") # sometimes helps, sometimes makes things worse
#s.setParameter("Newton.EnableAbsoluteResidualCriterion", "true")
#s.setParameter("Newton.SatisfyResidualAndShiftCriterion", "true")
s.setParameter("Problem.reactionExclusive", "0")  
s.setParameter("Newton.MaxRelativeShift","1e-10")
s.setParameter("Newton.Verbosity","0")
s.setParameter("Newton.EnableChop", "True")

s.setParameter( "Soil.css1Function", str(s.css1Function))

stf.setSoilParam(s,paramIndx)
stf.getBiochemParam(s,paramIndx,s.noAds)
stf.setBiochemParam(s)
stf.setIC(s,paramIndx)
s.win = -0.1; s.exudl_in = 0.
s.exuds_in = 0

stf.setBC1D(s)
s, __ = stf.setupOther(s, p_mean_)
#s.initializeProblem()
print('get length')
l = length_#s.segLength
# mol/d or g/d
s.Win = s.win * (2 * np.pi * r_in * length_)#s.segLength)
print('set default', l)
s.wilting_point = -15000
s.setCriticalPressure(s.wilting_point)  # for boundary conditions

#cell_volumes = comm.bcast(s.getCellVolumes() , root = 0) #cm3
#print('cell_volumesA',cell_volumes)
s.numFluidComp = 2

#outer_R_bc_sol =  sources_sol*0.          

s.ddt = 1.e-5

cell_volumes = s.getCellSurfacesCyl() * length_
#print('cell_volumesB',cell_volumes)
#raise Exception

#[0.00167912 0.00280094 0.00467225 0.00779378 0.01300081 0.02168667
# 0.03617554 0.06034444 0.10066059]
#[0.00167912 0.00280094 0.00467225 0.00779378 0.01300081 0.02168667
# 0.03617554 0.06034444 0.10066059]
#s.base.reportParams()


for i, dt in enumerate(np.diff(times)):
    s.ddt =min( 1.e-5,s.ddt)

    if rank == 0:
        print("*****", "external time step", dt, 
              " d, simulation time", s.simTime, 
              "d, internal time step", s.ddt, "d", s.css1Function)


    buWBefore_ =  s.getWaterContent()*cell_volumes
    buWBefore = sum(buWBefore_ )
        
    s.solve(dt)

    buWAfter_ =  s.getWaterContent()*cell_volumes
    buWAfter = sum(buWAfter_ )
    if True:
        Q_in_m, Q_out_m = s.getSavedBC(s.base.rIn*100, 
                    s.base.rOut*100, length_)   
        QflowIn_limited = Q_in_m[0]
        QflowOut_limited = Q_out_m[0]
                                    
    
        diffW = buWAfter - ( buWBefore + (QflowIn_limited + QflowOut_limited) * dt)
        diffW2 = buWAfter - ( buWBefore + (s.Win + QflowOut_limited) * dt)
        print('diffW',diffW,diffW2, buWAfter , buWBefore )
        print('QflowIn_limited',QflowIn_limited,s.Win , dt)
    else:
        print(buWAfter , buWBefore,buWAfter - buWBefore )
    buWBefore = buWAfter
    
if False:
    print('cell_volumes',cell_volumes)
    print('points',s.getPoints())
    flfl =s.getFlux_10c_()
    flfl_cs = np.array([ff[1] for ff in flfl])
    print('getFlux_10c_',flfl_cs)
    print('getFace2CellIds_', s.getFace2CellIds_())


