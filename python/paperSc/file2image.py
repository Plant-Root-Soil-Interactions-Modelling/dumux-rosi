""" 
    Maize using rhizosphere models  
"""
import matplotlib; matplotlib.use('agg')
import sys;
import os


sys.path.append("../modules_fpit/");
sys.path.append("../../../CPlantBox/");
sys.path.append("../../../CPlantBox/src")

from importlib import reload
import plantbox as pb  # CPlantBox
import visualisation.vtk_plot as vp
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()
import timeit
import numpy as np

import scenario_setup
scenario_setup = reload(scenario_setup)
import rhizo_modelsPlant  # Helper class for cylindrical rhizosphere models
rhizo_modelsPlant = reload(rhizo_modelsPlant)
from rhizo_modelsPlant import *
#import evapotranspiration as evap
#import cyl_exu
import cyl3plant as cyl3
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank(); max_rank = comm.Get_size()
import os
from scenario_setup import write_file_array, write_file_float, div0, div0f

#DumuxDune27/DUMUX/dumux-rosi/python/paperSc/results/newMucil4p/baseline_1476_17_10to25_20mn_0s_128/C3_17b.vtp
#results_dir="./results/newMucil4p/baseline_1476_17_10to25_20mn_0s_128/"
results_dir="./results/testspeed/thr4baseline_1476_5_10to25_20mn_0s_32/vtpvti/"
vp.plot_roots_and_soil_files(filename = "C3_25b", 
pname = "[C3] (mol/cm3)",path = results_dir, interactiveImage = False)  # VTK vizualisation
