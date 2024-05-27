""" 
    Maize using rhizosphere models  
"""
import matplotlib; matplotlib.use('agg')
import sys;
import os
absolutePath = False
if absolutePath:
    os.chdir( "/home/rbtlm640/DUMUXr/dumux-rosi/python/paperSc")
    sys.path.append( "/home/rbtlm640/DUMUXr/dumux-rosi/python/paperSc")
    sys.path.append("/home/rbtlm640/DUMUXr/dumux-rosi/python/fpit/data/");
    sys.path.append("/home/rbtlm640/DUMUXr/dumux-rosi/python/modules_fpit/");
    sys.path.append("/home/rbtlm640/DUMUXr/CPlantBox/");#absolut path work better when
    # using pycharm
    sys.path.append("/home/rbtlm640/DUMUXr/CPlantBox/src")
else:
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


rangecoa = [5.4e-8,9.2e-6]

results_dir="./results/vtivtp/"
vp.plot_roots_and_soil_files(filename = "C3_17e", 
pname = "[C3] (mol/cm3)",path = results_dir, myrange = rangecoa)  # VTK vizualisation
