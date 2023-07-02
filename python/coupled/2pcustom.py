import sys; 
#from joblib import Parallel, delayed
import math
import os
import numpy as np
import time
import psutil
start_time_ = time.time()

from uqrMaster_2p import launchUQR

    
    
simInit = int(sys.argv[1])#14
simEnd = int(sys.argv[2])#34
weather = sys.argv[3]#"wet"
phenotype = sys.argv[4]#"deep"
extraname = sys.argv[5]#"deep"
    
    
directoryN = "/"+os.path.basename(__file__)[:-3]+sys.argv[1]+weather+phenotype+sys.argv[2]+extraname+"/"

main_dir=os.environ['PWD']#dir of the file
results_dir = main_dir +"/results"+directoryN
if not os.path.exists(results_dir):
    os.makedirs(results_dir)
else:
    import shutil
    shutil.rmtree(results_dir)
    os.makedirs(results_dir)


launchUQR(directoryN,simInit,simEnd, weather,phenotype)
