import sys; 
#from joblib import Parallel, delayed
import math
import os
import numpy as np
import time
import psutil
start_time_ = time.time()

from uqrMaster import launchUQR

    
    
simDuration = sys.argv[1]
condition = sys.argv[2]
directoryN = "/"+os.path.basename(__file__)[:-3]+"/"+simDuration+condition+"/"

simDuration = float(simDuration)

main_dir=os.environ['PWD']#dir of the file
results_dir = main_dir +"/results"+directoryN
if not os.path.exists(results_dir):
    os.makedirs(results_dir)
else:
    import shutil
    shutil.rmtree(results_dir)
    os.makedirs(results_dir)


launchUQR(directoryN,simDuration, condition)
