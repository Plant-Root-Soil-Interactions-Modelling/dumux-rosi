""" 
    Maize using rhizosphere models  
"""
import matplotlib; matplotlib.use('agg')
import sys;
import os
absolutePath = False
if absolutePath:
    os.chdir( "/home/rbtlm640/DUMUXr/dumux-rosi/python/fpit")
    sys.path.append( "/home/rbtlm640/DUMUXr/dumux-rosi/python/fpit")
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

import multiprocessing
import time
import scenario_setup as scenario
scenario = reload(scenario)
import rhizo_modelsPlant  # Helper class for cylindrical rhizosphere models
rhizo_modelsPlant = reload(rhizo_modelsPlant)
from rhizo_modelsPlant import *
import evapotranspiration as evap
#import cyl_exu
import cyl3plant as cyl3
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank(); max_rank = comm.Get_size()
import os
import pandas as pd
from scenario_setup import write_file_array, write_file_float, div0, div0f


if __name__ == '__main__':

    
    dt = 1/2/24
    k_iter = 20
    simMax = 3.
    lightType =""#+- "nolight" # or ""
    extraName = "testerror3d"
    results_dir="./results/"+extraName+str(k_iter)+"k_"\
                    +"_"+str(int(dt*24*60))+"mn_"\
                    +str(int((dt*24*60 - int(dt*24*60))*60))+"s_"\
                    +str(max_rank)+"/"
    print('results_dir',results_dir)
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
    comm.barrier()
    """
     Cylindric rhizosphere models, C exudate from finite volumes
    """

    usemoles = True
    """ parameters   """
    soil_ = scenario.vg_SPP(0)

    min_b = [-5, -5, -10.] 
    max_b = [5, 5, 0.] 
    cell_number = [5,5,20]
    #min_b = [-5., -5, -5.] 
    #max_b = [5., 5, 0.] 
    #cell_number = [5, 5, 5]
    sim_time = 1 #154   #  [day]
    # dt_inner = 360 / (24 * 3600)  # time step [day] 20


    """ rhizosphere model parameters """
    wilting_point = -15000  # cm
    recreateComsol = False

    periodic = False
    nc = 10

    logbase = 0.5  # according to Mai et al. (2019)
    mode = "dumux_w"  

    """ initialize """
    p_mean = pd.read_csv("./error3d/pHead.txt",delimiter=",", header=None).iloc[0].to_numpy() 
    
    soil_type = "loam"

    s, soil = scenario.create_soil_model(soil_type, 2013, soil_,#comp, 
                min_b, max_b, cell_number, demoType = mode, times = None, net_inf = None,
                usemoles = usemoles, dirResults = results_dir, p_mean_ = p_mean)
    currentDay = 0.
    
    
    with open("./error3d/source0bis.txt") as f:
        j = f.readlines()[0]
        res = eval(j)

    # run a loop
    for i in res:
        res[i] = -abs(10000*res[i])
    print('source0',res)
    s.setSource(res.copy(), eq_idx = 0)  # [mol/day], in modules/richards.py

    # bar
    def bar(soil, dt):
        s.solve(dt, maxDt = 250/(3600*24))
    while currentDay < simMax: #for i, dt in enumerate(np.diff(times)):

        currentDay += dt
        print("Day", currentDay)
        comm.barrier()
        s.setParameter("Newton.MaxRelativeShift", str(1e-10))
        redoSolve = True
        while redoSolve:
            s.ddt = 1.e-5
            try:
                comm.barrier()
                print("entering the s.solve", rank)
                    
                s.solve(dt, maxDt = 20/(3600*24))  # in modules/solverbase.py
                print("leaving the s.solve", rank)
                comm.barrier()
                redoSolve = False
                # newton parameters are re-read at each 'solve()' calls
                s.setParameter("Newton.MaxRelativeShift", str(s.MaxRelativeShift))# reset value
            except Exception as err:
                print(rank, f"Unexpected {err=}, {type(err)=}")
                if k_soil_solve < 5:
                    print(rank,
                          'soil.solve() failed. Increase NewtonMaxRelativeShift from',
                          maxRelShift,'to',maxRelShift*10.)
                    maxRelShift *= 10.
                    # newton parameters are re-read at each 'solve()' calls
                    s.setParameter("Newton.MaxRelativeShift", str(maxRelShift))
                    s.reset()
                    k_soil_solve += 1
                else:
                    raise Exception
        water_content =comm.bcast( np.array(s.getWaterContent()), root = 0) 
        SolutionHead =comm.bcast( np.array(s.getSolutionHead()), root = 0) 
        write_file_array('getWaterContent',water_content, directory_ =results_dir, fileType = '.csv')
        write_file_array('getSolutionHead',SolutionHead, directory_ =results_dir, fileType = '.csv')
        

"""

                p = multiprocessing.Process(target=bar,args=(s,dt))
                p.start()
                p.join(10*60*60)
                
                # If thread is still active
                if p.is_alive():
                    print("running... let's kill it...")

                    # Terminate - may not work if process is stuck for good
                    # p.terminate()
                    # OR Kill - will work for sure, no chance for process to finish nicely however
                    p.kill()
                    raise Exception

"""