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
    # using pycharmSoil_solute_conc
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



"""
     * \brief Suggest a new number of time steps
     *
     * The default behavior is to suggest the old time-step number
     * scaled by the ratio between the target iterations and the
     * iterations required to actually solve the last time-step.
        // be aggressive increasing the time-step number but
        // conservative when decreasing it. the rationale is
        // that we want to avoid failing in the next iteration
        // nNew > nOld ==> dtnew < dtOld
"""
def suggestNumStepsChange(nOld, numIter_, targetIter_, results_dir):# taken from dumux
    if (numIter_ > targetIter_) :
        percent = float(numIter_ - targetIter_)/float(targetIter_)
        change = 1.0 + percent
    else:
        percent = float(targetIter_ - numIter_)/float(targetIter_)
        change = 1 /(1.0 + percent/1.2)
    #if not doMinimumPrint:
    #    write_file_array("suggestNumStepsChange",np.array([nOld, numIter_, targetIter_, percent,change, np.ceil(nOld * change)]), directory_ =results_dir, fileType = '.csv') 
    return int(np.ceil(nOld * change))
    

if __name__ == '__main__':

    #initsim =float(sys.argv[1])# initsim = 9.5
    #mode = sys.argv[2] #"dumux_w" "dumux_3c" "dumux_10c" 
    dt = 20/60/24
    #p_mean = -100
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
    css1Function_ = 9
    lightType =""#+- "nolight" # or ""
    mpiVerbose = True
    mpiVerboseInner = False
    noAds = False
    doSimple =False
    doMinimumPrint =  True
    if True:
        min_b = np.array([-3./2, -12./2, -40.])
        cell_number =np.array( [3,12,41]) # 1cm3 
        max_b =np.array( [3./2, 12./2, 0.])
    else:
        min_b = np.array([0,0,0])
        cell_number =np.array( [1,1,1]) # 1cm3 
        max_b =np.array( [1,1,1])
    
    results_dir="./results/soilshape/"+str(int(np.prod(cell_number)))+"bis/"
    
    #comm.barrier()
    if rank == 0:
        print('results_dir','DumuxDune27',results_dir, flush = True)
    #comm.barrier()
    if rank == 0:
        for extraText in [""]:
            if not os.path.exists(results_dir+extraText):
                os.makedirs(results_dir+extraText)
            else:
                test = os.listdir(results_dir+extraText)
                for item in test:
                    try:
                        os.remove(results_dir+extraText+item)
                    except:
                        pass
    #comm.barrier()
    """
     Cylindric rhizosphere models, C exudate from finite volumes
    """

    soil_ = scenario_setup.vg_SPP(0)
    

    """ initialize """

    s, soil = scenario_setup.create_soil_model("loam", 2013, soil_,#comp, 
                min_b, max_b, cell_number, demoType = None, times = None, net_inf = None,
                usemoles = True, dirResults = results_dir, p_mean_ = -100, 
                                         css1Function = 9,
                                        paramIndx=5,
                                        noAds = noAds)

    write_file_array("cellVol", np.array(s.getCellVolumes()), directory_ =results_dir) # cm3 
    write_file_array("cellIds", np.array(s.cellIndices), directory_ =results_dir) # cm3
    cellcenters = s.getCellCenters()
    cellcentersX = []
    cellcentersY = []
    cellcentersZ = []
    for sub_array in cellcenters:
        cellcentersX.append(sub_array[0])
        cellcentersY.append(sub_array[1])
        cellcentersZ.append(sub_array[2])

    write_file_array("cellcentersX", np.array(cellcentersX), directory_ =results_dir) # cm3
    write_file_array("cellcentersY", np.array(cellcentersY), directory_ =results_dir) # cm3
    write_file_array("cellcentersZ", np.array(cellcentersZ), directory_ =results_dir) # cm3
    allcconcent= [s.getContent(nc, nc<=2).mean() for nc in range(1,9)]
    print(allcconcent, sum(allcconcent) )
    
