
import numpy as np
import pandas as pd
import timeit
import pandas as pd
import matplotlib; matplotlib.use('agg')
import sys;
import os
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank(); max_rank = comm.Get_size()
import timeit

def suggestNumStepsChange(nOld, numIter_, targetIter_, results_dir):# taken from dumux
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
    if (numIter_ > targetIter_) :
        percent = float(numIter_ - targetIter_)/float(targetIter_)
        change = 1.0 + percent
    else:
        percent = float(targetIter_ - numIter_)/float(targetIter_)
        change = 1 /(1.0 + percent/1.2)
    return int(np.ceil(nOld * change))



def div0(a, b, c):   # function to avoid division by 0    
    if isinstance(a,(list, np.ndarray)):  
        return np.divide(a, b, out=np.full(len(a), c), where=b!=0)
    else:
        return div0f(a, b, c)
    
def div0f(a, b, c):    # function to avoid division by 0     
    if b != 0:
        return a/b
    else:
        return a/c
    
    

def write_file_float(name, data, directory_, allranks = False):
    if (rank == 0) or allranks:
        name2 = directory_+ name+ '.txt'
        #print('write_file_float',name2, allranks)
        modeOp = 'a'
        if  False:#'fpit' in name:
            modeOp = 'w'
        with open(name2, modeOp) as log:
            log.write(repr( data)  +'\n')
        
def write_file_array(name, data, space =",", directory_ ="./results/", fileType = '.txt',
                     allranks = False ):
    if (rank == 0) or allranks:
        try:

            np.array(data).reshape(-1)
            modeOp = 'a'
            if False:#'fpit' in name:
                modeOp = 'w'
            name2 = directory_+ name+ fileType
            #print('write_file_array',name2)
            with open(name2, modeOp) as log:
                log.write(space.join([num for num in map(str, data)])  +'\n')
        except:
            print('write_file_array',name, data,data.shape)
            raise Exception

def setupDir(results_dir):
    """if results directory already exists, make it empty"""
    if rank == 0:
        for extraText in ["","cyl_val/","printData/"]:
            if not os.path.exists(results_dir+extraText):
                os.makedirs(results_dir+extraText)
            else:
                test = os.listdir(results_dir+extraText)
                for item in test:
                    try:
                        os.remove(results_dir+extraText+item)
                    except:
                        pass