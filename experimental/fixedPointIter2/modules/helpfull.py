
import numpy as np
import pandas as pd
import timeit
import pandas as pd
import matplotlib; matplotlib.use('agg')
import sys;
import os
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank(); max_rank = comm.Get_size()
import timeit
import ctypes
import numbers
import tempfile
import signal

import threading
import multiprocessing


class TimeoutException(Exception):
    pass

def timeout_handler(signum, frame):
    raise TimeoutException

# Set the timeout handler
signal.signal(signal.SIGALRM, timeout_handler)


def check_os():
    if os.name == 'nt':
        return "Windows"
    elif os.name == 'posix':
        return "Linux/Unix"
    else:
        return "Unknown OS"


def run_with_timeout__(timeout, func, *args, **kwargs):
    class FuncThread(threading.Thread):
        def __init__(self):
            super().__init__()
            self.result = None
            self.exception = None

        def run(self):
            try:
                self.result = func(*args, **kwargs)
            except Exception as e:
                self.exception = e

    thread = FuncThread()
    thread.start()
    thread.join(timeout)

    if thread.is_alive():
        raise TimeoutException("Function call timed out")
    elif thread.exception:
        raise thread.exception
    else:
        return thread.result

def run_with_timeout(timeout, func, *args, **kwargs):
    # Start the timer
    signal.alarm(int(timeout))
    try:
        result = func(*args, **kwargs)
    finally:
        # Disable the alarm
        signal.alarm(0)
    return result

def is_number(obj):
    # Check for standard Python numeric types (int, float) and NumPy numeric types
    return isinstance(obj, (numbers.Number, np.number))


class StdoutRedirector:
    def __init__(self):
        self.libc = ctypes.CDLL(None)
        self.c_stdout = ctypes.c_void_p.in_dll(self.libc, "stdout")
        self.buffer = None
        self.filepath = None

    def __enter__(self):
        self.old_stdout_fd = os.dup(1)
        self.temp_file = tempfile.NamedTemporaryFile(delete=False)
        self.temp_file_fd = self.temp_file.fileno()
        os.dup2(self.temp_file_fd, 1)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        os.dup2(self.old_stdout_fd, 1)
        os.close(self.old_stdout_fd)
        self.temp_file.close()

        with open(self.temp_file.name, 'r') as temp_file:
            self.buffer = temp_file.read()

        os.remove(self.temp_file.name)

        if exc_type is not None:
            with open(self.filepath, 'w') as f:
                f.write(self.buffer)

        return False  # Do not suppress exceptions
        
def suggestNumStepsChange(dt, dt_inner, failedLoop, perirhizalModel, n_iter_inner_max):# taken from dumux
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
    dt : outer time step
    failedLoop: inner fixed pint iteration failed
    perirhizalModel: model object, stores usefull data
    n_iter_inner_max: max number of necessary iteration during the last fixed point iteration
    """
    targetIter_ = perirhizalModel.targetIter # objective max number of iteration for inner loop
    #dt_inner = perirhizalModel.dt_inner # time step of inner loop
    results_dir = perirhizalModel.results_dir
    nOld = int(dt/dt_inner) # old number of step for inner loop
    
    
    ## adapt inner loop time step:
    if perirhizalModel.doNestedFixedPointIter:
        NinnerOld = int(np.ceil(dt_inner/perirhizalModel.dt_inner2 ))  
    
    if failedLoop:# if the inner computation failed, we also need to decrease the timestep
        numIter_ = perirhizalModel.k_iter # k_iter == max accepted number of iteration for innerl loop
    else:
        numIter_ = n_iter_inner_max
    if (numIter_ > targetIter_) :
        percent = float(numIter_ - targetIter_)/float(targetIter_)
        change = 1.0 + percent
    else:
        percent = float(targetIter_ - numIter_)/float(targetIter_)
        change = 1 /(1.0 + percent/1.2)
    nNew = int(np.ceil(nOld * change))
    
    try:
        assert nOld == int(nOld)
        assert nNew == int(nNew)
    except:
        print('nNew iisue',nNew , int(nNew), nOld,dt,dt_inner,
              (nOld == int(nOld)), (nNew == int(nNew)))

    try:
        dt_inner = dt/float(nNew)
    except:
        write_file_array("suggestNumStepsChangeError", np.array([dt, dt_inner, failedLoop,n_iter_inner_max,
                                                                targetIter_,percent,change,nNew ,  nOld]), 
                         directory_ =results_dir, fileType = '.csv') 
        raise Exception
        
    write_file_array("suggestNumStepsChange", np.array([numIter_,n_iter_inner_max,targetIter_,percent,change,
                                                    nNew ,  nOld,
                                                 dt, dt_inner , failedLoop]), 
                         directory_ =results_dir, fileType = '.csv') 
    
    ## adapt inner loop time step:
    if perirhizalModel.doNestedFixedPointIter:  
        perirhizalModel.dt_inner2 = dt_inner/NinnerOld
    else:
        perirhizalModel.dt_inner2 = dt_inner
        
    
    return dt_inner


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
        for extraText in ["","cyl_val/","printData/", "vtpvti/", "fpit/","fpit2/", "fpit2/cyl_val/"]:
            if not os.path.exists(results_dir+extraText):
                os.makedirs(results_dir+extraText)
            else:
                test = os.listdir(results_dir+extraText)
                for item in test:
                    try:
                        os.remove(results_dir+extraText+item)
                    except:
                        pass


def sinusoidal3(t, dt):
    """ sinusoidal functgion from 5:00 - 23:00, 0 otherwise (max is 1)"""
    return ((np.maximum(0.,  0.75+((np.cos(2 * np.pi * (t - 0.5)) + np.cos(2 * np.pi * ((t + dt) - 0.5)) )/ 2))))/1.75

 

def getPsiAir(RH, TairC):#constants are within photosynthesys.h
    Mh2o = 18 #molar weight, g mol-1, mg mmol-1
    R_ph = 83.143 #perfect gas constant, hPa cm3K−1mmol−1
    rho_h2o = 1000 #water density mg cm-3
    return np.log(RH) * rho_h2o * R_ph * (TairC + 237.3)/Mh2o * (1/0.9806806)  ; #in cm
 
def continueLoop(perirhizalModel,n_iter, dt_inner: float,failedLoop: bool,
                         real_dt : float,
                 name="continueLoop", isInner = False,doPrint =
                 True,
                         fileType = '.csv' , plant = None, FPIT_id = 2):
                         
    gaveUp = (n_iter >= perirhizalModel.k_iter)
    cL = None
    if rank ==0:
        results_dir = perirhizalModel.results_dir
        plant.time_rhizo_cumul += plant.time_rhizo_i
        plant.time_3ds_cumul += plant.time_3ds_i
        plant.time_rhizo_i = 0
        plant.time_3ds_i = 0
        plant.time_plant_cumul = plant.time_plant_cumulW + plant.time_plant_cumulS 

        write_file_array("totalComputetime",
                         np.array([timeit.default_timer() - plant.time_start_global,
                        plant.time_plant_cumul,plant.time_rhizo_cumul ,plant.time_3ds_cumul]) , 
                         directory_ = results_dir)

        # stop if converged or gave up
        if FPIT_id ==2:
            cL = ((np.floor(perirhizalModel.err) > perirhizalModel.max_err) or
                   perirhizalModel.solve_gave_up or failedLoop
                    or (np.floor(perirhizalModel.diff1d3dCurrant_rel*1000.)/1000.>0.01)
                    or (np.floor(perirhizalModel.maxdiff1d3dCurrant_rel*1000.)/1000.>0.01)
                    or (min(perirhizalModel.new_soil_solute.reshape(-1)) < 0)
                    or ((n_iter < perirhizalModel.minIter) and (isInner)))  and (n_iter < perirhizalModel.k_iter)
        else:
            cL = ((np.floor(perirhizalModel.err2) > perirhizalModel.max_err) or
                   perirhizalModel.solve_gave_up or failedLoop
                    or ((n_iter < perirhizalModel.minIter)
                        and (isInner)))  and (n_iter < perirhizalModel.k_iter)


        if (rank == 0) and isInner:
            if FPIT_id == 2:
                print(f'continue loop? {bool(cL)}\n\tn_iter: {n_iter}, non-convergence and error metrics: {perirhizalModel.err:.2e}\n\ttotal relative 1d-3d difference added at this time step: {perirhizalModel.diff1d3dCurrant_rel:.2e}\n\tmax relative 1d-3d difference added at this time step: {perirhizalModel.maxdiff1d3dCurrant_rel:.2e}')
            if FPIT_id == 3:
                print(f'\t\t\t\tcontinue INNER loop? {bool(cL)}\n\t\t\t\tn_iter: {n_iter}, '
                      f'non-convergence and error metrics: {perirhizalModel.err2:.2e}, '
                      f'solver gave up for 1ds: {bool(perirhizalModel.solve_gave_up)}')
                print('\t\t\t\tD(psi_rsi)',f'{perirhizalModel.errWrsi:.2e},',
                      "psi_{rsi,in} vs psi_{rsi,out}:",f'{perirhizalModel.errWrsiRealInput:.2e}')

        if doPrint:
            if not os.path.isfile(results_dir+name+fileType):
                write_file_array(name, np.array(['n_iter', 'err', 
                                                 'diff1d3dCurrant_rel','maxdiff1d3dCurrant_rel',
                                                 'solve_gave_up', 'min__soil_solute',
                                                         'dt_inner','real_dt ',
                                                 'failedLoop','cL']), directory_ =results_dir, fileType = fileType)
                if not perirhizalModel.doMinimumPrint:
                    write_file_array(name+"Bool", np.array(['n_iter',  'err', 

                                                            'diff1d3dCurrant_rel','maxdiff1d3dCurrant_rel',
                                                 'solve_gave_up',  'min__soil_solute',
                                                             'dt_inner','failedLoop','cL']), directory_ =results_dir, fileType = fileType)

            write_file_array(name, np.array([n_iter, perirhizalModel.err, 
                                             perirhizalModel.diff1d3dCurrant_rel,perirhizalModel.maxdiff1d3dCurrant_rel,
                                             perirhizalModel.solve_gave_up, min(perirhizalModel.new_soil_solute.reshape(-1)),
                                                         dt_inner, real_dt ,failedLoop,
                                             cL]), directory_ =results_dir, fileType = fileType)
            if not perirhizalModel.doMinimumPrint:
                write_file_array(name+"2", 
                                 np.concatenate((perirhizalModel.sumDiff1d3dCW_rel,perirhizalModel.maxDiff1d3dCW_rel)),  
                                 directory_ =results_dir, fileType = fileType)
                write_file_array(name+"Bool", np.array([n_iter, (np.floor(perirhizalModel.err) > perirhizalModel.max_err), 
                                                        (np.floor(perirhizalModel.diff1d3dCurrant_rel*10.)/10. > 0.1),
                                                        (np.floor(perirhizalModel.maxdiff1d3dCurrant_rel) >1),
                                                        #(abs(perirhizalModel.rhizoMassWError_abs) > 1e-13), (abs(perirhizalModel.rhizoMassCError_abs) > 1e-9), 
                                                        #(max(abs(perirhizalModel.errDiffBCs*0)) > 1e-5), 
                                                        perirhizalModel.solve_gave_up, min(perirhizalModel.new_soil_solute.reshape(-1)) < 0.,
                                                             dt_inner,failedLoop,cL]), directory_ =results_dir, fileType = fileType)

    cL = comm.bcast(cL,root = 0)
        
    if perirhizalModel.debugMode:
        failedLoop_ = np.array( comm.bcast(comm.gather(failedLoop,root = 0),root = 0))
        assert (failedLoop_ ==failedLoop_[0]).all() # all true or all false
    
    return cL, gaveUp

def checkseg2cellMapping(seg2cell_old, plantModel):
    # make sure that, once a segment is created, it stays in the same soil voxel
    for segid in seg2cell_old.keys():
        assert seg2cell_old[segid] == plantModel.seg2cell[segid]

    # also segs are in only one voxel
    cell2segVals = np.concatenate((list(plantModel.cell2seg.values()))).flatten()
    if len(cell2segVals) != len(set(cell2segVals)):#make sure all values are unique
        print(plantModel.seg2cell)
        print(plantModel.cell2seg)
        print(cell2segVals)
        print(len(cell2segVals), len(set(cell2segVals)))
        raise Exception

def resetAndSaveData(perirhizalModel):
    perirhizalModel.rhizoMassWError_abs = 1.
    perirhizalModel.rhizoMassCError_abs = 1.
    perirhizalModel.errDiffBCs = 1.
    perirhizalModel.err = 1.
    perirhizalModel.max_err = 1.
    perirhizalModel.diff1d3dCurrant_rel = 1e6
    perirhizalModel.maxdiff1d3dCurrant_rel = 1e6

    n_iter = 0 # number of iteration in the loop
    keepGoing = True # stay in the fixed point iteration
    
    
    return n_iter, keepGoing


def resetPlantWFlux(plantModel,perirhizalModel):
    perirhizalModel.err2 = 1.
    plantModel.seg_fluxes0Cumul_inner = 0
    plantModel.seg_fluxes1Cumul_inner = 0
    plantModel.seg_fluxes2Cumul_inner = 0
    plantModel.TranspirationCumul_inner = 0 # reset transpiration of inner time step to 0
    plantModel.AnCumul_inner = 0 # reset transpiration of inner time step to 0


def saveData(plantModel, perirhizalModel, s):
    resetPlantWFlux(plantModel,perirhizalModel)

    # save data before entering iteration loop
    s.saveManual()
    perirhizalModel.saveManual()
    perirhizalModel.leftSpellBU = perirhizalModel.leftSpell
    perirhizalModel.enteredSpellBU = perirhizalModel.enteredSpell

def resetData(plantModel, perirhizalModel, s):

    s.resetManual()
    perirhizalModel.resetManual()
    perirhizalModel.leftSpell = perirhizalModel.leftSpellBU
    perirhizalModel.enteredSpell = perirhizalModel.enteredSpellBU

    # to get the error values for the old solution vector
    perirhizalModel.check1d3dDiff()
    


def getCumulativeTranspirationAg(plantModel, perirhizalModel, dt):
    if rank == 0:
        deltalen = len(plantModel.seg_fluxes0Cumul_inner)-len(plantModel.seg_fluxes0Cumul)# plant grew?
        if deltalen > 0:
            plantModel.seg_fluxes0Cumul = np.concatenate((plantModel.seg_fluxes0Cumul, np.zeros(deltalen))) 
            plantModel.seg_fluxes1Cumul = np.concatenate((plantModel.seg_fluxes1Cumul, np.zeros(deltalen))) 
            plantModel.seg_fluxes2Cumul = np.concatenate((plantModel.seg_fluxes2Cumul, np.zeros(deltalen))) 

        plantModel.seg_fluxes0Cumul += plantModel.seg_fluxes0Cumul_inner 
        plantModel.seg_fluxes1Cumul += plantModel.seg_fluxes1Cumul_inner
        plantModel.seg_fluxes2Cumul += plantModel.seg_fluxes2Cumul_inner

        plantModel.seg_fluxes0 = plantModel.seg_fluxes0Cumul_inner/dt
        plantModel.seg_fluxes1 = plantModel.seg_fluxes1Cumul_inner/dt
        plantModel.seg_fluxes2 = plantModel.seg_fluxes2Cumul_inner/dt
        
        plantModel.TranspirationCumul += plantModel.TranspirationCumul_inner 
        if perirhizalModel.doPhotosynthesis:
            if perirhizalModel.enteredSpell and (not perirhizalModel.leftSpell):
                plantModel.TranspirationCumul_eval += plantModel.TranspirationCumul_inner
            elif perirhizalModel.leftSpell:
                plantModel.TranspirationCumul_eval = 0.
            
            plantModel.An = plantModel.AnCumul_inner/(dt*24*3600)
        
        