
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

class TimeoutException(Exception):
    pass

def timeout_handler(signum, frame):
    raise TimeoutException

# Set the timeout handler
signal.signal(signal.SIGALRM, timeout_handler)


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

        if exc_type is not None:
            self.buffer = open(self.temp_file.name, 'r').read()
            with open(self.filepath, 'w') as f:
                f.write(self.buffer)
        os.remove(self.temp_file.name)
        
        
def suggestNumStepsChange(dt, failedLoop, perirhizalModel, n_iter_inner_max):# taken from dumux
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
    dt_inner = perirhizalModel.dt_inner # time step of inner loop
    results_dir = perirhizalModel.results_dir
    nOld = int(dt/dt_inner) # old number of step for inner loop
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


    dt_inner = dt/float(nNew)
    write_file_array("suggestNumStepsChange", np.array([numIter_,n_iter_inner_max,targetIter_,percent,change,
                                                    nNew ,  nOld,
                                                 dt, dt_inner , failedLoop]), 
                         directory_ =results_dir, fileType = '.csv') 
    
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
        for extraText in ["","cyl_val/","printData/", "vtpvti/", "fpit/", "fpit/cyl_val/"]:
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
 
def continueLoop(perirhizalModel,s,n_iter, dt_inner: float,failedLoop: bool,
                         real_dtinner: float,name="continueLoop", isInner = False,doPrint = True, 
                         fileType = '.csv' , plant = None):
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


    cL = ((np.floor(perirhizalModel.err) > perirhizalModel.max_err) or  perirhizalModel.solve_gave_up or failedLoop
            #or (s.bulkMassErrorWater_rel > s.bulkMassErrorWater_relLim*10)
            #or (perirhizalModel.rhizoMassWError_rel >
            #perirhizalModel.rhizoMassWError_relLim*10)
            #or (perirhizalModel.errWrsi > 1.)#already added in 'err' metric
            or (np.floor(perirhizalModel.diff1d3dCurrant_rel*1000.)/1000.>0.01)
            or (np.floor(perirhizalModel.maxdiff1d3dCurrant_rel*1000.)/1000.>0.01)
            or (min(perirhizalModel.new_soil_solute.reshape(-1)) < 0)  
            or ((n_iter < perirhizalModel.minIter) and (isInner)))  and (n_iter < perirhizalModel.k_iter)

    if (rank == 0) and isInner:
        print(f'continue loop? {bool(cL)}\n\t\tn_iter: {n_iter}, non-convergence and error metrics: {perirhizalModel.err:.2e}, solver gave up for 1ds: {bool(perirhizalModel.solve_gave_up)}\n\t\ttotal relative 1d-3d difference added at this time step: {perirhizalModel.diff1d3dCurrant_rel:.2e}\n\t\tmax relative 1d-3d difference added at this time step: {perirhizalModel.maxdiff1d3dCurrant_rel:.2e}')
        print('\t\tD(psi_rsi)',f'{perirhizalModel.errWrsi:.2e},',
              "psi_{rsi,in} vs psi_{rsi,out}:",f'{perirhizalModel.errWrsiRealInput:.2e}')

    cL = comm.bcast(cL,root = 0)
    failedLoop_ = np.array( comm.bcast(comm.gather(failedLoop,root = 0),root = 0))
    assert (failedLoop_ ==failedLoop_[0]).all() # all true or all false

    if doPrint:
        if not os.path.isfile(results_dir+name+fileType):
            write_file_array(name, np.array(['n_iter', 'err', 
                                             'diff1d3dCurrant_rel','maxdiff1d3dCurrant_rel',
                                             'solve_gave_up', 'min__soil_solute',
                                                     'dt_inner','real_dtinner','failedLoop','cL']), directory_ =results_dir, fileType = fileType)
            if not perirhizalModel.doMinimumPrint:
                write_file_array(name+"Bool", np.array(['n_iter',  'err', 

                                                        'diff1d3dCurrant_rel','maxdiff1d3dCurrant_rel',
                                             'solve_gave_up',  'min__soil_solute',
                                                         'dt_inner','failedLoop','cL']), directory_ =results_dir, fileType = fileType)

        write_file_array(name, np.array([n_iter, perirhizalModel.err, 
                                         perirhizalModel.diff1d3dCurrant_rel,perirhizalModel.maxdiff1d3dCurrant_rel,
                                         perirhizalModel.solve_gave_up, min(perirhizalModel.new_soil_solute.reshape(-1)),
                                                     dt_inner, real_dtinner,failedLoop,cL]), directory_ =results_dir, fileType = fileType)
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

    return cL

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

def resetAndSaveData1(perirhizalModel):
    perirhizalModel.rhizoMassWError_abs = 1.
    perirhizalModel.rhizoMassCError_abs = 1.
    perirhizalModel.errDiffBCs = 1.
    perirhizalModel.err = 1.
    perirhizalModel.max_err = 1.
    perirhizalModel.diff1d3dCurrant_rel = 1e6
    perirhizalModel.maxdiff1d3dCurrant_rel = 1e6

    n_iter = 0 # number of iteration in the loop
    failedLoop = False # had non-convergence error in dumux
    keepGoing = True # stay in the fixed point iteration
    
    return n_iter, failedLoop, keepGoing

def resetAndSaveData2(plantModel, perirhizalModel, s):
    plantModel.TranspirationCumul_inner = 0 # reset transpiration of inner time step to 0
    plantModel.AnCumul_inner = 0 # reset transpiration of inner time step to 0

    # save data before entering iteration loop
    s.saveManual()
    perirhizalModel.saveManual()
    perirhizalModel.leftSpellBU = perirhizalModel.leftSpell
    perirhizalModel.enteredSpellBU = perirhizalModel.enteredSpell

def resetAndSaveData3(plantModel, perirhizalModel, s):
    plantModel.TranspirationCumul_inner = 0 # reset transpiration of inner time step to 0
    plantModel.AnCumul_inner = 0 # reset transpiration of inner time step to 0
    s.resetManual()
    perirhizalModel.resetManual()
    perirhizalModel.leftSpell = perirhizalModel.leftSpellBU
    perirhizalModel.enteredSpell = perirhizalModel.enteredSpellBU

    # to get the error values for the old solution vector
    perirhizalModel.check1d3dDiff()
    
    

def getCumulativeTranspirationAg(plantModel, perirhizalModel, dt):
    plantModel.TranspirationCumul += plantModel.TranspirationCumul_inner 
    if perirhizalModel.doPhotosynthesis:
        if perirhizalModel.enteredSpell and (not perirhizalModel.leftSpell):
            plantModel.TranspirationCumul_eval += plantModel.TranspirationCumul_inner
        elif perirhizalModel.leftSpell:
            plantModel.TranspirationCumul_eval = 0.
        
        plantModel.An =comm.bcast( plantModel.AnCumul_inner/(dt*24*3600) , root=0)
        
        


def adapt_values( val_new, minVal, maxVal, volumes, divideEqually, verbose=False):
    """
    Update the concentration to get specific total content and gradient.
    @param val_new: total concentration of the element in the cylinder cells (cm3/cm3 water or mol/cm3 solute)
    @param minVal: min acceptable content or concentration (float or list)
    @param maxVal: max acceptable content or concentration (float or list)
    @param volumes: volume of each cell of the cylinder (cm3)
    @param divideEqually : divid excess or deficit element equally (True) or incrementally from cell nearest to the root surface outwards (False)
    @param verbose
    """
    if isinstance(minVal, numbers.Number):
        minVal=np.full(len(val_new), minVal)
    if isinstance(maxVal, numbers.Number):
        maxVal=np.full(len(val_new), maxVal)

    def print_verbose(*args):
        if verbose:
            print(rank, *args)
               
    def redistribute_excess(val_new):
        # Copy the initial values
        tempContent = val_new.copy()

        # Identify and handle excess values
        excess_indices = np.where(val_new > maxVal)
        excess_elements = (val_new[excess_indices] - maxVal[excess_indices]) * volumes[excess_indices]
        extraElement = np.sum(excess_elements)
        val_new[excess_indices] = maxVal[excess_indices]

        # Redistribution loop
        n_iter = 0
        while extraElement > 1e-20 and n_iter < 5:
            n_iter += 1
            canAdd = np.maximum((maxVal - val_new) * volumes, 0.) # max(,0.) to avoid rounding error issues
            if np.sum(canAdd) < extraElement:
                raise ValueError("Not enough capacity to redistribute excess elements.")
            
            addable_indices = np.where(canAdd > 0)[0]

            if divideEqually:
                # Method 1: Divide equally
                addable_amount = extraElement / len(addable_indices)# mean amount to add
                to_addList = np.minimum(addable_amount, canAdd)
                val_new -= to_addList / volumes
                extraElement -= sum(to_addList) 
            else:
                # Method 2: Incrementally add to the lowest index first
                for idx in addable_indices:
                    if extraElement <= 0:
                        break
                    to_add = min(extraElement, canAdd[idx])
                    val_new[idx] += to_add / volumes[idx]
                    extraElement -= to_add

        return val_new
        
    def redistribute_deficit(val_new):
        # Copy the initial values
        tempContent = val_new.copy()

        # Identify and handle deficit values
        deficit_indices = np.where(val_new < minVal)
        deficit_elements = (minVal[deficit_indices] - val_new[deficit_indices]) * volumes[deficit_indices]
        missingElement = np.sum(deficit_elements)
        val_new[deficit_indices] = minVal[deficit_indices]

        # Redistribution loop
        n_iter = 0
        while missingElement > 1e-20 and n_iter < 5:
            n_iter += 1
            print_verbose('helpfull::redistribute_deficit: add missing water', missingElement, val_new, minVal, maxVal)
            canTake = np.maximum((val_new - minVal) * volumes, 0)
            if np.sum(canTake) < missingElement:
                raise ValueError("Not enough capacity to redistribute deficit elements.")
            
            removable_indices = np.where(canTake > 0)[0]

            if divideEqually:
                # Method 1: Divide equally
                removable_amount = missingElement / len(removable_indices) # mean amount to take out
                to_removeList = np.minimum(removable_amount, canTake)
                val_new -= to_removeList / volumes
                missingElement -= sum(to_removeList) 
            else:
                # Method 2: Incrementally remove from the lowest index first
                for idx in removable_indices:
                    if missingElement <= 0:
                        break
                    to_remove = min(missingElement, canTake[idx])
                    print_verbose('helpfull::redistribute_deficit: not equally deficit', missingElement, val_new, idx, to_remove, removable_indices)
                    val_new[idx] -= to_remove / volumes[idx]
                    missingElement -= to_remove

        return val_new

    def handle_excess( val_new):
        if (val_new - maxVal > 0).any():
            if max(val_new - maxVal) < 1e-20:
                val_new = np.minimum(val_new, maxVal)
            else:
                val_new = redistribute_excess(val_new)
        return val_new

    def handle_deficit(val_new ):
        if (minVal - val_new > 0).any():
            if max(minVal - val_new) < 1e-20:
                val_new = np.maximum(val_new, minVal)
            else:
                val_new =  redistribute_deficit(val_new)
        return val_new
        
    val_newBU = val_new
    n_iter = 0
    while (max(val_new - maxVal) > 0 or max(minVal - val_new) > 0) and n_iter <= 5:
        n_iter += 1
        if (val_new - maxVal > 0).any():
            val_new = handle_excess(val_new)
        if (minVal - val_new > 0).any():
            val_new = handle_deficit(val_new)
    print_verbose('helpfull:adapt_values',repr(val_new),repr(val_newBU),repr(volumes), 
            sum(val_new) - sum(val_newBU), n_iter)
    assert abs(sum(val_new*volumes) - sum(val_newBU*volumes)) < 1e-20
    assert (val_new >= minVal).all()
    assert (val_new <= maxVal).all()
    return val_new


def distributeValSolute_(seg_values_content, volumes, source, dt, verbose: bool):

    # adapt solute sink if needed
    sourceInit = source
    source = max(-sum(seg_values_content)/dt,source)
    if abs((sourceInit - source)/sourceInit)*100 > 1.:
        verbose =True
        print("distributeVals sourceInit > source", 'dt',dt,
              'old source',sourceInit ,'new source',source,
              'Qsolute',seg_values_content,'cylVol',
             volumes)
        
    if (sum(seg_values_content) == 0.):# should normally only happen with source >= 0
        weightVals =np.full(len(seg_values_content), 1 /len(seg_values_content))
    elif source < 0:# goes away from the 1d models

        if min(seg_values_content)<0:
            seg_values_content = adapt_values(
                seg_values_content/volumes,
                        0,
                        np.Inf, 
                        volumes, divideEqually = False, 
                        verbose= verbose) *volumes

        weightVals = seg_values_content /sum(seg_values_content)
        if verbose:
            print("sum(seg_values[segIds])", seg_values, weightVals)
            
    else:# goes toward  the 1d models
        if min(abs(seg_values_content)) == 0.:# at least 1 seg_values = 0 but not all
            seg_values_content = np.maximum(seg_values_content,1.e-14)
            assert min(abs(seg_values_content)) != 0.
        weightVals = (1 / seg_values_content) / sum(1/seg_values_content)
    return weightVals * source, source
    
    
def distributeValWater_(seg_values_perVol_, volumes, source, dt, vg_soil, theta_wilting_point, verbose):    
    if verbose:            
        print('rhichardnompi::distribVals: before adapt',seg_values_perVol_, sum(seg_values_perVol_*cylVol))
    
    if (source > 0):# space available
        availableSpaceOrWater = (vg_soil.theta_S-seg_values_perVol_)*volumes
        # adapt source to not add more than the 1ds can manage:
        source = min(sum(availableSpaceOrWater)/dt,source)
    else:
        availableSpaceOrWater = (seg_values_perVol_- theta_wilting_point)*volumes
        # adapt source to not take more than the 1ds can manage:
        sourceInit = source
        source = max(-sum(availableSpaceOrWater)/dt,source)
        if abs((sourceInit - source)/sourceInit)*100 > 1.:
            verbose =True
            print("distributeVals sourceInit > source", 'dt',dt,
                  'old source',sourceInit ,'new source',source,
                  'theta',seg_values_perVol_,'cylVol',
                 volumes,'availableSpaceOrWater',availableSpaceOrWater,
                  'vg_soil.theta_S',vg_soil.theta_S,
                  'theta_wilting_point',theta_wilting_point)
            #raise Exception

    if sum(availableSpaceOrWater) > 0:
        if verbose:
            print('helpfull::distributeValWater_: after availableSpaceOrWater',
                  repr(availableSpaceOrWater))                   
        if min(availableSpaceOrWater)<0:
            availableSpaceOrWater = adapt_values(availableSpaceOrWater/volumes,
                        0,
                        np.Inf,#/dt/abs(source), 
                        volumes, divideEqually = False, 
                        verbose= verbose) *volumes


        if verbose:
            print('rhichardnompi::distribVals: after weightVals ', 
                        sum(availableSpaceOrWater /sum(availableSpaceOrWater)),
                  repr(availableSpaceOrWater /sum(availableSpaceOrWater)))

        return (availableSpaceOrWater /sum(availableSpaceOrWater)) * source, source # mol/day
    else: # no water or space available
        return np.full(len(availableSpaceOrWater),0.), source