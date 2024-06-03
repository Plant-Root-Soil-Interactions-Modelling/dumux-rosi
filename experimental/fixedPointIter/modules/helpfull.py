
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

class StdoutRedirector:
    def __init__(self, filepath):
        self.filepath = filepath
        self.libc = ctypes.CDLL(None)
        self.c_stdout = ctypes.c_void_p.in_dll(self.libc, "stdout")

    def __enter__(self):
        self.old_stdout_fd = os.dup(1)
        self.file = open(self.filepath, 'w')
        self.file_fd = self.file.fileno()
        os.dup2(self.file_fd, 1)

    def __exit__(self, exc_type, exc_val, exc_tb):
        os.dup2(self.old_stdout_fd, 1)
        os.close(self.old_stdout_fd)
        self.file.close()

def suggestNumStepsChange(dt, dt_inner, failedLoop,  targetIter_, results_dir):# taken from dumux
    """
    adapted from dumux function
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
    nOld = int(dt/dt_inner)
    if failedLoop:# if the inner computation failed, we also need to decrease the timestep
        numIter_ = perirhizalModel.k_iter
    else:
        numIter_ = nOld
        
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

        write_file_array("nNew_error", np.array([nNew , int(nNew), nOld,dt,dt_inner,
                                                 (nOld == int(nOld)), (nNew == int(nNew)),
                                                 real_dtinner ,dt, dt_inner , failedLoop,
                                                 abs((real_dtinner - dt)/dt*100.),rs_age]), 
                         directory_ =results_dir, fileType = '.csv') 

    dt_inner = dt/float(nNew)
    
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
        for extraText in ["","cyl_val/","printData/", "vtpvti/"]:
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

                    
def continueLoop(rs,n_iter, dt_inner: float,failedLoop: bool,
                         real_dtinner: float,name="continueLoop", isInner = False,doPrint = True, 
                         fileType = '.csv' , plant = None):
    results_dir = rs.results_dir
    plant.time_rhizo_cumul += plant.time_rhizo_i
    plant.time_3ds_cumul += plant.time_3ds_i
    plant.time_rhizo_i = 0
    plant.time_3ds_i = 0
    plant.time_plant_cumul = plant.time_plant_cumulW + plant.time_plant_cumulS 

    write_file_array("totalComputetime",
                     np.array([timeit.default_timer() - plant.time_start_global,
                    plant.time_plant_cumul,plant.time_rhizo_cumul ,plant.time_3ds_cumul]) , 
                     directory_ = results_dir)


    cL = ((np.floor(rs.err) > rs.max_err) or  rs.solve_gave_up or failedLoop
            or (np.floor(rs.diff1d3dCurrant_rel*1000.)/1000.>0.001) 
            or (np.floor(rs.maxdiff1d3dCurrant_rel*1000.)/1000.>0.001) 
            or (min(rs.new_soil_solute.reshape(-1)) < 0)  
            or ((n_iter < rs.minIter) and (isInner)))  and (n_iter < rs.k_iter)

    if (rank == 0) and isInner:
        print('continue loop?',rank,'n_iter',n_iter,'cL',cL,'failedLoop',failedLoop, ' np.floor(rs.err)',
        np.floor(rs.err),  'solve_gave_up',rs.solve_gave_up,
        'diff1d3dCurrant_rel',rs.diff1d3dCurrant_rel, 
        'maxdiff1d3dCurrant_rel',rs.maxdiff1d3dCurrant_rel, 
        'k_iter',rs.k_iter)

    cL = comm.bcast(cL,root = 0)
    failedLoop_ = np.array( comm.bcast(comm.gather(failedLoop,root = 0),root = 0))
    assert (failedLoop_ ==failedLoop_[0]).all() # all true or all false

    if doPrint:
        if not os.path.isfile(results_dir+name+fileType):
            write_file_array(name, np.array(['n_iter', 'err', 
                                             'diff1d3dCurrant_rel','maxdiff1d3dCurrant_rel',
                                             'solve_gave_up', 'min__soil_solute',
                                                     'dt_inner','real_dtinner','failedLoop','cL']), directory_ =results_dir, fileType = fileType)
            if not rs.doMinimumPrint:
                write_file_array(name+"Bool", np.array(['n_iter',  'err', 

                                                        'diff1d3dCurrant_rel','maxdiff1d3dCurrant_rel',
                                             'solve_gave_up',  'min__soil_solute',
                                                         'dt_inner','failedLoop','cL']), directory_ =results_dir, fileType = fileType)

        write_file_array(name, np.array([n_iter, rs.err, 
                                         rs.diff1d3dCurrant_rel,rs.maxdiff1d3dCurrant_rel,
                                         rs.solve_gave_up, min(rs.new_soil_solute.reshape(-1)),
                                                     dt_inner, real_dtinner,failedLoop,cL]), directory_ =results_dir, fileType = fileType)
        if not rs.doMinimumPrint:
            write_file_array(name+"2", 
                             np.concatenate((rs.sumDiff1d3dCW_rel,rs.maxDiff1d3dCW_rel)),  
                             directory_ =results_dir, fileType = fileType)
            write_file_array(name+"Bool", np.array([n_iter, (np.floor(rs.err) > rs.max_err), 
                                                    (np.floor(rs.diff1d3dCurrant_rel*10.)/10. > 0.1),
                                                    (np.floor(rs.maxdiff1d3dCurrant_rel) >1),
                                                    #(abs(rs.rhizoMassWError_abs) > 1e-13), (abs(rs.rhizoMassCError_abs) > 1e-9), 
                                                    #(max(abs(rs.errDiffBCs*0)) > 1e-5), 
                                                    rs.solve_gave_up, min(rs.new_soil_solute.reshape(-1)) < 0.,
                                                         dt_inner,failedLoop,cL]), directory_ =results_dir, fileType = fileType)

    return cL

def checkseg2cellMapping(seg2cell_old, plantModel):
    """
        some checks to do during the troubleshooting period
    """
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
    """
        save sumulative transpiration and assimilated carbon during inner fixed point iteration
    """
    plantModel.TranspirationCumul += plantModel.TranspirationCumul_inner 
    if perirhizalModel.doPhotosynthesis:
        if perirhizalModel.enteredSpell and (not perirhizalModel.leftSpell):
            plantModel.TranspirationCumul_eval += plantModel.TranspirationCumul_inner
        elif perirhizalModel.leftSpell:
            plantModel.TranspirationCumul_eval = 0.
        
        plantModel.An =comm.bcast( plantModel.AnCumul_inner/(dt*24*3600) , root=0)
