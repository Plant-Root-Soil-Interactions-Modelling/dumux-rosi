
import numpy as np
import pandas as pd
import timeit
import pandas as pd
import matplotlib
import os
import timeit
import ctypes
import numbers
import tempfile
import signal
import sys
import threading
import multiprocessing
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
max_rank = comm.Get_size()


class TimeoutException(Exception):
    pass


def timeout_handler(signum, frame):
    raise TimeoutException


# Set the timeout handler
signal.signal(signal.SIGALRM, timeout_handler)

R_ph = 83.143  # perfect gas constant, hPa cm3K−1mmol−1


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
    except Exception as e:
        raise e  # halpe passe the error ?
    finally:
        # Disable the alarm
        signal.alarm(0)
    return result


def is_number(obj):
    # Check for standard Python numeric types (int, float) and NumPy numeric types
    return isinstance(obj, (numbers.Number, np.number))


class StdoutRedirector:

    def __init__(self, suppress_cerr = False):
        self.libc = ctypes.CDLL(None)
        self.c_stdout = ctypes.c_void_p.in_dll(self.libc, "stdout")
        self.buffer = None
        self.filepath = None
        if suppress_cerr:
            self.fd = 2
        else:
            self.fd = 1  # Do not suppress error messages

    def __enter__(self):
        self.old_stdout_fd = os.dup(self.fd)
        self.temp_file = tempfile.NamedTemporaryFile(delete = False)
        self.temp_file_fd = self.temp_file.fileno()
        os.dup2(self.temp_file_fd, self.fd)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        os.dup2(self.old_stdout_fd, self.fd)
        os.close(self.old_stdout_fd)
        self.temp_file.close()

        with open(self.temp_file.name, 'r') as temp_file:
            self.buffer = temp_file.read()

        os.remove(self.temp_file.name)

        if exc_type is not None:
            with open(self.filepath, 'w') as f:
                f.write(self.buffer)

        return False


def suggestNumStepsChange(dt, dt_inner, failedLoop, perirhizalModel, n_iter_inner_max):  # taken from dumux
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
    targetIter_ = perirhizalModel.targetIter  # objective max number of iteration for inner loop
    # dt_inner = perirhizalModel.dt_inner # time step of inner loop
    results_dir = perirhizalModel.results_dir
    nOld = int(dt / dt_inner)  # old number of step for inner loop

    # # adapt inner loop time step:
    if perirhizalModel.doNestedFixedPointIter:
        NinnerOld = int(np.ceil(dt_inner / perirhizalModel.dt_inner2))

    if failedLoop:  # if the inner computation failed, we also need to decrease the timestep
        numIter_ = perirhizalModel.k_iter  # k_iter == max accepted number of iteration for innerl loop
    else:
        numIter_ = n_iter_inner_max
    if (numIter_ > targetIter_):
        percent = float(numIter_ - targetIter_) / float(targetIter_)
        change = 1.0 + percent
    else:
        percent = float(targetIter_ - numIter_) / float(targetIter_)
        change = 1 / (1.0 + percent / 1.2)
    nNew = int(np.ceil(nOld * change))

    try:
        assert nOld == int(nOld)
        assert nNew == int(nNew)
    except:
        print('nNew iisue', nNew , int(nNew), nOld, dt, dt_inner,
              (nOld == int(nOld)), (nNew == int(nNew)))

    try:
        dt_inner = dt / float(nNew)
    except:
        write_file_array("suggestNumStepsChangeError", np.array([dt, dt_inner, failedLoop, n_iter_inner_max,
                                                                targetIter_, percent, change, nNew , nOld]),
                         directory_ = results_dir, fileType = '.csv')
        raise Exception

    write_file_array("suggestNumStepsChange", np.array([numIter_, n_iter_inner_max, targetIter_, percent, change,
                                                    nNew , nOld,
                                                 dt, dt_inner , failedLoop]),
                         directory_ = results_dir, fileType = '.csv')

    # # adapt inner loop time step:
    if perirhizalModel.doNestedFixedPointIter:
        perirhizalModel.dt_inner2 = dt_inner / NinnerOld
    else:
        perirhizalModel.dt_inner2 = dt_inner

    return dt_inner


def div0(a, b, c):  # function to avoid division by 0
    if isinstance(a, (list, np.ndarray)):
        return np.divide(a, b, out = np.full(len(a), c), where = b != 0)
    else:
        return div0f(a, b, c)


def div0f(a, b, c):  # function to avoid division by 0
    if b != 0:
        return a / b
    else:
        return a / c


def write_file_float(name, data, directory_, allranks = False):
    if (rank == 0) or allranks:
        os.makedirs(directory_, exist_ok = True)
        name2 = directory_ + name + '.txt'
        # print('write_file_float',name2, allranks)
        modeOp = 'a'
        if  False:  # 'fpit' in name:
            modeOp = 'w'
        with open(name2, modeOp) as log:
            # log.write(repr( data)  +'\n')
            np.savetxt(log, [data], delimiter = ',')


def write_file_array(name, data, space = ",", directory_ = "./results/", fileType = '.txt',
                     allranks = False):
    if (rank == 0) or allranks:
        try:
            os.makedirs(directory_, exist_ok = True)
            data = np.array(data, dtype = str).ravel()  # .reshape(-1)
            file_path = directory_ + name + fileType
            with open(file_path, "a") as log:
                np.savetxt(log, [data], delimiter = space, fmt = "%s")

        except:
            print('write_file_array', name, data, data.shape)
            raise Exception


def sinusoidal3(t, dt):
    """ sinusoidal function from 5:00 - 23:00, 0 otherwise (max is 1)"""
    return ((np.maximum(0., 0.75 + ((np.cos(2 * np.pi * (t - 0.5)) + np.cos(2 * np.pi * ((t + dt) - 0.5))) / 2)))) / 1.75


def getPsiAir(RH, TairC):  # constants are within photosynthesys.h
    Mh2o = 18  # molar weight, g mol-1, mg mmol-1

    rho_h2o = 1000  # water density mg cm-3
    return np.log(RH) * rho_h2o * R_ph * (TairC + 237.3) / Mh2o * (1 / 0.9806806)  ;  # in cm


def setDefault(s):
    """ Defined some usefull default parameters
    """
    s.setParameter("Problem.verbose", "0")
    s.setParameter("Newton.Verbosity", "0")

    # force solute mole fraction > 0 and water in possible pressure ranges
    s.setParameter("Newton.EnableChop", "true")

    # UpwindWeight = 1, better when we have high solute gradient.
    # UpwindWeight = 0.5, better when have high water flow and low solute gradient
    s.setParameter("Flux.UpwindWeight", "1")  # very important because we get high solute gradient.

    s.EnableResidualCriterion = False
    s.setParameter("Newton.EnableResidualCriterion",
                     str(s.EnableResidualCriterion))
    s.EnableAbsoluteResidualCriterion = False
    s.setParameter("Newton.EnableAbsoluteResidualCriterion",
                     str(s.EnableAbsoluteResidualCriterion))
    s.SatisfyResidualAndShiftCriterion = False
    s.setParameter("Newton.SatisfyResidualAndShiftCriterion",
                     str(s.SatisfyResidualAndShiftCriterion))
    s.MaxTimeStepDivisions = 100
    s.setParameter("Newton.MaxTimeStepDivisions",
                     str(s.MaxTimeStepDivisions))
    s.MaxSteps = 100
    s.setParameter("Newton.MaxSteps",
                     str(s.MaxSteps))
    # low MaxRelativeShift == higher precision in dumux
    s.MaxRelativeShift = 1e-10
    s.setParameter("Newton.MaxRelativeShift", str(s.MaxRelativeShift))

    return s
