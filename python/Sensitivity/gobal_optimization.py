import scipy
import os
import time
import numpy as np

from bayes_opt import BayesianOptimization
from bayes_opt.logger import JSONLogger
from bayes_opt.event import Events
from bayes_opt import acquisition

import run_SA


def start_objective(type_str, enviro_type, sim_time, job, run_local):
    """ starts a simulation on the cluster with sbatch (run_local = True), 
    or runs the Python3 simulation script on the local machine (run_local = False) """
    p = params_as_tuple(job, type_str, enviro_type, sim_time)
    file_name = type_str + "_" + str(hash(p))
    jobs = [p]
    root_type = 0  # unused
    if not finished_objective(type_str, enviro_type, sim_time, job):
        run_SA.start_jobs(type_str, file_name, root_type, enviro_type, sim_time, jobs, run_local)


def finished_objective(type_str, enviro_type, sim_time, job):
    """ checks if the computation is finished """
    p = params_as_tuple(job, type_str, enviro_type, sim_time)
    file_name = type_str + "_" + str(hash(p)) + str(0)
    found = os.path.isfile("results/" + file_name + ".npz")
    return found


def get_objective(type_str, enviro_type, sim_time, job):
    """ retrieves results from the output file """
    p = params_as_tuple(job, type_str, enviro_type, sim_time)
    file_name = type_str + "_" + str(hash(p)) + str(0)

    alldata = np.load("results/" + file_name + ".npz")

    times = alldata["times"]
    act_trans = alldata["act_trans"]
    cu = -scipy.integrate.simpson(act_trans, x = times)  # cumulative uptake over the whole simualtion period

    collar_pot_ = alldata["collar_pot"]
    mean_pot = np.mean(collar_pot_)  # mean collar potential over time

    return mean_pot  # cu, mean_pot


def run_optimizer(optimizer, type_str, enviro_type, sim_time, run_local):
    """
    starts the optimizer
    """
    iterations = 10
    points_per_iteration = 64

    for i in range(0, iterations):

        points = [optimizer.suggest() for j in range(0, points_per_iteration)]  # new sampling points

        print(points)
        for p in points:
            start_objective(type_str, enviro_type, sim_time, p, run_local)

        finished = False;
        while not finished:
            finished = True
            for p in points:
                finished = finished & finished_objective(type_str, enviro_type, sim_time, p)
                if not finished_objective(type_str, enviro_type, sim_time, p):
                    waiting_for = p.copy()
                    print("unfinished", finished, p, finished_objective(type_str, enviro_type, sim_time, p))
            if not finished:
                time.sleep(1)  # wait another minute

        results = []
        for p in points:
            results.append(get_objective(type_str, enviro_type, sim_time, p))

        for r, p in zip(results, points):
            optimizer.register(params = p, target = r)

# def run_otpimizer_async (if one is finished, start next suggestion)


def run_otpimizer_sbatch():
    pass


def params_as_tuple(p, type_str, enviro_type, sim_time):
    if type_str == "singleroot_conductivities10":
        return tuple([0., p["kr_young"], p["kr_old"], p["kx_young"], p["kx_old"], 0., 0., 0., 0., 0., 0., enviro_type, sim_time])  # first parameter is id in run_SA
    else:
        raise


def experiment_singleroot():

    type_str = "singleroot_conductivities10"  # change experimental setup in run_sra.py
    enviro_type = 0
    sim_time = 10.

    # acquisition_function = acquisition.UpperConfidenceBound(kappa = 1.)

    pbounds = {
        'kr_young': (1.e-6, 1.),
        'kr_old': (1.e-6, 1.),
        'kx_young': (1.e-6, 1.),
        'kx_old': (1.e-6, 1.),
        # "kr_age_young":(1., 30.),
        # "kx_age_young":(1., 30.)
        }

    optimizer = BayesianOptimization(
        f = None,  # dummy (never called)
        # acquisition_function = acquisition_function,
        pbounds = pbounds,
        verbose = 2,  # verbose = 1 prints only when a maximum is observed, verbose = 0 is silent
        random_state = 1,
    )

    logger = JSONLogger(path = "results/" + type_str + ".log")
    optimizer.subscribe(Events.OPTIMIZATION_STEP, logger)

    # job = {'kr_old': np.float64(1), 'kr_young': np.float64(2), 'kx_old': np.float64(3), 'kx_young': np.float64(4)}
    # start_objective(type_str, enviro_type, sim_time, job, True)
    # print(finished_objective(type_str, enviro_type, sim_time, job))

    run_optimizer(optimizer, type_str, enviro_type, sim_time, run_local = False)

    print("\n ... and the answer is")
    print(optimizer.max)

    p = optimizer.max["params"]
    print(p["kr_young"], p["kr_old"], p["kx_young"], p["kx_old"])


if __name__ == "__main__":

    experiment_singleroot()

