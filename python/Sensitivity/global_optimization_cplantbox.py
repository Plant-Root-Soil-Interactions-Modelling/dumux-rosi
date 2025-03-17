"""
    Bayesian optimization of the static root hydraulic model
"""

import sys; sys.path.append("../../../BayesianOptimization")

import scipy
import copy
import os
import time
import json
import timeit
import hashlib
import numpy as np

import asyncio

from bayes_opt import BayesianOptimization
from bayes_opt.logger import JSONLogger
from bayes_opt.event import Events
from bayes_opt import acquisition

import run_cplantbox


def hash_(name):
    hash_obj = hashlib.sha256(name.encode('utf-8'))
    r = hash_obj.hexdigest()
    return r


def background(f):

    def wrapped(*args, **kwargs):
        return asyncio.get_event_loop().run_in_executor(None, f, *args, **kwargs)

    return wrapped


@background
def start_objective(type_str, enviro_type, sim_time, initial, job):
    """ starts a simulation """
    file_name = type_str + "_" + hash_(json.dumps(job, sort_keys = True))
    print("start_objective()", file_name)
    if not finished_objective(type_str, enviro_type, sim_time, job):  # check if it already exists
        params = copy.deepcopy(job)
        params.update(initial)
        run_cplantbox.run_soybean(file_name, enviro_type, sim_time, params, False)
    else:
        print("start_objective()", "found", file_name)


def finished_objective(type_str, enviro_type, sim_time, job):
    """ checks if the computation is finished """
    file_name = type_str + "_" + hash_(json.dumps(job, sort_keys = True))
    found = os.path.isfile("results_cplantbox/" + file_name + ".npz")
    return found


def get_objective(type_str, enviro_type, sim_time, job):
    """ retrieves results from the output file """
    # file_name = type_str + "_" + hash_(json.dumps(job, sort_keys = True))
    # alldata = np.load("results_cplantbox/" + file_name + ".npz")
    # krs = alldata["krs"]
    return 0.  # krs[-1]  # return 0. for uniform(?) sampling


def run_optimizer(optimizer, type_str, enviro_type, sim_time, initial):
    """
    starts the optimizer
    """
    iterations = 5000
    points_per_iteration = 10

    for i in range(0, iterations):

        print("\nIteration", i, "\n")

        points = [optimizer.suggest() for j in range(0, points_per_iteration)]  # new sampling points

        for p in points:
            start_objective(type_str, enviro_type, sim_time, initial, p)

        finished = False;
        while not finished:
            finished = True
            for p in points:
                finished = finished & finished_objective(type_str, enviro_type, sim_time, p)
                if not finished:
                    waiting_for = p.copy()
                    print("run_optimizer() wating for", finished, p, finished_objective(type_str, enviro_type, sim_time, p))
                    print(type_str + "_" + str(hash(json.dumps(p, sort_keys = True))))
            if not finished:
                time.sleep(1)  # wait another minute

        results = []
        for p in points:
            results.append(get_objective(type_str, enviro_type, sim_time, p))

        # for r, p in zip(results, points):
        #     optimizer.register(params = p, target = r)


def cplantbox_tests():

    type_str = "soybean_test"  # change experimental setup in run_sra.py
    enviro_type = 0
    sim_time = 87.5

    initial = {
        "output_times": [40],
        "conductivity_mode": "scale",
        "scale_kr": 0.1,
        "scale_kx": 0.1
        }

    pbounds = {
        'src_a': (3, 11),
        'src_first_a': (3, 14),
        'src_delay_a': (3, 14),
         #
        'a1_a': (0.025, 0.25),
        'lmax1_a': (50, 150),
        'ln1_a': (0.2, 5.),
        'r1_a': (0.2, 7.),
        'tropismN1_a': (0., 7.),
        'hairsLength_a': (0.1, 2.),
        'hairsZone_a': (1., 10.),
        'hairsElongation_a': (0.1, 2.),
        }

    optimizer = BayesianOptimization(
        f = None,  # dummy (never called)
        pbounds = pbounds,
        verbose = 2,  # verbose = 1 prints only when a maximum is observed, verbose = 0 is silent
        random_state = 1,
    )

    logger = JSONLogger(path = "results_cplantbox/" + type_str + ".log")
    optimizer.subscribe(Events.OPTIMIZATION_STEP, logger)

    run_optimizer(optimizer, type_str, enviro_type, sim_time, initial)


def cplantbox_all14():  # 30 parameters

    type_str = "soybean_all14"
    enviro_type = 0
    sim_time = 87.5

    initial = {"output_times": [45, 50, 60],
              "conductivity_mode": "from_mecha",
              "mecha_path": "/home/daniel/Dropbox/granar/mecha_results",
              'hairsElongation1_a': 1.,
              'hairsElongation2_a': 1.,
              'hairsElongation3_a': 1.,
               }

    # acquisition_function = acquisition.UpperConfidenceBound(kappa = 1.)
    max_ind = 151
    pbounds = {

        'conductivity_index1': (0, max_ind),
        'conductivity_index2': (0, max_ind),
        'conductivity_index3': (0, max_ind),
        "conductivity_age1": (1, 21),
        "conductivity_age2": (1, 21),
        "conductivity_age3": (1, 21),

        'src_a': (3, 11),
        'src_first_a': (3, 14),
        'src_delay_a': (3, 14),

        'lmax145_a': (50, 150),
        'ln145_a': (0.5, 10.),
        'r145_a': (0.2, 7.),
        'theta145_a': (np.pi / 8., np.pi / 2.),
        'tropismN145_a': (0., 3.5),
        'hairsLength145_a': (0.1, 1.),
        'hairsZone145_a': (1., 10.),

        'lmax2_a': (5., 50.),
        'ln2_a': (0.5, 10.),
        'r2_a': (0.2, 7.),
        'theta2_a': (np.pi / 8., np.pi / 2.),
        'tropismN2_a': (0., 3.5),
        'hairsLength2_a': (0.1, 1.),
        'hairsZone2_a': (1., 10.),

        'lmax3_a': (5., 50.),
        'r3_a': (0.2, 7.),
        'theta3_a': (np.pi / 8., np.pi / 2.),
        'tropismN3_a': (0., 3.5),
        'hairsLength3_a': (0.1, 1.),
        'hairsZone3_a': (1., 10.),

        }

    optimizer = BayesianOptimization(
        f = None,  # dummy (never called)
        # acquisition_function = acquisition_function,
        pbounds = pbounds,
        verbose = 2,  # verbose = 1 prints only when a maximum is observed, verbose = 0 is silent
        random_state = 1,
    )

    logger = JSONLogger(path = "results_cplantbox/" + type_str + ".log")
    optimizer.subscribe(Events.OPTIMIZATION_STEP, logger)

    # job = {}
    # start_objective(type_str, enviro_type, sim_time, initial, job, True)
    # print(finished_objective(type_str, enviro_type, sim_time, job))

    run_optimizer(optimizer, type_str, enviro_type, sim_time, initial)

    print("\n ... and the answer is")
    print(optimizer.max)

    p = optimizer.max["params"]


def cplantbox_length14():

    type_str = "soybean_length14"
    enviro_type = 0
    sim_time = 87.5

    initial = {
        "output_times": [45, 50, 60],
        "conductivity_mode": "scale",
        "scale_kr": 0.1,
        "scale_kx": 0.1
        }

    pbounds = {
        'src_a': (3, 11),
        'src_delay_a': (3, 14),

        'lmax145_a': (50, 150),
        'ln145_a': (0.5, 10.),
        'r145_a': (0.2, 7.),

        'lmax2_a': (5., 50.),
        'ln2_a': (0.5, 10.),
        'r2_a': (0.2, 7.),

        'lmax3_a': (5., 50.),
        'r3_a': (0.2, 7.),
        }

    optimizer = BayesianOptimization(
        f = None,  # dummy (never called)
        pbounds = pbounds,
        verbose = 2,  # verbose = 1 prints only when a maximum is observed, verbose = 0 is silent
        random_state = 1,
    )

    logger = JSONLogger(path = "results_cplantbox/" + type_str + ".log")
    optimizer.subscribe(Events.OPTIMIZATION_STEP, logger)
    run_optimizer(optimizer, type_str, enviro_type, sim_time, initial)


if __name__ == "__main__":

    start_time = timeit.default_timer()

    # cplantbox_tests()
    cplantbox_all14()
    # cplantbox_length14()

    print ("\nThat took", timeit.default_timer() - start_time, " s")
