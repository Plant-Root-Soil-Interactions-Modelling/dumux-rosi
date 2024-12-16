import sys; sys.path.append("../../../BayesianOptimization")

import scipy
import os
import time
import numpy as np

from bayes_opt import BayesianOptimization
from bayes_opt.logger import JSONLogger
from bayes_opt.event import Events
from bayes_opt import acquisition
from bayes_opt.util import load_logs

import run_SA


def open_singleroot():

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

    load_logs(optimizer, logs = ["results/" + type_str + ".log"]);

    print("\n ... and the answer is")
    print(optimizer.max)

    p = optimizer.max["params"]
    print(p["kr_young"], p["kr_old"], p["kx_young"], p["kx_old"])

    # res = optimizer.space.res()
    # for r in res:
    #     print(r)


if __name__ == "__main__":

    open_singleroot()

