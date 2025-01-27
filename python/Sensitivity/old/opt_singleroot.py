from mpi4py import MPI
from scipy.optimize import differential_evolution
import numpy as np
import timeit

import run_sra

# MPI setup
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

sim_time = 10.
enviro_type = 0
file_name = "opt_test"


# Objective function
def objective_function(x):
    """ Define objective function"""
    print("*********************************************")
    print("Rank", rank, ":", x)
    print("*********************************************")

    kr = np.zeros((3,))
    kr_old = np.zeros((2,))
    kx = np.zeros((3,))
    kx_old = np.zeros((2,))

    kr[0] = float(x[0])
    kr_old[0] = float(x[1])
    kx[0] = float(x[2])
    kx_old[0] = float(x[3])

    mods = {
    "filename": "data/Glycine_max_Moraes2020_singleroot.xml",
    "initial_age": 1.,
    "initial_totalpotential":-500
    # "domain_size": [10., 10., 200.]
    }
    try:
        cu = run_sra.run_soybean(file_name, enviro_type, sim_time, mods, kr, kx, kr_old, kx_old, save_all = False)
    except:
        cu = 1.e6  # bad

    return cu


def get_bounds():
    """ Define bounds"""
    # kx = np.array([0.1])
    # kx_old = np.array([0.35])
    # kr = np.array([1.e-3])
    # kr_old = np.array([5e-4])
    kr = [1.e-6, 1.]
    kr_old = [1.e-6, 1.]
    kx = [1.e-6, 1.]
    kx_old = [1.e-6, 1.]

    return [kr, kr_old, kx, kx_old]


# Split work among ranks
def parallel_differential_evolution():

    # Local optimization for each rank
    result = differential_evolution(
        objective_function,
        get_bounds(),
        popsize = 16,
        workers = 8,
        maxiter = 10
    )

    # Gather results from all ranks
    all_results = comm.gather(result, root = 0)

    if rank == 0:
        # Rank 0 selects the best solution from all results
        best_result = min(all_results, key = lambda res: res.fun)
        print(f"Best solution: {best_result.x}")
        print(f"Objective value: {best_result.fun}")


if __name__ == "__main__":

    start_time = timeit.default_timer()

    parallel_differential_evolution()

    print ("Optimization problem solved in ", timeit.default_timer() - start_time, " s")
