from mpi4py import MPI
from scipy.optimize import differential_evolution
import numpy as np

# MPI setup
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()


# Objective function
def objective_function(x):
    print(rank, ":", x)
    return np.sum(x ** 2)


# Define bounds
bounds = [(-5, 5)] * 10  # 10-dimensional problem


# Split work among ranks
def parallel_differential_evolution():
    # Each rank optimizes with a subset of the population
    popsize_per_rank = 15 // size  # Adjust population size per rank

    # Local optimization for each rank
    result = differential_evolution(
        objective_function,
        bounds,
        popsize = popsize_per_rank,
        seed = rank,
        disp = (rank == 0)  # Only display messages from rank 0
    )

    # Gather results from all ranks
    all_results = comm.gather(result, root = 0)

    if rank == 0:
        # Rank 0 selects the best solution from all results
        best_result = min(all_results, key = lambda res: res.fun)
        print(f"Best solution: {best_result.x}")
        print(f"Objective value: {best_result.fun}")


if __name__ == "__main__":
    parallel_differential_evolution()
