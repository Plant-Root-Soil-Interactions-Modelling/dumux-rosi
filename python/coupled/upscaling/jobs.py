"""
                 job managment
"""
import sys; sys.path.append("../../../../CPlantBox");  sys.path.append("../../../../CPlantBox/src")
sys.path.append("../../modules"); sys.path.append("../../../build-cmake/cpp/python_binding/");

from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank(); size = comm.Get_size()
import numpy as np
import os

import scenario_sra as sra
import scenario_sra_old as sra_old


def run_jobs(jobs, sim_time):
    """ distribute jobs to MPI ranks """

    jobs = np.array(jobs)
    b = jobs.shape[0]
    if rank == 0:
        print("Total number of jobs", b)

    num_per_rank = b // size
    num_per_rank = max(1, num_per_rank)

    lower_bound = rank * num_per_rank
    upper_bound = (rank + 1) * num_per_rank
    if rank == size - 1:
        upper_bound = b

    print("Rank", rank, "does job numbers from", lower_bound, "to", upper_bound - 1, flush = True)

    for i in range(lower_bound, upper_bound):  # run jobs

        method, plant, dim, soil, outer_method = jobs[i]

        if method == "sra":
            sra.run_sra(sim_time, *jobs[i])
        elif method == "sraOld":
            sra_old.run_sraOld(sim_time, *jobs[i])
        else:
            raise("Unknown method" + method)

    # if rank == 0:  # save index to values
    #     np.savetxt("results/" + file_name)

    comm.barrier()
    MPI.Finalize()


def make_list():
    jobs = []

    method = ['sra']
    plant = ['soybean', 'maize']
    dim = ['1D', '3D']  # 1D, 3D
    soil = ['hydrus_loam', 'hydrus_clay', 'hydrus_sandyloam']
    outer_radius = ['voronoi']

    print("Creating", len(method) * len(plant) * len(dim) * len(soil) * len(outer_radius), "simulations")
    print()

    for m in method:
        for p in plant:
            for d in dim:
                for s in soil:
                    for o in outer_radius:
                        jobs.append([m, p, d, s, o])
    return jobs


if __name__ == "__main__":

    sim_time = 7.5  # days

    if rank == 0:
        jobs = make_list()
    else:
        jobs = None

    jobs = comm.bcast(jobs, root = 0)
    run_jobs(jobs, sim_time)

