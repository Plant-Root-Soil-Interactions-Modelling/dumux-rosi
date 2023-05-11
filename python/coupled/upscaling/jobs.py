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


def run_jobs(jobs):
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
            sra.run_sra(*jobs[i])
        elif method == "sraOld":
            sra_old.run_sraOld(*jobs[i])
        else:
            raise("Unknown method" + method)

    # if rank == 0:  # save index to values
    #     np.savetxt("results/" + file_name)

    comm.barrier()
    MPI.Finalize()


def make_list():
    jobs = []

    # jobs.append(['sra', 'maize', "3D", "hydrus_loam", "surface"])
    # jobs.append(['sra', 'soybean', "3D", "hydrus_loam", "surface"])
    # jobs.append(['sra', 'maize', "3D", "hydrus_clay", "surface"])
    # jobs.append(['sra', 'soybean', "3D", "hydrus_clay", "surface"])
    # jobs.append(['sra', 'maize', "3D", "hydrus_sandyloam", "surface"])
    # jobs.append(['sra', 'soybean', "3D", "hydrus_sandyloam", "surface"])

    jobs.append(['sra', 'maize', "1D", "hydrus_loam", "voronoi"])
    jobs.append(['sra', 'soybean', "1D", "hydrus_loam", "voronoi"])
    jobs.append(['sra', 'maize', "1D", "hydrus_clay", "voronoi"])
    jobs.append(['sra', 'soybean', "1D", "hydrus_clay", "voronoi"])
    jobs.append(['sra', 'maize', "1D", "hydrus_sandyloam", "voronoi"])
    jobs.append(['sra', 'soybean', "1D", "hydrus_sandyloam", "voronoi"])
    jobs.append(['sraOld', 'maize', "1D", "hydrus_loam", "voronoi"])
    jobs.append(['sraOld', 'soybean', "1D", "hydrus_loam", "voronoi"])
    jobs.append(['sraOld', 'maize', "1D", "hydrus_clay", "voronoi"])
    jobs.append(['sraOld', 'soybean', "1D", "hydrus_clay", "voronoi"])
    jobs.append(['sraOld', 'maize', "1D", "hydrus_sandyloam", "voronoi"])
    jobs.append(['sraOld', 'soybean', "1D", "hydrus_sandyloam", "voronoi"])

    return jobs


if __name__ == "__main__":

    jobs = make_list()
    run_jobs(jobs)

