"""
    job managment
    
    just put your plan into make_list()
    then run file __main__
    
    use start_jobs() for multiple slurm jobs (for cluster calling sbatch) 
    or run_jobs() for MPI                            
"""
import sys; sys.path.append("../../../../CPlantBox");  sys.path.append("../../../../CPlantBox/src")
sys.path.append("../../modules"); sys.path.append("../../../build-cmake/cpp/python_binding/");

from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank(); size = comm.Get_size()
import numpy as np
import os

import scenario_Axx as sra
import scenario_sra_old as sra_old
import scenario_Bxx as agg
import scenario_Cxx as par


def start_jobs(jobs):
    """ send as individual jobs """

    jobs = np.array(jobs)

    for job in jobs:

        method, plant, dim, soil, outer_method = job

        if method == "sra":
            py_name = "scenario_Axx.py"
        elif method == "sraOld":
            py_name = "scenario_sra_old.py"
        elif method == "agg":
            py_name = "scenario_Bxx.py"
        elif method == "par":
            py_name = "scenario_Cxx.py"
        else:
            raise("Unknown method" + method)

        job_name = method + "_" + plant + "_" + dim + "_" + soil + "_" + outer_method
        print(job_name)
        job_file = os.path.join("jobs", job_name + ".sh")

        with open(job_file, 'w') as fh:

            fh.writelines("#!/bin/bash\n")
            fh.writelines("#SBATCH --job-name={:s}\n".format(job_name))
            fh.writelines("#SBATCH --ntasks=1\n")
            fh.writelines("#SBATCH --nodes=1\n")
            fh.writelines("#SBATCH --time=48:00:00\n")
            fh.writelines("#SBATCH --mem=200G\n")
            fh.writelines("#SBATCH --partition=cpu256\n")
            # fh.writelines("#SBATCH --mail-type=BEGIN,TIME_LIMIT_50,END\n")
            # fh.writelines("#SBATCH --mail-user=d.leitner@fz-juelich.de\n")
            # fh.writelines("module load openmpi/4.1.4\n")
            fh.writelines("python3 {:s} {:s} {:s} {:s} {:s}".format(py_name, plant, dim, soil, outer_method))

        os.system("sbatch {:s}".format(job_file))

        # os.system("python3 run_sra.py {:s} {:g} {:g} {:g} {:g} {:g} {:g} {:g} {:g} {:g} {:g} {:g} {:g}\n".
        #                   format(job_name, enviro_type, sim_time, *job[1:]))


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
        elif method == "agg":
            agg.run_agg(sim_time, *jobs[i])
        elif method == "par":
            par.run_par(sim_time, *jobs[i])
        else:
            raise("Unknown method" + method)

    # if rank == 0:  # save index to values
    #     np.savetxt("results/" + file_name)

    comm.barrier()
    MPI.Finalize()


def make_list():
    jobs = []

    # all springbarley
    method = ["sra"]  # 'sra', sraOld, agg, par
    plant = ['maize']  # 'springbarley', 'soybean', 'maize'
    dim = ["3D"]  # "1D", "2D"
    soil = ['hydrus_loam', 'hydrus_clay', 'hydrus_sandyloam']  # 'hydrus_loam', 'hydrus_clay', 'hydrus_sandyloam'
    outer_radius = ['surface']  # 'length', 'surface', 'volume', 'voronoi'

    # method = ['agg']  # 'sra', sraOld, agg
    # plant = ['maize', 'springbarley']  # 'springbarley', 'soybean', 'maize'
    # dim = ["1D", "3D"]  # 1D, 3D
    # soil = ['hydrus_loam', 'hydrus_clay', 'hydrus_sandyloam']  # , 'hydrus_clay'
    # outer_radius = ['surface']

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

    sim_time = 14.5  # days

    if rank == 0:
        jobs = make_list()
    else:
        jobs = None

    jobs = comm.bcast(jobs, root = 0)
    start_jobs(jobs)  # sim_time is hardcoded in the __main__ parts
    # run_jobs(jobs, sim_time)

