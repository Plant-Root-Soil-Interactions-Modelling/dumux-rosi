"""
    starts a multiple sra simulations of soybean or maize and distribute over MPI ranks, 
    freezing fixed parameters, 
    and passing parameters for steady state analysis
"""
import sys; sys.path.append("../modules/"); sys.path.append("../../../CPlantBox/");  sys.path.append("../../../CPlantBox/src/python_modules")
sys.path.append("../../build-cmake/cpp/python_binding/")  # DUMUX solver

from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank(); size = comm.Get_size()
import numpy as np
import os

import run_sra


def mid_(l):
    """ mid element of the list (floor)"""
    return l[len(l) // 2]


def make_lists(kr_, kx_, lmax1_, lmax2_, lmax3_, theta1_, r1_, r2_, a_, src_):
    """ turns single values into lists containing the value """
    if isinstance(kr_, float):
        kr_ = [kr_]
    if isinstance(kx_, float):
        kx_ = [kx_]
    if isinstance(lmax1_, float):
        lmax1_ = [lmax1_]
    if isinstance(lmax2_, float):
        lmax2_ = [lmax2_]
    if isinstance(lmax3_, float):
        lmax3_ = [lmax3_]
    if isinstance(theta1_, float):
        theta1_ = [theta1_]
    if isinstance(r1_, float):
        r1_ = [r1_]
    if isinstance(r2_, float):
        r2_ = [r2_]
    if isinstance(a_, float):
        a_ = [a_]
    if isinstance(src_, int):
        src_ = [src_]
    return kr_, kx_, lmax1_, lmax2_, lmax3_, theta1_, r1_, r2_, a_, src_


def run_jobs(file_name, root_type, enviro_type, sim_time, jobs):
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
        if (i < len(jobs)):
            if root_type == "soybean":
                run_sra.run_soybean(file_name + str(int(jobs[i][0])), enviro_type, sim_time, *jobs[i, 1:])
            else:
                raise("Unknown root type " + root_type)

    if rank == 0:  # save index to values
        np.savetxt("results/" + file_name, jobs[:, 1:])

    comm.barrier()
    MPI.Finalize()


def mkdir_p(dir):
    '''make a directory (dir) if it doesn't exist'''

    if not os.path.exists(dir):
        os.mkdir(dir)


def start_jobs(file_name, root_type, enviro_type, sim_time, jobs):
    """ send as individual jobs """

    job_directory = os.path.join(os.getcwd(), file_name)
    mkdir_p(job_directory)
    print(job_directory)
    jobs = np.array(jobs)

    for job in jobs:

        job_name = file_name + str(int(job[0]))
        print("Job", int(job[0]), ":", job_name, enviro_type, sim_time, *job[1:])
        job_file = os.path.join(job_directory, job_name + ".job")

        with open(job_file, 'w') as fh:

            fh.writelines("#!/bin/bash\n")
            fh.writelines("#SBATCH --job-name={:s}.job\n".format(job_name))
            fh.writelines("#SBATCH --ntasks=1\n")
            fh.writelines("#SBATCH --nodes=1\n")
            fh.writelines("#SBATCH --time=5:00:00\n")
            fh.writelines("#SBATCH --mem=2G\n")
            fh.writelines("#SBATCH --partition=cpu256\n")
            # fh.writelines("#SBATCH --mail-type=BEGIN,TIME_LIMIT_50,END\n")
            # fh.writelines("#SBATCH --mail-user=d.leitner@fz-juelich.de\n")
            # fh.writelines("module load openmpi/4.1.4\n")
            fh.writelines("python3 run_sra.py {:s} {:g} {:g} {:g} {:g} {:g} {:g} {:g} {:g} {:g} {:g} {:g} {:g}\n".
                          format(job_name, enviro_type, sim_time, *job[1:]))

        os.system("sbatch {:s}".format(job_file))
        # os.system("python3 run_sra.py {:s} {:g} {:g} {:g} {:g} {:g} {:g} {:g} {:g} {:g} {:g} {:g} {:g}\n".
        #                   format(job_name, enviro_type, sim_time, *job[1:]))


def make_local(kr_, kx_, lmax1_, lmax2_, lmax3_, theta1_, r1_, r2_, a_, src_):
    """ creates the jobs for a local sensitivity analysis  """

    kr_, kx_, lmax1_, lmax2_, lmax3_, theta1_, r1_, r2_, a_, src_ = make_lists(kr_, kx_, lmax1_, lmax2_, lmax3_, theta1_, r1_, r2_, a_, src_)

    # create local sa jobs
    jobs = []
    i = 1
    jobs.append([i, mid_(kr_), mid_(kx_), mid_(lmax1_), mid_(lmax2_), mid_(lmax3_), mid_(theta1_), mid_(r1_), mid_(r2_), mid_(a_), mid_(src_)])
    if len(kr_) > 1:
        for kr in kr_:
            i += 1
            jobs.append([i, kr, mid_(kx_), mid_(lmax1_), mid_(lmax2_), mid_(lmax3_), mid_(theta1_), mid_(r1_), mid_(r2_), mid_(a_), mid_(src_)])
    if len(kx_) > 1:
        for kx in kx_:
            i += 1
            jobs.append([i, mid_(kr_), kx, mid_(lmax1_), mid_(lmax2_), mid_(lmax3_), mid_(theta1_), mid_(r1_), mid_(r2_), mid_(a_), mid_(src_)])
    if len(lmax1_) > 1:
        for lmax1 in lmax1_:
            i += 1
            jobs.append([i, mid_(kr_), mid_(kx_), lmax1, mid_(lmax2_), mid_(lmax3_), mid_(theta1_), mid_(r1_), mid_(r2_), mid_(a_), mid_(src_)])
    if len(lmax2_) > 1:
        for lmax2 in lmax2_:
            i += 1
            jobs.append([i, mid_(kr_), mid_(kx_), mid_(lmax1_), lmax2, mid_(lmax3_), mid_(theta1_), mid_(r1_), mid_(r2_), mid_(a_), mid_(src_)])
    if len(lmax3_) > 1:
        for lmax3 in lmax3_:
            i += 1
            jobs.append([i, mid_(kr_), mid_(kx_), mid_(lmax1_), mid_(lmax2_), lmax3, mid_(theta1_), mid_(r1_), mid_(r2_), mid_(a_), mid_(src_)])
    if len(theta1_) > 1:
        for theta1 in theta1_:
            i += 1
            jobs.append([i, mid_(kr_), mid_(kx_), mid_(lmax1_), mid_(lmax2_), mid_(lmax3_), theta1, mid_(r1_), mid_(r2_), mid_(a_), mid_(src_)])
    if len(r1_) > 1:
        for r1 in r1_:
            i += 1
            jobs.append([i, mid_(kr_), mid_(kx_), mid_(lmax1_), mid_(lmax2_), mid_(lmax3_), mid_(theta1_), r1, mid_(r2_), mid_(a_), mid_(src_)])
    if len(r2_) > 1:
        for r2 in r2_:
            i += 1
            jobs.append([i, mid_(kr_), mid_(kx_), mid_(lmax1_), mid_(lmax2_), mid_(lmax3_), mid_(theta1_), mid_(r1_), r2, mid_(a_), mid_(src_)])
    if len(a_) > 1:
        for a in a_:
            i += 1
            jobs.append([i, mid_(kr_), mid_(kx_), mid_(lmax1_), mid_(lmax2_), mid_(lmax3_), mid_(theta1_), mid_(r1_), mid_(r2_), a, mid_(src_)])
    if len(src_) > 1:
        for src in src_:
            i += 1
            jobs.append([i, mid_(kr_), mid_(kx_), mid_(lmax1_), mid_(lmax2_), mid_(lmax3_), mid_(theta1_), mid_(r1_), mid_(r2_), mid_(a_), src])

    return jobs


def make_global(kr_, kx_, lmax1_, lmax2_, lmax3_, theta1_, r1_, r2_, a_, src_):
    """ creates the jobs for a global sensitivity analysis  """

    kr_, kx_, lmax1_, lmax2_, lmax3_, theta1_, r1_, r2_, a_, src_ = make_lists(kr_, kx_, lmax1_, lmax2_, lmax3_, theta1_, r1_, r2_, a_, src_)

    # create global sa jobs
    jobs = []
    i = 0
    if global_:
        for kr in kr_:
            for kx in kx_:
                for lmax1 in lmax1_:
                    for lmax2 in lmax2_:
                        for lmax3 in lmax3_:
                            for theta1 in theta1_:
                                for r1 in r1_:
                                    for r2 in r2_:
                                        for a in a_:
                                            for src in src_:
                                                i += 1
                                                jobs.append([i, kr, kx, lmax1, lmax2, lmax3, theta1, r1, r2, a, src])
    return jobs


if __name__ == "__main__":

    # root_type = "soybean"
    # file_name = "test"
    # enviro_type = 0
    # sim_time = 25.
    # p = np.array([1.* 2 ** x for x in np.linspace(-2., 2., 9)])
    # # p = np.array([1. + x for x in np.linspace(-0.5, 0.5, 7)])
    # print(p)
    # kr = 1.e-5
    # kx = 1.e-3
    # jobs = make_global(kr * 1. , kx * 1. , 1., 1., 1., 1., 1., 1., 1., [4])
    # print("fin")

    root_type = "soybean"
    file_name = "local_SA_test"
    enviro_type = 0
    sim_time = 87.5
    if rank == 0:
        p2 = np.array([1.* 2 ** x for x in np.linspace(-2., 2., 9)])
        kr = 1.e-4
        kx = 1.e-3
        p1 = np.array([1.* 2 ** x for x in np.linspace(-1., 1., 9)])
        theta_ = np.linspace(0, np.pi / 2, 9)
        jobs = make_local(kr * p2 , kx * p2 , p1, p1, p1, theta_, p1, p1, p1, [2, 3, 4, 5])
    else:
        jobs = None
    jobs = comm.bcast(jobs, root = 0)
    # start_jobs(file_name, root_type, enviro_type, sim_time, jobs)
    run_jobs(file_name, root_type, enviro_type, sim_time, jobs)
