"""
    starts a multiple sra simulations of soybean or maize and distribute over MPI ranks, 
    freezing fixed parameters, 
    and passing parameters for steady state analysis
"""
import sys; sys.path.append("../modules/"); sys.path.append("../../../CPlantBox/");  sys.path.append("../../../CPlantBox/src/python_modules")
sys.path.append("../../build-cmake/cpp/python_binding/")  # DUMUX solver

from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank(); size = comm.Get_size()
import numpy as np

import run_sra


def mid_(l):
    return l[len(l) // 2]


def run_all(global_, file_name, root_type, enviro_type, sim_time, kr_, kx_, lmax0_, lmax1_, lmax2_, theta0_, r0_, r1_, a_, src_):

    # make lists out of values
    if isinstance(kr_, float):
        kr_ = [kr_]
    if isinstance(kx_, float):
        kx_ = [kx_]
    if isinstance(lmax0_, float):
        lmax0_ = [lmax0_]
    if isinstance(lmax1_, float):
        lmax1_ = [lmax1_]
    if isinstance(lmax2_, float):
        lmax2_ = [lmax2_]
    if isinstance(theta0_, float):
        theta0_ = [theta0_]
    if isinstance(r0_, float):
        r0_ = [r0_]
    if isinstance(r1_, float):
        r1_ = [r1_]
    if isinstance(a_, float):
        a_ = [a_]
    if isinstance(src_, int):
        src_ = [src_]

    # create jobs
    jobs = []
    i = 0
    if global_:  # global sensitivity analysis
        for kr in kr_:
            for kx in kx_:
                for lmax0 in lmax0_:
                    for lmax1 in lmax1_:
                        for lmax2 in lmax2_:
                            for theta0 in theta0_:
                                for r0 in r0_:
                                    for r1 in r1_:
                                        for a in a_:
                                            for src in src_:
                                                i += 1
                                                jobs.append([i, kr, kx, lmax0, lmax1, lmax2, theta0, r0, r1, a, src])
    else:  # local sensitivity analysis
        i += 1
        jobs.append([i, mid_(kr_), mid_(kx_), mid_(lmax0_), mid_(lmax1_), mid_(lmax2_), mid_(theta0_), mid_(r0_), mid_(r1_), mid_(a_), mid_(src_)])
        if len(kr_) > 1:
            for kr in kr_:
                i += 1
                jobs.append([i, kr, mid_(kx_), mid_(lmax0_), mid_(lmax1_), mid_(lmax2_), mid_(theta0_), mid_(r0_), mid_(r1_), mid_(a_), mid_(src_)])
        if len(kx_) > 1:
            for kx in kx_:
                i += 1
                jobs.append([i, mid_(kr_), kx, mid_(lmax0_), mid_(lmax1_), mid_(lmax2_), mid_(theta0_), mid_(r0_), mid_(r1_), mid_(a_), mid_(src_)])
        if len(lmax0_) > 1:
            for lmax0 in lmax0_:
                i += 1
                jobs.append([i, mid_(kr_), mid_(kx_), lmax0, mid_(lmax1_), mid_(lmax2_), mid_(theta0_), mid_(r0_), mid_(r1_), mid_(a_), mid_(src_)])
        if len(lmax1_) > 1:
            for lmax1 in lmax1_:
                i += 1
                jobs.append([i, mid_(kr_), mid_(kx_), mid_(lmax0_), lmax1, mid_(lmax2_), mid_(theta0_), mid_(r0_), mid_(r1_), mid_(a_), mid_(src_)])
        if len(lmax2_) > 1:
            for lmax2 in lmax2_:
                i += 1
                jobs.append([i, mid_(kr_), mid_(kx_), mid_(lmax0_), mid_(lmax1_), lmax2, mid_(theta0_), mid_(r0_), mid_(r1_), mid_(a_), mid_(src_)])
        if len(theta0_) > 1:
            for theta0 in theta0_:
                i += 1
                jobs.append([i, mid_(kr_), mid_(kx_), mid_(lmax0_), mid_(lmax1_), mid_(lmax2_), theta0, mid_(r0_), mid_(r1_), mid_(a_), mid_(src_)])
        if len(r0_) > 1:
            for r0 in r0_:
                i += 1
                jobs.append([i, mid_(kr_), mid_(kx_), mid_(lmax0_), mid_(lmax1_), mid_(lmax2_), mid_(theta0_), r0, mid_(r1_), mid_(a_), mid_(src_)])
        if len(r1_) > 1:
            for r1 in r1_:
                i += 1
                jobs.append([i, mid_(kr_), mid_(kx_), mid_(lmax0_), mid_(lmax1_), mid_(lmax2_), mid_(theta0_), mid_(r0_), r1, mid_(a_), mid_(src_)])
        if len(a_) > 1:
            for a in a_:
                i += 1
                jobs.append([i, mid_(kr_), mid_(kx_), mid_(lmax0_), mid_(lmax1_), mid_(lmax2_), mid_(theta0_), mid_(r0_), mid_(r1_), a, mid_(src_)])
        if len(src_) > 1:
            for src in src_:
                i += 1
                jobs.append([i, mid_(kr_), mid_(kx_), mid_(lmax0_), mid_(lmax1_), mid_(lmax2_), mid_(theta0_), mid_(r0_), mid_(r1_), mid_(a_), src])

    # distribute jobs
    jobs = np.array(jobs)
    b = jobs.shape[0]
    num_per_rank = b // size

    lower_bound = rank * num_per_rank
    upper_bound = (rank + 1) * num_per_rank
    if rank == size - 1:
        upper_bound = b

    if rank == 0:
        print("Total number of jobs", b)
    print("Rank", rank, "does job numbers from", lower_bound, "to", upper_bound - 1, flush = True)

    # run jobs
    for i in range(lower_bound, upper_bound):
        if root_type == "soybean":
            run_sra.run_soybean(file_name + str(int(jobs[i][0])), enviro_type, sim_time, *jobs[i, 1:])
        else:
            raise("Unknown root type " + root_type)

    # save index to values
    if rank == 0:
        np.savetxt("results/" + file_name, jobs[:, 1:])


if __name__ == "__main__":

    root_type = "soybean"
    file_name = "test"
    enviro_type = 0
    sim_time = 25.
    p = np.array([1.* 2 ** x for x in np.linspace(-2., 2., 9)])
    # p = np.array([1. + x for x in np.linspace(-0.5, 0.5, 7)])
    print(p)
    kr = 1.e-5
    kx = 1.e-3
    run_all(True, file_name, root_type, enviro_type, sim_time, kr * 1. , kx * 1. , 1., 1., 1., 1., 1., 1., 1., [4])
    print("fin")

