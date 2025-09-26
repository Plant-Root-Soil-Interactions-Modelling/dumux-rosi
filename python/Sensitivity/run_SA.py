"""
    Dynamic:

    Sensitivity Analysis,
    starting multiple dynamic simulations of soybean or maize, either running locally (crashes for many jobs), or creating slurm files for the cluster
    see run_sra.py, 
    see __main__, arguments are passed over command line (freezing fixed parameters, and passing sensitive parameters) 
        
    On Cluster:
    source myenv/bin/activate
    pyhton3 run_SA.py (make sure target directory is empty)
    then use a bash file to send jobs to cluster (e.g. run_file.sh) 
        
    TODO json parameter file, and passing the hash, would be nicer 

    Daniel Leitner, 2025   
"""
import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src");

import numpy as np
import os
import random

import run_sra


def mid_(l):
    """ mid element of the list (ceiling)"""
    return l[len(l) // 2]


def make_lists(*args):
    """ turns single values into lists containing a single value """
    l = list(args)
    for i, v in enumerate(l):
        if isinstance(v, float):
            l[i] = [v]
    return (*l,)


def mkdir_p(dir):
    '''make a directory (dir) if it doesn't exist'''
    if not os.path.exists(dir):
        os.mkdir(dir)


def make_local(*args, ref = "mid"):
    """ creates the jobs for a local sensitivity analysis: 
        
        the input *args are lists representing value ranges, 
        values are varied for each range, using the mid values for the other ranges         
    """
    ranges = make_lists(*args)  # puts floats into single valued lists
    jobs = []

    c = 1  # job counter
    if ref == "mid":
        mids = [mid_(r) for r in ranges]
    elif ref == "left":
        mids = [r[0] for r in ranges]
    elif ref == "right":
        mids = [r[-1] for r in ranges]
    else:
        raise "make_local(): ref value unknown"

    jobs.append([c, *mids])

    for i, r in enumerate(ranges):
        if len(r) > 1:
            for v in r:
                c += 1
                mids_copy = mids.copy()
                mids_copy[i] = v
                jobs.append([c, *mids_copy])

    return np.array(jobs)


def start_jobs(type_str, file_name, root_type, enviro_type, sim_time, jobs, run_local):
    """starts the jobs calling run_sra.py with the arguments 
    type_str                 type of sensitivity analysis to let run_sra know how to interpet the data passed to it
    file_name                name of the sensitivity analysis, simulation number is added
    root_type                currently unused TODO to switch between soybean and maize
    enviro_type              choose parametrisation  
    sim_time                 simulation time
    jobs                     parameters for each simulation run (see make_local), 11 dimensions? 
    run_local                if True calls pyhton3, else for cluster calls sbatch
    """

    if type_str == "file":
        job_directory = "file"
    else:
        job_directory = os.path.join(os.getcwd(), file_name)

    if not run_local:
        print("Creating job files in folder:", job_directory, ", use bash script to send jobs to cluster with sbatch (e.g. run_file.sh)")
        mkdir_p(job_directory)

    for i, job in enumerate(jobs):

        try:
            job_name = file_name + str(int(job[0]))
        except:  # in case of "file" the filename is passed in job[0]
            job_name = file_name + str(job[0])

        # print("starting job <{:s}>:".format(type_str), job_name, enviro_type, sim_time, *job[1:])
        print("creating job <{:s}>:".format(type_str), "python3 run_sra.py {:s} {:s} {:g} {:g} {:g} {:g} {:g} {:g} {:g} {:g} {:g} {:g} {:g} {:g}".format(type_str, job_name, enviro_type, float(sim_time), *job[1:]))  # , sim_time)  #

        if run_local:
            # print("len", len(job[1:]))
            os.system("python3 run_sra.py {:s} {:s} {:g} {:g} {:g} {:g} {:g} {:g} {:g} {:g} {:g} {:g} {:g} {:g} & \n".format(type_str, job_name, enviro_type, float(sim_time), *job[1:]))  # , sim_time)  #

        else:
            if type_str == "file":  # for "file" add an index to the file name, for all others job id is stored in job[0] and added in L89
                job_file = os.path.join(job_directory, job_name + "_" + str(enviro_type) + ".job")
            else:
                job_file = os.path.join(job_directory, job_name + ".job")

            with open(job_file, 'w') as fh:
                fh.writelines("#!/bin/bash\n")
                fh.writelines("#SBATCH --job-name={:s}.job\n".format(job_name))
                fh.writelines("#SBATCH --ntasks=1\n")
                fh.writelines("#SBATCH --cpus-per-task=1\n")
                fh.writelines("#SBATCH --nodes=1\n")
                fh.writelines("#SBATCH --time=24:00:00\n")
                fh.writelines("#SBATCH --mem=16G\n")
                fh.writelines("#SBATCH --partition=cpu256\n")
                fh.writelines("\n")
                fh.writelines("cd ..\n")
                fh.writelines("source /etc/profile.d/modules.sh\n")
                fh.writelines("module load openmpi/4.1.4\n")
                fh.writelines("source myenv/bin/activate\n")
                fh.writelines("python3 run_sra.py {:s} {:s} {:g} {:g} {:g} {:g} {:g} {:g} {:g} {:g} {:g} {:g} {:g} {:g}\n".
                              format(str(type_str), str(job_name), int(enviro_type), float(sim_time), *job[1:]))
            # s.system("sbatch {:s}".format(job_file)) # DOES NOT WORK ANYMORE


def write_ranges(file_name, names, ranges):
    """ writes variable names and parameter ranges of the sensitivity analysis """
    assert len(names) == len(ranges), "run_SA.write_ranges(): len(names) != len(ranges) {:g}, {:g}".format(len(names), len(ranges))
    with open(file_name + "_range", 'w') as file:
        for i in range(0, len(ranges)):
            l = [str(r) for r in ranges[i]]
            range_ = ', '.join(l)
            file.write(names[i] + ", " + str(len(ranges[i])) + ", " + range_ + "\n")


def read_ranges(file_name):
    """ reads variable names and parameter ranges of the sensitivity analysis """
    names, ranges = [], []
    with open(file_name + "_range", 'r') as file:
       for line in file:
            entries = line.rstrip().split(", ")
            names.append(entries[0])
            n = int(entries[1])
            ranges.append([])
            for i in range(0, n):
                ranges[-1].append(float(entries[2 + i]))

    return names, ranges


def local_soybean():
    """ constructs a local sensitivity analysis presented in the preproject for values 
        "kr", "kx", "lmax1", "lmax2", "lmax3", "theta1", "a", "src"
    """
    print("local_soybean")
    type_str = "original"  # the 'original' sa analysis from the pre project
    root_type = "soybean"
    file_name = "local_soybean_new_"
    enviro_type = 0
    sim_time = 87.5  # 87.5  # 87.5  # days

    p1 = np.array([1.* 2 ** x for x in np.linspace(-1., 1., 9)])
    p2 = np.array([1.* 2 ** x for x in np.linspace(-2., 2., 9)])
    theta_ = np.linspace(0, np.pi / 2, 9)
    write_ranges("results/" + file_name,
                 ["kr", "kx", "lmax1", "lmax2", "lmax3", "theta1", "a", "src"],
                 [p2, p2, p1, p1, p1, theta_, p1, [2., 3, 4, 5]])
    jobs = make_local(p2 , p2, p1, p1, p1, theta_, 1., 1., p1, [2., 3, 4, 5])  #

    start_jobs(type_str, file_name, root_type, enviro_type, sim_time, jobs, run_local = False)


def local_soybean_new2():
    """ constructs a local sensitivity  """

    print("new local_soybean")
    type_str = "original_new2"
    root_type = "soybean"
    file_name = "local_soybean_new2_"
    enviro_type = 0
    sim_time = 1  # 87.5  # 87.5  # 87.5  # days

    p1 = np.array([1.* 2 ** x for x in np.linspace(-1., 1., 9)])
    write_ranges("results/" + file_name,
                 ["r", "r145", "r2", "r3", "ln", "ln145", "ln2", "a145", "a2", "a3"],
                 [p1, p1, p1, p1, p1, p1, p1, p1, p1, p1])
    jobs = make_local(p1, p1, p1, p1, p1, p1, p1, p1, p1, p1)

    print("number of jobs", len(jobs))

    start_jobs(type_str, file_name, root_type, enviro_type, sim_time, jobs, run_local = False)


def local_soybean_conductivities():
    """ constructs a local sensitivity analysis of hydraulic conductivities 
       for young and old parts, dependent on root order (145, 2, 3)
    """
    print("local_soybean_conductivities")
    type_str = "conductivities10"  # varying age dependent conductiviies
    root_type = "soybean"
    file_name = "local_soybean_conductivities_"
    enviro_type = 0
    sim_time = 87.5  # days

    kx = [0.1, 1.e-3, 1.e-3]
    kx_old = [0.35, 0.015]

    kr = [1.e-3, 4.e-3, 4.e-3]
    kr_old = [5e-4, 0.0015]

    # p2 = np.array([1.* 2 ** x for x in np.linspace(-1., 1., 9)])
    p2 = np.array([1.* 2 ** x for x in np.linspace(-4., 4., 17)])
    write_ranges("results/" + file_name,
                 ["ykr1", "okr1", "ykr2", "okr2", "kr3_", "ykx1", "okx1", "ykx2", "okx2", "kx3_"],
                 [p2 * kr[0], p2 * kr_old[0], p2 * kr[1], p2 * kr_old[1], p2 * kr[2], p2 * kx[0], p2 * kx_old[0], p2 * kx[1], p2 * kx_old[1], p2 * kx[2]])
    jobs = make_local(p2 * kr[0], p2 * kr_old[0], p2 * kr[1], p2 * kr_old[1], p2 * kr[2], p2 * kx[0], p2 * kx_old[0], p2 * kx[1], p2 * kx_old[1], p2 * kx[2])

    start_jobs(type_str, file_name, root_type, enviro_type, sim_time, jobs, run_local = False)


def local_singleroot_conductivities():
    """ constructs a local sensitivity analysis of hydraulic conductivities 
       for young and old parts, dependent on root order (145, 2, 3)
    """
    print("local_singleroot_conductivities")
    type_str = "singleroot_conductivities10"  # varying age dependent conductiviies
    root_type = "soybean"
    file_name = "local_singleroot_conductivities64_"
    enviro_type = 0
    sim_time = 87.5  # days

    kx = np.array([0.1]) / 16
    kx_old = np.array([0.35]) / 16 / 4

    kr = np.array([1.e-3]) * 16 * 4
    kr_old = np.array([5e-4]) * 16 * 4

    # p2 = np.array([1.* 2 ** x for x in np.linspace(-1., 1., 9)])
    p2 = np.array([1.* 2 ** x for x in np.linspace(-2., 2., 35)])
    write_ranges("results/" + file_name,
                 ["ykr1", "okr1", "ykx1", "okx1"],
                 [p2 * kr[0], p2 * kr_old[0], p2 * kx[0], p2 * kx_old[0]])
    jobs = make_local(p2 * kr[0], p2 * kr_old[0], p2 * kx[0], p2 * kx_old[0], 0., 0., 0., 0., 0., 0.)

    start_jobs(type_str, file_name, root_type, enviro_type, sim_time, jobs, run_local = False)


def local_soybean_tropisms():
    """ constructs a local sensitivity analysis of the influence of tropisms
       alters number of trials n and mean alteration sigma, dependent on root order (145, 2, 3)
    """
    print("local_soybean_tropisms")
    type_str = "tropisms"
    root_type = "soybean"
    file_name = "local_soybean_tropisms_"
    enviro_type = 0
    sim_time = 87.5  # 87.5  # days

    sigma_ = np.linspace(0.1, 0.5, 5)
    n_ = np.linspace(0., 5., 9)

    write_ranges("results/" + file_name,
                 ["n145", "n2", "n3", "sigma145", "sigma2", "sigma3"],
                 [n_, n_, n_, sigma_, sigma_, sigma_])

    jobs = make_local(n_, n_, n_, sigma_, sigma_, sigma_, 0., 0., 0., 0.)  # currently we always pass 10 valeus to run_sra
    start_jobs(type_str, file_name, root_type, enviro_type, sim_time, jobs, run_local = False)


def local_soybean_radii():
    """ constructs a local sensitivity analysis of root radii and root hairs 
       radii, dependent on root order (145, 2, 3)
       root hair zone length, dependent on root order (145, 2, 3)
       root hair length, dependent on root order (145, 2, 3)
       root elongation zone is fixed to 0.5
    """
    print("local_soybean_radii")
    type_str = "radii"
    root_type = "soybean"
    file_name = "local_soybean_radii_"
    enviro_type = 0
    sim_time = 87.5  # 87.5  # days

    a145 = np.linspace(0.01, .5, 9)
    a2 = np.linspace(0.01, .1, 9)
    a3 = np.linspace(0.01, .1, 9)

    hz = np.linspace(0., 4., 9)
    hairsZone145, hairsZone2, hairsZone3 = hz, hz, hz

    hl = np.linspace(0.01, 0.99, 9)
    hairsLength145, hairsLength2, hairsLength3 = hl, hl, hl

    write_ranges("results/" + file_name,
                 ["a145", "a2", "a3", "hairsZone145", "hairsZone2", "hairsZone3", "hairsLength145", "hairsLength2", "hairsLength3"],
                 [a145, a2, a3, hz, hz, hz, hl, hl, hl])
    jobs = make_local(a145, a2, a3, hz, hz, hz, hl, hl, hl, 0., ref = "left")  # currently we always pass 10 values to run_sra
    jobs = np.array(jobs)

    start_jobs(type_str, file_name, root_type, enviro_type, sim_time, jobs, run_local = False)


def simulate_list():
    print("simulate_list")
    root_type = "soybean"  # unused
    list_filename = "data/my_pick.txt"

    with open(list_filename, "r", encoding = "utf-8") as file:
        lines = file.readlines()
    lines = [line.strip() for line in lines]
    print("number of lines in my_pick.txt:", len(lines))
    jobs = []
    for line in lines:  # line = experiment name
        jobs.append([line, float(0.), 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.])  # <- put lines here

    sim_time = 87.5  # 87.5  # 87.5  # days
    enviro_type = 0
    start_jobs("file", "", "", enviro_type, sim_time, jobs, run_local = False)
    # run_sra.run_soybean_exp(lines[1], enviro_type, sim_time, save_all = False)

    # enviro_type = 1
    # start_jobs(type_str, "", root_type, enviro_type, sim_time, jobs, run_local = False)
    # enviro_type = 5
    # start_jobs(type_str, "", root_type, enviro_type, sim_time, jobs, run_local = False)
    # enviro_type = 36
    # start_jobs(type_str, "", root_type, enviro_type, sim_time, jobs, run_local = False)
    # enviro_type = 59
    # start_jobs(type_str, "", root_type, enviro_type, sim_time, jobs, run_local = False)

    # TODO
    # run_experiment(exp_name, enviro_type, sim_time) using a params file
    #


if __name__ == "__main__":

    i = 1

    if i == 1:
        local_soybean()
    elif i == 2:
        local_soybean_conductivities()
    elif i == 3:
        local_soybean_tropisms()  # did not work, try again with larger sigma? and insertion angles?
    elif i == 4:
        local_soybean_radii()  # did not work, check plausibility
    elif i == 5:
        local_singleroot_conductivities()
    elif i == 6:
        simulate_list()
    elif i == 7:
        local_soybean_new2()
    else:
        raise("on no")

# def make_global(kr_, kx_, lmax1_, lmax2_, lmax3_, theta1_, r1_, r2_, a_, src_):
#     """ creates the jobs for a global sensitivity analysis  """
#
#     kr_, kx_, lmax1_, lmax2_, lmax3_, theta1_, r1_, r2_, a_, src_ = make_lists(kr_, kx_, lmax1_, lmax2_, lmax3_, theta1_, r1_, r2_, a_, src_)
#
#     # create global sa jobs
#     jobs = []
#     i = 1
#     for kr in kr_:
#         for kx in kx_:
#             for lmax1 in lmax1_:
#                 for lmax2 in lmax2_:
#                     for lmax3 in lmax3_:
#                         for theta1 in theta1_:
#                             for r1 in r1_:
#                                 for r2 in r2_:
#                                     for a in a_:
#                                         for src in src_:
#                                             i += 1
#                                             jobs.append([i, kr, kx, lmax1, lmax2, lmax3, theta1, r1, r2, a, src])
#     return jobs

# def local_maize():
#     root_type = "maize"
#     file_name = "local_maize"
#     enviro_type = 0
#     sim_time = 95
#
#     if rank == 0:
#         p1 = np.array([1.* 2 ** x for x in np.linspace(-1., 1., 9)])
#         p2 = np.array([1.* 2 ** x for x in np.linspace(-2., 2., 9)])
#         # p3 = np.array([1.* 2 ** x for x in np.linspace(-3., 3., 19)])
#         theta_ = np.linspace(0, np.pi / 2, 9)
#         write_ranges("results/" + file_name,
#                      ["kr", "kx", "lmax1", "lmax2", "lmax3", "theta1", "a", "delaySB"],
#                      [p2, p2, p1, p1, p1, theta_, p1, p1])
#         jobs = make_local(p2 , p2 , p1, p1, p1, theta_, 1., 1., p1, p1)
#
#     else:
#         jobs = None
#
#     jobs = comm.bcast(jobs, root = 0)
#     run_jobs(file_name, root_type, enviro_type, sim_time, jobs)

# def local1_maize():
#     root_type = "maize"
#     file_name = "local_timing"
#     enviro_type = 0
#     sim_time = 30  # 95
#
#     if rank == 0:
#         timings = np.linspace(1., 20., 20)
#         theta_ = theta = 85. / 180 * np.pi
#         write_ranges("results/" + file_name,
#                      ["kr"],
#                      [timings])
#         jobs = make_local(timings , 1. , 1., 1., 1., theta_, 1., 1., 1., 1.)
#
#     else:
#         jobs = None
#
#     jobs = comm.bcast(jobs, root = 0)
#     run_jobs(file_name, root_type, enviro_type, sim_time, jobs)
#
# def global1_soybean():
#     root_type = "soybean"
#     file_name = "global_soybean"
#     enviro_type = 0
#     sim_time = 30
#
#     if rank == 0:
#         theta_ = theta = 85. / 180 * np.pi
#         kx = np.array([10 ** x for x in np.linspace(-1., 1., 20)])
#         kr = np.array([10 ** x for x in np.linspace(-1., 1., 25)])
#         # kx = np.array([10 ** x for x in np.linspace(-2., 0., 25)])
#         # kr = np.array([10 ** x for x in np.linspace(-2., 0., 25)])
#         write_ranges("results/" + file_name,
#                      ["kr", "kx"],
#                      [ kr , kx ])
#         jobs = make_global(kr , kx , 1., 1., 1., theta_, 1., 1., 1., 1.)
#     else:
#         jobs = None
#     jobs = comm.bcast(jobs, root = 0)
#     run_jobs(file_name, root_type, enviro_type, sim_time, jobs)
#
#
# def global1_maize():
#     root_type = "maize"
#     file_name = "global_maize_inc2"
#     enviro_type = 0
#     sim_time = 30
#
#     if rank == 0:
#         theta_ = theta = 85. / 180 * np.pi
#         # kx = np.array([10 ** x for x in np.linspace(-1., 1., 20)])
#         # kr = np.array([10 ** x for x in np.linspace(-1., 1., 25)])
#         kr = np.array([10 ** x for x in np.linspace(0., 4., 50)])
#         kx = np.array([10 ** x for x in np.linspace(-1., 1., 10)])  # _inc
#         write_ranges("results/" + file_name,
#                      ["kr", "kx"],
#                      [ kr , kx ])
#         jobs = make_global(kr , kx , 1., 1., 1., theta_, 1., 1., 1., 1.)
#     else:
#         jobs = None
#     jobs = comm.bcast(jobs, root = 0)
#     run_jobs(file_name, root_type, enviro_type, sim_time, jobs)
