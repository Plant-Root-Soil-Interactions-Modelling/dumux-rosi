"""
    job managment
    
    just put your plan into make_list()
    then run file __main__
    
    use start_jobs() for multiple slurm jobs (for cluster calling sbatch)                     
"""
import sys; sys.path.append("../../../../CPlantBox");  sys.path.append("../../../../CPlantBox/src")
sys.path.append("../../modules"); sys.path.append("../../../build-cmake/cpp/python_binding/");
import numpy as np
import os


def start_jobs(jobs):
    """ send as individual jobs """
    jobs = np.array(jobs)

    for job in jobs:
        
        soil_type, sorption, SWP_ini, trans = job
        job_name = soil_type + "_" + sorption + "_"+ SWP_ini + "_"+ trans 
        print(job_name)
        job_file = os.path.join("jobs", job_name + ".sh")

        with open(job_file, 'w') as fh:

            fh.writelines("#!/bin/bash\n")
            fh.writelines("#SBATCH --job-name={:s}\n".format(job_name))
            fh.writelines("#SBATCH --ntasks=4\n")
            fh.writelines("#SBATCH --nodes=1\n")
            fh.writelines("#SBATCH --time=72:00:00\n")
            fh.writelines("#SBATCH --mem=50G\n")
            fh.writelines("#SBATCH --partition=cpu256\n")
            fh.writelines("source /etc/profile.d/modules.sh\n")
            fh.writelines("module load openmpi/4.1.4\n")
            fh.writelines("\n")
            fh.writelines("cd $HOME/Dumux_ExudatePaper\n")
            fh.writelines("source cpbenv_exudPaper/bin/activate\n")
            fh.writelines("cd $HOME/Dumux_ExudatePaper/dumux/dumux-rosi/experimental/fixedPointIter2/scripts\n")
           fh.writelines("mpiexec -n 1 python3 mainExudate.py {:s} {:s} {:s} {:s}".format(soil_type, sorption, SWP_ini, trans))
        os.system("sbatch {:s}".format(job_file))


def make_list():
    jobs = []

    # scenario static
    soil_type = ["loam", "sand"] 
    sorption = ["low", "medium", "high"]
    SWP_ini = ['100']
    trans = ['0']


    print("Creating", len(soil_type) * len(sorption) * len(SWP_ini)* len(trans), "simulations")
    print()
    
    if not os.path.exists('jobs'):
        os.makedirs('jobs')

    for s in soil_type:
        for p in sorption:
            for h in SWP_ini:
                for t in trans:
                    jobs.append([s, p, h, t])

    return jobs


if __name__ == "__main__":

    # good luck
    jobs = make_list()
    start_jobs(jobs)  