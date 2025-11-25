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
        
        soil_type, res, simInit, simMax, exudate, decay = job
        job_name = soil_type + "_" + res + "_"+ simInit + "_"+ simMax +"_" + exudate+"_"+decay
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
            fh.writelines("cd $HOME/Dumux\n")
            fh.writelines("source cpbenv/bin/activate\n")
            fh.writelines("cd $HOME/Dumux/dumux/dumux-rosi/experimental/fixedPointIter2/scripts\n")
            fh.writelines("mpiexec -n 1 python3 mainExudate.py {:s} {:s} {:s} {:s} {:s}".format(soil_type, res, simInit, simMax, exudate, decay))
        os.system("sbatch {:s}".format(job_file))


def make_list():
    jobs = []

    # scenario static
    soil_type = ["loam"] #, "sand"] 
    res = ["5"] #, "4", "2", "1"] 
    simInit = ['10']
    simMax = ['60']
    exudate = ['True'] #, 'False']
    decay = ['True', 'False']

    print("Creating", len(soil_type) * len(res) * len(simInit)* len(simMax)* len(exudate)*len(decay), "simulations")
    print()
    
    if not os.path.exists('jobs'):
        os.makedirs('jobs')

    for s in soil_type:
        for r in res:
            for i in simInit:
                for m in simMax:
                    for e in exudate:
                        for d in decay:
                            jobs.append([s, r, i, m, e, d])

    return jobs


if __name__ == "__main__":

    # good luck
    jobs = make_list()
    start_jobs(jobs)  