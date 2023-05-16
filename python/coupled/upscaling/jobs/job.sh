#!/bin/bash                                                                                                                                                
#                                                                                                                                                          
#SBATCH --job-name=upscaling_scenarios
#SBATCH --ntasks=128
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --exclude=node02
#SBATCH --partition=cpu256
#SBATCH --time=48:00:00
#SBATCH --mem=0
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=d.leitner@fz-juelich.de
cd ..
module load openmpi/4.1.4
mpirun -n 48 python3 jobs.py