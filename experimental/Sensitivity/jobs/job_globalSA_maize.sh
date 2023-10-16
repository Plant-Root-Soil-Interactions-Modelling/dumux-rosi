#!/bin/bash                                                                                                                                                
#                                                                                                                                                          
#SBATCH --job-name=global_maize
#SBATCH --ntasks=250
#SBATCH --cpus-per-task=1
#SBATCH --exclude=node02
#SBATCH --partition=cpu256
#SBATCH --time=48:00:00
#SBATCH --mem=0
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=d.leitner@fz-juelich.de
 
cd ..
module load openmpi/4.1.4
mpirun -n 250 python3 run_SA.py 4