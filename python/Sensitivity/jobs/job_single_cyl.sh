#!/bin/bash                                                                                                                                                
#                                                                                                                                                          
#SBATCH --job-name=parallel-job
#SBATCH --ntasks=512
#SBATCH --nodes=2
#SBATCH --partition=cpu256
#SBATCH --time=10:00
#SBATCH --mem=2G
#SBATCH --mail-type=BEGIN,TIME_LIMIT_50,END
#SBATCH --mail-user=h.mustermann@fz-juelich.de
 
cd ..
module load openmpi/4.1.4
mpirun python3 maize_nitrate.py