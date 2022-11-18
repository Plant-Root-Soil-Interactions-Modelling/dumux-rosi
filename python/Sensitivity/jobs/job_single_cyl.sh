#!/bin/bash                                                                                                                                                
#                                                                                                                                                          
#SBATCH --job-name=single_rhizosphere
#SBATCH --ntasks=6
#SBATCH --nodes=2
#SBATCH --partition=cpu256
#SBATCH --time=10:00:00
#SBATCH --mem=16G
#SBATCH --mail-type=BEGIN,TIME_LIMIT_50,END
#SBATCH --mail-user=d.leitner@fz-juelich.de
 
cd ..
module load openmpi/4.1.4
mpirun python3 maize_nitrate.py