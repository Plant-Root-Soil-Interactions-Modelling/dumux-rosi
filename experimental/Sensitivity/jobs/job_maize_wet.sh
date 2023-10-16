#!/bin/bash                                                                                                                                                
#                                                                                                                                                          
#SBATCH --job-name=single_maize_wet
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --exclude=node02
#SBATCH --partition=cpu256
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH --mail-type=BEGIN,TIME_LIMIT_50,END
#SBATCH --mail-user=d.leitner@fz-juelich.de
 
cd ..
module load openmpi/4.1.4
mpirun python3 maize_water_wet.py