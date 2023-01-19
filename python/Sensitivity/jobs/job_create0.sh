#!/bin/bash                                                                                                                                                
#                                                                                                                                                          
#SBATCH --job-name=create_envirotype0
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --exclude=node02
#SBATCH --partition=cpu256
#SBATCH --time=48:00:00
#SBATCH --mem=16G
#SBATCH --mail-type=BEGIN,TIME_LIMIT_50,END
#SBATCH --mail-user=d.leitner@fz-juelich.de
 
cd ..
module load openmpi/4.1.4
mpirun python3 create_envirotypes_tables.py