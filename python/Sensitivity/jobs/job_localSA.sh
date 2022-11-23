#!/bin/bash                                                                                                                                                
#                                                                                                                                                          
#SBATCH --job-name=local_sensitivity_analysis
#SBATCH --ntasks=100
#SBATCH --nodes=1
#SBATCH --nodelist=node10
#SBATCH --partition=cpu256
#SBATCH --time=5:00:00
#SBATCH --mem=0
#SBATCH --mail-type=BEGIN,TIME_LIMIT_50,END
#SBATCH --mail-user=d.leitner@fz-juelich.de
 
cd ..
module load openmpi/4.1.4
mpirun -n 100 python3 run_SA.py