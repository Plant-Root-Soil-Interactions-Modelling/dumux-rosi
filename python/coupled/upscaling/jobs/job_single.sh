#!/bin/bash                                                                                                                                                
#                                                                                                                                                          
#SBATCH --job-name=test_scenarios
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --exclude=node02
#SBATCH --partition=cpu256
#SBATCH --time=48:00:00
#SBATCH --mem=10G
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=d.leitner@fz-juelich.de
cd ..
python3 scenario_sra.py
