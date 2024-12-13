#!/bin/bash                                                                                                                                                
                                                                                                                                                       #SBATCH --job-name=optimize
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --exclude=node02
#SBATCH --cpus-per-task=1
#SBATCH --partition=cpu256
#SBATCH --time=48:00:00
#SBATCH --mem=8G
 
cd ..
python3 global_optimization.py

