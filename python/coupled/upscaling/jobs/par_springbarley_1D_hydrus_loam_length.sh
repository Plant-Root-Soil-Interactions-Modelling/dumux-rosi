#!/bin/bash
#SBATCH --job-name=par_springbarley_1D_hydrus_loam_length
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=48:00:00
#SBATCH --mem=200G
#SBATCH --partition=cpu256
python3 scenario_par2.py springbarley 1D hydrus_loam length