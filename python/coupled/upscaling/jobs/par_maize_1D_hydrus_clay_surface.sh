#!/bin/bash
#SBATCH --job-name=par_maize_1D_hydrus_clay_surface
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=48:00:00
#SBATCH --mem=200G
#SBATCH --partition=cpu256
python3 scenario_par2.py maize 1D hydrus_clay surface