#!/bin/bash
#SBATCH --job-name=agg_maize_1D_hydrus_clay_surface
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=5:00:00
#SBATCH --mem=100G
#SBATCH --partition=cpu256
python3 scenario_agg.py maize 1D hydrus_clay surface