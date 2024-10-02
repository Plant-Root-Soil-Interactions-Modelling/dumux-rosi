#!/bin/bash
#SBATCH --job-name=sra_springbarley_1D_hydrus_sandyloam_surface
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=48:00:00
#SBATCH --mem=200G
#SBATCH --partition=cpu256
python3 scenario_Axx.py springbarley 1D hydrus_sandyloam surface