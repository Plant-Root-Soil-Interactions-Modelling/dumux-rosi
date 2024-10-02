#!/bin/bash
#SBATCH --job-name=sra_springbarley_3D_hydrus_loam_surface
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=48:00:00
#SBATCH --mem=200G
#SBATCH --partition=cpu256
python3 scenario_Axx.py springbarley 3D hydrus_loam surface