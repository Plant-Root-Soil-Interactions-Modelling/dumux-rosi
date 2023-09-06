#!/bin/bash
#SBATCH --job-name=agg_maize_3D_hydrus_sandyloam_surface
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=5:00:00
#SBATCH --mem=100G
#SBATCH --partition=cpu256
python3 scenario_agg.py maize 3D hydrus_sandyloam surface