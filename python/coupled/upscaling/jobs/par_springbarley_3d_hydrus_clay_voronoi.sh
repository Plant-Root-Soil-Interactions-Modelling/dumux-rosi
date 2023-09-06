#!/bin/bash
#SBATCH --job-name=par_springbarley_3d_hydrus_clay_voronoi
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=48:00:00
#SBATCH --mem=200G
#SBATCH --partition=cpu256
python3 scenario_par2.py springbarley 3d hydrus_clay voronoi