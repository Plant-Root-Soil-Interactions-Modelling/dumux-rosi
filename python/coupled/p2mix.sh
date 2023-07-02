#!/bin/bash
#
#SBATCH --job-name=p2mix
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --partition=cpu256
#SBATCH --time=20-00:00:00
#SBATCH --mem=10G
#SBATCH --mail-type=BEGIN,TIME_LIMIT_50,END,FAIL,ALL
#SBATCH --mail-user=m.giraud@fz-juelich.de


cd $HOME/DUMUX/dumux-rosi/python/coupled

python3 2pcustom.py 14 30 wet mix ${SLURM_JOB_ID}