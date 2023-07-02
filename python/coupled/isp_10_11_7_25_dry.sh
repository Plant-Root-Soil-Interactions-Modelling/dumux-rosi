#!/bin/bash
#
#SBATCH --job-name=isp_10_11_7_25_dry
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --partition=cpu256
#SBATCH --time=20-00:00:00
#SBATCH --mem=10G
#SBATCH --mail-type=BEGIN,TIME_LIMIT_50,END,FAIL,ALL
#SBATCH --mail-user=m.giraud@fz-juelich.de


cd $HOME/DUMUX/dumux-rosi/python/coupled

#python3 2pshallow.py
python3 uqrMaster2102_starchDynamic.py 10 11 7 25 dry ${SLURM_JOB_ID}