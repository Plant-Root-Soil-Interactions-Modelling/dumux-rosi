#!/bin/bash
#
#SBATCH --job-name=10cGr95
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=60
#SBATCH --nodes=1
#SBATCH --partition=cpu256
#SBATCH --time=20-00:00:00
#SBATCH --mem=200G
#SBATCH --mail-type=BEGIN,TIME_LIMIT_50,END,FAIL,ALL
#SBATCH --mail-user=m.giraud@fz-juelich.de


cd $HOME/DUMUXexudRelease/DUMUX/dumux-rosi/python/fpit 


#export NUMEXPR_MAX_THREADS=${SLURM_CPUS_PER_TASK}
#export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
DUMUX_NUM_THREADS=60 mpirun -n 60 python3 XcGrowth.py 9.5 dumux_10c