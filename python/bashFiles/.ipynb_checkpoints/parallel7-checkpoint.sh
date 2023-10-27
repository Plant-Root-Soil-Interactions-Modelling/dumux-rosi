#!/bin/bash
#
#SBATCH --job-name=parallel7
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=7
#SBATCH --nodes=1
#SBATCH --partition=cpu256
#SBATCH --time=20-00:00:00
#SBATCH --mem=10G
#SBATCH --mail-type=BEGIN,TIME_LIMIT_50,END,FAIL,ALL
#SBATCH --mail-user=m.giraud@fz-juelich.de


cd $HOME/DUMUXexud/dumux-rosi/python/toymodel3_mpi


#export NUMEXPR_MAX_THREADS=${SLURM_CPUS_PER_TASK}
#export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
mpirun -n 7 python3 parallel.py