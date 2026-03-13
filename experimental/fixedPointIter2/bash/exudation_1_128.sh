#!/bin/bash
#
#SBATCH --job-name=fpitEx
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --partition=cpu128
#SBATCH --time=20-00:00:00
#SBATCH --mem=10G
#SBATCH --mail-type=BEGIN,TIME_LIMIT_50,END,FAIL,ALL
#SBATCH --mail-user=m.giraud@fz-juelich.de
#SBATCH --output=./slurmOut/slurm-%j.out
#SBATCH --error=./slurmErr/slurm-%j.err


export OMPI_MCA_btl=^openib
module load openmpi/4.1.4

cd $HOME/dumuxMagda/dumux/dumux-rosi/experimental/fixedPointIter2/scripts

source $HOME/cpbenv/bin/activate


export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
mpiexec -n 1 python3 mainExudate.py loam 2 5 60 True True True
