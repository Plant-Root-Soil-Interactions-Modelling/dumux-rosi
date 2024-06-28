#!/bin/bash
#
#SBATCH --job-name=checksri
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=64
#SBATCH --nodes=1
#SBATCH --partition=cpu128
#SBATCH --time=20-00:00:00
#SBATCH --mem=100G
#SBATCH --mail-type=BEGIN,TIME_LIMIT_50,END,FAIL,ALL
#SBATCH --mail-user=m.giraud@fz-juelich.de
#SBATCH --output=./slurmOut/slurm-%j.out
#SBATCH --error=./slurmErr/slurm-%j.err

export OMPI_MCA_btl=^openib
module load openmpi/4.1.4

cd $HOME/DumuxDune27/DUMUX/dumux-rosi/python/paperSc 


#export NUMEXPR_MAX_THREADS=${SLURM_CPUS_PER_TASK}
#export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
DUMUX_NUM_THREADS=2 mpirun -n 64 python3 XcGrowth.py $1 dumux_10c $2 $3 $4 thr2

# $1 : start, $2: end, $3: param ind, $4 scenario, $5 noAds or nothing

# sbatch fpit1_10cind.sh 9.5 10.5 1640 none noAds
# sbatch fpit1_10cind.sh 9 25 1640 baseline noAds