#!/bin/bash
#
#SBATCH --job-name=10cind
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --partition=cpu256
#SBATCH --time=20-00:00:00
#SBATCH --mem=10G
#SBATCH --mail-type=BEGIN,TIME_LIMIT_50,END,FAIL,ALL
#SBATCH --mail-user=m.giraud@fz-juelich.de


cd $HOME/DumuxDune27/DUMUX/dumux-rosi/python/paperSc 


#export NUMEXPR_MAX_THREADS=${SLURM_CPUS_PER_TASK}
#export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
DUMUX_NUM_THREADS=1 mpirun -n 1 python3 XcGrowth.py $1 dumux_10c $2 $3 $4 $5
#DUMUX_NUM_THREADS=1 mpirun -n 1 python3 XcGrowth.py 9.5 dumux_10c 10.5 $1 noAds
# $1 : start, $2: end, $3: param ind, $4 scenario, $5 noAds or nothing

# sbatch fpit1_10cind.sh 9.5 10.5 1640 none noAds
# sbatch fpit1_10cind.sh 9.5 10.5 1640 baseline noAds