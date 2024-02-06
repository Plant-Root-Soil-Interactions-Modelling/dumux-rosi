#!/bin/bash
#
#SBATCH --job-name=10cind
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=5
#SBATCH --nodes=1
#SBATCH --partition=cpu256
#SBATCH --time=20-00:00:00
#SBATCH --mem=10G
#SBATCH --mail-type=BEGIN,TIME_LIMIT_50,END,FAIL,ALL
#SBATCH --mail-user=m.giraud@fz-juelich.de

module load openmpi/4.1.4
cd $HOME/DumuxDune27/DUMUX/dumux-rosi/python/paperSc 


#export NUMEXPR_MAX_THREADS=${SLURM_CPUS_PER_TASK}
#export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
DUMUX_NUM_THREADS=5 mpirun -n 5 python3 XcGrowth.py $1 dumux_10c $2 $3 $4 $5

# $1 : start, $2: end, $3: param ind, $4 scenario, $5 noAds or nothing, $6 start spell, $7 end spell

# python3 XcGrowth.py 9 dumux_10c 9.06 0 customDry noAds 9.02 0.02
# sbatch fpit2_10cind.sh 9 10 0 customDry noAds 9.02 0.02
# sbatch fpit1_10cind.sh 9 25 1640 baseline noAds