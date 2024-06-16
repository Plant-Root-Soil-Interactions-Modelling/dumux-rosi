#!/bin/bash
#
#SBATCH --job-name=10cind
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=128
#SBATCH --nodes=1
#SBATCH --partition=cpu256-highmem
#SBATCH --time=20-00:00:00
#SBATCH --mem=200G
#SBATCH --mail-type=BEGIN,TIME_LIMIT_50,END,FAIL,ALL
#SBATCH --mail-user=m.giraud@fz-juelich.de
#SBATCH --output=./slurmOut/slurm-%j.out
#SBATCH --error=./slurmErr/slurm-%j.err


export OMPI_MCA_btl=^openib
module load openmpi/4.1.4

cd $HOME/DumuxDune27/DUMUX/dumux-rosi/python/paperSc 


#export NUMEXPR_MAX_THREADS=${SLURM_CPUS_PER_TASK}
#export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
DUMUX_NUM_THREADS=128 mpirun -n 128 python3 XcGrowth.py $1 dumux_10c $2 $3 $4 xx $5 $6

# $1 : start, $2: end, $3: param ind, $4 scenario, 
# optional : $5: spellStart, $6: spellduration

# sbatch fpit1_10cind.sh 9.5 10.5 1640 none noAds
# sbatch fpit1_10cind.sh 9 25 1640 baseline noAds
# sbatch fpit128cpu256.sh 9 9.06 76 customDry 9.02 0.02
# sbatch fpit128cpu256.sh 10 10.1 76 customDry 10.04 0.04