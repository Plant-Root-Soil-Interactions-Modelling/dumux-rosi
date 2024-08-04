#!/bin/bash
#
#SBATCH --job-name=4_32trR
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=32
#SBATCH --nodes=1
#SBATCH --partition=cpu256
#SBATCH --time=20-00:00:00
#SBATCH --mem=50G
#SBATCH --mail-type=BEGIN,TIME_LIMIT_50,END,FAIL,ALL
#SBATCH --mail-user=m.giraud@fz-juelich.de
#SBATCH --output=./slurmOut/slurm-%j.out
#SBATCH --error=./slurmErr/slurm-%j.err


export OMPI_MCA_btl=^openib
module load openmpi/4.1.4

cd $HOME/dumux38TraiRhizo/dumux/dumux-rosi/experimental/fixedPointIter2/scripts 

#echo $1 $2 $3 $4 $5 $6
#export NUMEXPR_MAX_THREADS=${SLURM_CPUS_PER_TASK}
#export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
DUMUX_NUM_THREADS=4 mpirun -n 32 python3 mainTraiRhizo.py $1 $2 $3 $4 $5 $6

# $1 : start, $2: end, $3: param ind, $4 scenario, 
# optional : $5: spellStart, $6: spellduration

# sbatch fpit1_10cind.sh 9.5 10.5 1640 none noAds
# sbatch fpit1_10cind.sh 9 25 1640 baseline noAds
# sbatch fpit128cpu256.sh 9 9.06 76 customDry 9.02 0.02
# sbatch fpit128cpu256.sh 10 10.1 76 customDry 10.04 0.04