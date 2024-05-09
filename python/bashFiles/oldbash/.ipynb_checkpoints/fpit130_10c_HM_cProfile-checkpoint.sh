#!/bin/bash
#
#SBATCH --job-name=cProf
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=130
#SBATCH --nodes=1
#SBATCH --partition=cpu256-highmem
#SBATCH --time=20-00:00:00
#SBATCH --mem=10G
#SBATCH --mail-type=BEGIN,TIME_LIMIT_50,END,FAIL,ALL
#SBATCH --mail-user=m.giraud@fz-juelich.de


cd $HOME/DumuxDune27/DUMUX/dumux-rosi/python/paperSc 


#export NUMEXPR_MAX_THREADS=${SLURM_CPUS_PER_TASK}
#export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
DUMUX_NUM_THREADS=130 mpirun -n 130 python3 XcGrowth.py 10 dumux_10c 11 98 none cProfile

# $1 : start, $2: end, $3: param ind, $4 scenario, $5 noAds or nothing, $6 start spell, $7 end spell

# python3 XcGrowth.py 9 dumux_10c 10 0 customDry noAds 9.02 0.02