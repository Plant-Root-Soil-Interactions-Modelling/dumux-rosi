#!/bin/bash
#SBATCH --job-name=selectPset
#SBATCH -A esmtst
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=130
#SBATCH --cpus-per-task=1 
#SBATCH --output=mpi-out.%j
#SBATCH --error=mpi-err.%j
#SBATCH --time=10:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=m.giraud@fz-juelich.de


module --force purge
module load Stages/2023
module load GCC/10.3.0
module load ParaStationMPI/5.4.10-1
module load GCCcore/.10.3.0
module load Python/3.10
module load OpenCV
module load CMake
module load SciPy-Stack
module load mpi4py

cd $HOME/DumuxDune27/DUMUX/dumux-rosi/python/paperSc 

echo ${SLURM_CPUS_PER_TASK}
echo "Starting at `date`"
echo "Running on hosts: $SLURM_NODELIST"
echo "Running on $SLURM_NNODES nodes."
echo "Running $SLURM_NTASKS tasks."
cd /p/project/cjicg41/mgiraud/DuMux/DUMUX/dumux-rosi/python/paperSc

echo "Current working directory is `pwd`"

DUMUX_NUM_THREADS=130 srun -n 130 python3 XcGrowth.py $1 dumux_10c $2 $3 $4 $5

# $1 : start, $2: end, $3: param ind, $4 scenario, $5 noAds or nothing

# sbatch fpit1_10cind.sh 9.5 10.5 1640 none noAds
# sbatch fpit1_10cind.sh 9 25 1640 baseline noAds