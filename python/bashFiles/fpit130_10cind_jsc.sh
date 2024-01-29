#!/bin/bash
#SBATCH --job-name=selectPset
#SBATCH -A esmtst
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=96
#SBATCH --cpus-per-task=1 
#SBATCH --output=mpi-%j.out
#SBATCH --error=mpi-%j.err
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=m.giraud@fz-juelich.de

module --force purge
module load Stages/2023
module load GCC
module load ParaStationMPI
# module load GCCcore/.10.3.0
module load Python 
module load OpenCV
module load CMake
module load SciPy-Stack
module load VTK


echo ${SLURM_CPUS_PER_TASK}
echo "Starting at `date`"
echo "Running on hosts: $SLURM_NODELIST"
echo "Running on $SLURM_NNODES nodes."
echo "Running $SLURM_NTASKS tasks."
cd /p/project/cjicg41/mgiraud/DuMux/DUMUX/dumux-rosi/python/paperSc

echo "Current working directory is `pwd`"

DUMUX_NUM_THREADS=192 srun -n 192 python3 XcGrowth.py $1 dumux_10c $2 $3 $4 $5

# $1 : start, $2: end, $3: param ind, $4 scenario, $5 noAds or nothing

# sbatch fpit1_10cind.sh 9.5 10.5 1640 none noAds
# sbatch fpit1_10cind.sh 9 25 1640 baseline noAds


#    dirPart1= "/p/scratch/cjicg41/mgiraud/results1d3d/addedPrints"
#    results_dir=dirPart1 +"/"+extraName+str(spellData['scenario']) \
#    +str(paramIndx_)+str(int(mpiVerbose))+l_ks+mode\
#                +"_"+str(initsim)+"to"+str(simMax)\
#                    +"_"+str(int(dt*24*60))+"mn_"\
#                    +str(int((dt*24*60 - int(dt*24*60))*60))+"s_"\
#                    +str(max_rank)+"_"+str(abs(p_mean))+"/"