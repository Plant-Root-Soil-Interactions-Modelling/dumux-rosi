#!/bin/bash
#
#SBATCH --job-name=CPB20R3contr
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --partition=cpu256
#SBATCH --time=20-00:00:00
#SBATCH --mem=10G
#SBATCH --mail-type=BEGIN,TIME_LIMIT_50,END,FAIL,ALL
#SBATCH --mail-user=m.giraud@fz-juelich.de
#SBATCH --output=./slurmOut/slurm-%j.out
#SBATCH --error=./slurmErr/slurm-%j.err


export OMPI_MCA_btl=^openib
module load openmpi/4.1.4

cd $HOME/dumux38bis/dumux/dumux-rosi/experimental/CPB20R3sim

#sbatch runSimulation.sh 20 21 getExud
python3 runSimulationControl2.py $1 $2 $3