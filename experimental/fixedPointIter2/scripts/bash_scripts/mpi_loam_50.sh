#!/bin/bash
#SBATCH --job-name=mpi_loam_test
#SBATCH --ntasks=4
#SBATCH --nodes=1
#SBATCH --time=48:00:00
#SBATCH --mem=200G
#SBATCH --partition=cpu256
source /etc/profile.d/modules.sh
module load openmpi/4.1.4

cd $HOME/Dumux
source cpbenv/bin/activate
cd $HOME/Dumux/dumux/dumux-rosi/experimental/fixedPointIter2/scripts
mpiexec -n 1 python3 launchExudate.py