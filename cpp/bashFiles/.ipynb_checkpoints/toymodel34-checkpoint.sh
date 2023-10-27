#!/bin/bash
#
#SBATCH --job-name=10cRcyl
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=2
#SBATCH --nodes=1
#SBATCH --partition=cpu256
#SBATCH --time=20-00:00:00
#SBATCH --mem=10G
#SBATCH --mail-type=BEGIN,TIME_LIMIT_50,END,FAIL,ALL
#SBATCH --mail-user=m.giraud@fz-juelich.de


cd $HOME/DUMUXexud/dumux-rosi/python/comparision
#mpirun -np 2 python3 toyModel3.4.py
mpiexec -n 2 python3 toyModel3.4.py 