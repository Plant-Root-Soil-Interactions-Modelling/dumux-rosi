#!/bin/bash
#
#SBATCH --job-name=3derr30
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=31
#SBATCH --nodes=1
#SBATCH --partition=cpu256
#SBATCH --time=20-00:00:00
#SBATCH --mem=10G
#SBATCH --mail-type=BEGIN,TIME_LIMIT_50,END,FAIL,ALL
#SBATCH --mail-user=m.giraud@fz-juelich.de


cd $HOME/DUMUXexudRelease/DUMUX/dumux-rosi/build-cmake/cpp/soil_richardsnc


DUMUX_NUM_THREADS=31 mpirun -np 31 ./richards10c_PyBi ./input/b3d_10c.input