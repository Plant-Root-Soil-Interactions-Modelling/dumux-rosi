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


cd $HOME/DUMUXexud/dumux-rosi/build-cmake/cpp/soil_richardsnc
./richards10c1d_cyl ./input/10cRcyl.input
#mpirun -n 2 richards10c1d_cyl ./input/10cRcyl.input
#mpirun -n 2 gdb -ex run --args richards10c1d_cyl ./input/10cRcyl.input