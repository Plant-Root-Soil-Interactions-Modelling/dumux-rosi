#!/bin/bash
#
#SBATCH --job-name=sobolPhlOd
#SBATCH --ntasks=256
#SBATCH --nodes=1
#SBATCH --partition=cpu256
#SBATCH --time=20-00:00:00
#SBATCH --mem=450G
#SBATCH --mail-type=BEGIN,TIME_LIMIT_50,END,FAIL,ALL
#SBATCH --mail-user=m.giraud@fz-juelich.de


cd $HOME/DUMUX/dumux-rosi/python/coupled

python3 sobolParallel_phloem.py 18 "dry"
#zip AllAuxC1 AllAuxC1/*.txt

#sbatch --nodelist=node03 All_C1.sh
