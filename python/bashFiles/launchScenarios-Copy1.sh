#!/bin/bash


# pset 7
sbatch fpit128cpu128.sh 10 25 7 baseline
sbatch fpit128cpu128.sh 10 25 7 earlyDry
sbatch fpit128cpu128.sh 10 25 7 lateDry
# pset 21
sbatch fpit128cpu256mem.sh 10 25 21 baseline
sbatch fpit128cpu256mem.sh 10 25 21 earlyDry
sbatch fpit128cpu256mem.sh 10 25 21 lateDry 
# pset 47 
sbatch fpit128cpu128.sh 10 25 47 baseline
sbatch fpit128cpu256mem.sh 10 25 47 earlyDry #--nodelist=node06
sbatch fpit128cpu256.sh 10 25 47 lateDry --nodelist=node10
# pset 85
sbatch fpit128cpu256.sh 10 25 85 baseline  --nodelist=node10
sbatch fpit128cpu256.sh 10 25 85 earlyDry --nodelist=node09
sbatch fpit128cpu256.sh 10 25 85 lateDry  --nodelist=node09

# sbatch fpit128cpu128.sh 10 25 17 none
# sbatch --nodelist=node02 fpit128cpu256.sh 10 25 47 earlyDry
# sbatch fpit128cpu128.sh 10 25 85 earlyDry 