#!/bin/bash


# pset 17
sbatch fpit128cpu128.sh 10 25 17 baseline
sbatch fpit128cpu128.sh 10 25 17 earlyDry
sbatch fpit128cpu128.sh 10 25 17 lateDry
# pset 38
sbatch fpit128cpu256mem.sh 10 25 38 baseline
sbatch fpit128cpu256mem.sh 10 25 38 earlyDry
sbatch fpit128cpu256mem.sh 10 25 38 lateDry 
# pset 44 
sbatch fpit128cpu128.sh 10 25 44 baseline
sbatch fpit128cpu256.sh 10 25 44 earlyDry
sbatch fpit128cpu256.sh 10 25 44 lateDry
# pset 85
sbatch fpit128cpu256.sh 10 25 85 baseline 
sbatch fpit128cpu256.sh 10 25 85 earlyDry
sbatch fpit128cpu256.sh 10 25 85 lateDry  

# sbatch fpit128cpu128.sh 10 25 17 none
# sbatch --nodelist=node02 fpit128cpu256.sh 10 25 44 earlyDry
# sbatch fpit128cpu128.sh 10 25 85 earlyDry 