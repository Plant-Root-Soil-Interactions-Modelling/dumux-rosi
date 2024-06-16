#!/bin/bash


# pset 5
sbatch fpit128cpu256.sh 10 25 5 earlyDry --nodelist=node08
sbatch fpit128cpu256.sh 10 25 5 lateDry --nodelist=node08
# pset 44
sbatch fpit128cpu256.sh 10 25 44 earlyDry --nodelist=node09
sbatch fpit128cpu256.sh 10 25 44 lateDry --nodelist=node09 

#sbatch fpit256cpu256.sh 10 25 44 earlyDry  
#sbatch fpit256cpu256.sh 10 25 44 lateDry  