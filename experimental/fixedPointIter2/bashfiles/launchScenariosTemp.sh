#!/bin/bash


# pset 5
sbatch fpit64cpu256.sh 10 25 5 earlyDry --nodelist=node04
sbatch fpit64cpu256.sh 10 25 5 lateDry --nodelist=node04
# pset 44
sbatch fpit64cpu128.sh 10 25 44 earlyDry 
sbatch fpit64cpu128.sh 10 25 44 lateDry 
# pset 61
sbatch fpit64cpu256.sh 10 25 61 earlyDry --nodelist=node06
sbatch fpit64cpu256.sh 10 25 61 lateDry --nodelist=node06

#sbatch fpit256cpu256.sh 10 25 44 earlyDry  
#sbatch fpit256cpu256.sh 10 25 44 lateDry  