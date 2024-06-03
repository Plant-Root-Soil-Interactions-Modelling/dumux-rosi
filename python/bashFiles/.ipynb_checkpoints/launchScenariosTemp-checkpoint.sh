#!/bin/bash


# pset 44
sbatch fpit128cpu128.sh 10 25 44 baseline  
sbatch fpit128cpu256.sh 10 25 44 earlyDry --nodelist=node10
sbatch fpit128cpu256.sh 10 25 44 lateDry --nodelist=node10
# pset 61
sbatch fpit128cpu256mem.sh 10 25 61 baseline  
sbatch fpit128cpu256mem.sh 10 25 61 earlyDry 
sbatch fpit128cpu256mem.sh 10 25 61 lateDry  

