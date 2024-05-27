#!/bin/bash


# pset 5
sbatch fpit128cpu128.sh 10 25 5 baseline  
sbatch fpit128cpu128.sh 10 25 5 earlyDry
sbatch fpit128cpu128.sh 10 25 5 lateDry
# pset 61
#sbatch fpit128cpu128.sh 10 25 61 baseline  
#sbatch fpit128cpu256mem.sh 10 25 61 earlyDry 
#sbatch fpit128cpu256mem.sh 10 25 61 lateDry  