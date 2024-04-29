#!/bin/bash


# pset 19
sbatch fpit128cpu256.sh 10 25 19 baseline
sbatch fpit128cpu256.sh 10 25 19 earlyDry
sbatch fpit128cpu256.sh 10 25 19 lateDry
# pset 47
sbatch fpit128cpu256.sh 10 25 47 baseline
sbatch fpit128cpu256.sh 10 25 47 earlyDry
sbatch fpit128cpu256.sh 10 25 47 lateDry 
# pset 76 
sbatch fpit128cpu256.sh 10 25 76 baseline
sbatch fpit128cpu256.sh 10 25 76 earlyDry
sbatch fpit128cpu256.sh 10 25 76 lateDry
# pset 83
sbatch fpit128cpu256.sh 10 25 83 baseline 
sbatch fpit128cpu256.sh 10 25 83 earlyDry
sbatch fpit128cpu256.sh 10 25 83 lateDry  