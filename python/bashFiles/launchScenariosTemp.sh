#!/bin/bash


sbatch fpit128cpu128.sh 10 25 19 earlyDry
sbatch fpit128cpu256.sh 10 25 47 earlyDry
sbatch fpit128cpu256.sh 10 25 76 earlyDry
sbatch fpit128cpu256.sh 10 25 83 earlyDry


# pset 83
#sbatch fpit128cpu256mem.sh 10 25 83 baseline 