#!/bin/bash


# pset 76 
sbatch fpit128cpu128.sh 10 25 76 baseline
sbatch fpit128cpu128.sh 10 25 76 earlyDry
sbatch fpit128cpu128.sh 10 25 76 lateDry
# pset 83
sbatch fpit128cpu256mem.sh 10 25 83 baseline 
sbatch fpit128cpu256mem.sh 10 25 83 earlyDry
sbatch fpit128cpu256mem.sh 10 25 83 lateDry  