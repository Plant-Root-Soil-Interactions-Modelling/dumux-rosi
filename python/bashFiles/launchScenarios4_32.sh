#!/bin/bash


# pset 5 fpit4_32cpu256 fpit4_32cpu128
sbatch fpit4_32cpu128.sh 10 25 5 baseline  
sbatch fpit4_32cpu128.sh 10 25 5 earlyDry
sbatch fpit4_32cpu128.sh 10 25 5 lateDry
# pset 44
sbatch fpit4_32cpu128.sh 10 25 44 baseline  
sbatch fpit4_32cpu128.sh 10 25 44 earlyDry  
sbatch fpit4_32cpu128.sh 10 25 44 lateDry 
# pset 49 
#sbatch fpit128cpu128.sh 10 25 49 baseline
#sbatch fpit128cpu256mem.sh 10 25 49 earlyDry #--nodelist=node06
#sbatch fpit128cpu256.sh 10 25 49 lateDry --nodelist=node10
# pset 61
sbatch fpit4_32cpu128.sh 10 25 61 baseline  
sbatch fpit4_32cpu128.sh 10 25 61 earlyDry 
sbatch fpit4_32cpu128.sh 10 25 61 lateDry  

# sbatch fpit128cpu128.sh 10 25 17 none
# sbatch --nodelist=node02 fpit128cpu256.sh 10 25 49 earlyDry
# sbatch fpit128cpu128.sh 10 25 61 earlyDry 
# sbatch fpit128cpu128.sh 10 25 61 earlyDry

# sbatch fpit4_64cpu128.sh 10 25 5 lateDry