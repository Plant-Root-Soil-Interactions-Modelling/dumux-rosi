#!/bin/bash


# pset 5
sbatch --nodelist=node10 rhizo32cpu256.sh 10 25 5 baseline   
sbatch rhizo32cpu256.sh 10 25 5 earlyDry --nodelist=node10
sbatch rhizo32cpu128.sh 10 25 5 lateDry
# pset 44
sbatch rhizo32cpu256.sh 10 25 44 baseline --nodelist=node10 
sbatch rhizo32cpu256.sh 10 25 44 earlyDry --nodelist=node10  
sbatch rhizo32cpu128.sh 10 25 44 lateDry 
# pset 61
sbatch rhizo32cpu256.sh 10 25 61 baseline --nodelist=node10
sbatch rhizo32cpu256.sh 10 25 61 earlyDry --nodelist=node10
sbatch rhizo32cpu256.sh 10 25 61 lateDry --nodelist=node10

#sbatch fpit256cpu256.sh 10 25 44 earlyDry  
#sbatch fpit256cpu256.sh 10 25 44 lateDry  