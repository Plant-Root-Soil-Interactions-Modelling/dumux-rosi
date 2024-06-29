#!/bin/bash


# pset 5
sbatch rhizo32cpu128.sh 10 25 5 baseline  
sbatch rhizo32cpu128.sh 10 25 5 earlyDry
sbatch rhizo32cpu128.sh 10 25 5 lateDry
# pset 44
sbatch rhizo32cpu128.sh 10 25 44 baseline  
sbatch rhizo32cpu128.sh 10 25 44 earlyDry  
sbatch rhizo32cpu128.sh 10 25 44 lateDry 
# pset 61
sbatch rhizo32cpu128.sh 10 25 61 baseline
sbatch rhizo32cpu128.sh 10 25 61 earlyDry
sbatch rhizo32cpu128.sh 10 25 61 lateDry

#sbatch rhizo4cpu256.sh 10 25 44 none