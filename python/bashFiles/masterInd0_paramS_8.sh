#!/bin/bash

for (( COUNTER=0; COUNTER < 99; COUNTER+=1 )); do
    sbatch fpit8cpu256.sh 10 14 $COUNTER none
done
# $1 : start, $2: end, $3: scenario, $4 noAds or nothing

# sbatch fpit2_10cind.sh 10 11 0 none
#https://ibg3113.ibg.kfa-juelich.de/hub/user-redirect/lab/tree/DumuxDune27/DUMUX/dumux-rosi/python/bashFiles/masterInd0_paramS.sh