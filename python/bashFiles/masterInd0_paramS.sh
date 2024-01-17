#!/bin/bash

for (( COUNTER=0; COUNTER < 99; COUNTER+=1 )); do
    sbatch fpit2_10cind.sh 10 11 $COUNTER none
done
# $1 : start, $2: end, $3: scenario, $4 noAds or nothing

# sbatch fpit2_10cind.sh 10 11 0 none