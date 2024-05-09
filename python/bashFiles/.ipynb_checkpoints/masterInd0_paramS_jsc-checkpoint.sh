#!/bin/bash

for (( COUNTER=0; COUNTER < 99; COUNTER+=1 )); do
    sbatch fpit130_10cind_jsc.sh 10 11 $COUNTER none
done
# $1 : start, $2: end, $3: scenario, $4 noAds or nothing

# bash masterInd0_1640_by164.sh 9.5 10.5 none noAds