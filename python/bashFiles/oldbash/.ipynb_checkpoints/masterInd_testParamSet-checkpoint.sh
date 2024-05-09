#!/bin/bash

for (( COUNTER=0; COUNTER<=1640; COUNTER+=164 )); do
    sbatch fpit1_10cind.sh $1 $2 $COUNTER $3 $4
done
# $1 : start, $2: end, $3: scenario, $4 noAds or nothing

# bash masterInd0_1640_by164.sh 9.5 10.5 none noAds