#!/bin/bash

for (( COUNTER=0; COUNTER<=1640; COUNTER+=164 )); do
    sbatch fpit1_10cind.sh 9.5 10.5 $COUNTER none noAds
done