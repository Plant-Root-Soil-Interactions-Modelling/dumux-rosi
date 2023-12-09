#!/bin/bash

for (( COUNTER=0; COUNTER<=1640; COUNTER+=164 )); do
    sbatch fpit1_10c95_105ind.sh $COUNTER
done