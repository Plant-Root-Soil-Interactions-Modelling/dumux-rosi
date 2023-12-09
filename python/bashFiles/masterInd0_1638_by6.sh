#!/bin/bash

for (( COUNTER=0; COUNTER<=1638; COUNTER+=6 )); do
    sbatch fpit1_10c95_105ind.sh $COUNTER
done