#!/bin/bash

for i in {0..0}; do
    echo $i
    sbatch fpit1_10c95ind.sh $i
done