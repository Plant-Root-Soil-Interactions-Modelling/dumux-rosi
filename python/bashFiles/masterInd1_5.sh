#!/bin/bash

for i in {1..5}; do
    echo $i
    sbatch fpit1_10c95ind.sh $i
done