#!/bin/bash

for i in {1640..1646}; do
    echo $i
    sbatch fpit1_10c95ind.sh $i
done