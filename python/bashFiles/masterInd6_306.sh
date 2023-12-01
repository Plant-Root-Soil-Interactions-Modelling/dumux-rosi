#!/bin/bash

for i in {6..306}; do
    echo $i
    sbatch fpit1_10c95ind.sh $i
done