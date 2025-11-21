#!/bin/bash

sbatch till1_4cpu128.sh 7.25 7.75 12
sbatch till1_13cpu128.sh 14.25 14.75 12
sbatch till1_33cpu128.sh 25.25 25.75 12

sbatch till1_4cpu128.sh 7.25 7.75 250
sbatch till1_13cpu128.sh 14.25 14.75 250
sbatch till1_33cpu128.sh 25.25 25.75 250

sbatch till1_4cpu128.sh 7.25 7.75 400
sbatch till1_13cpu128.sh 14.25 14.75 400
sbatch till1_33cpu128.sh 25.25 25.75 400

#sbatch till1_20cpu128.sh 7.25 7.75 12
#sbatch till1_60cpu128.sh 14.25 14.75 12
#sbatch till1_150cpu128.sh 25.25 25.75 12