#!/bin/bash


sbatch fpit100_10cind.sh 10 25 0 none noAds
sbatch fpit101_10cind.sh 10 25 0 none
sbatch fpit101_10cind.sh 10 25 0 baseline noAds 	
sbatch fpit101_10cind.sh 10 25 0 baseline	
sbatch fpit101_10cind.sh 10 25 0 earlyDry noAds 	
sbatch fpit128_10cind.sh 10 25 0 earlyDry 		
sbatch fpit128_10cind.sh 10 25 0 lateDry noAds 		
sbatch fpit101_10cind.sh 10 25 0 lateDry 	