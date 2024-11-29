#!/bin/bash


# pset 5 fpit2_32cpu128TraiRhizo fpit4_32cpu128
#sbatch fpit2_32cpu128TraiRhizo.sh 10 25 5 baseline --nodelist=node13  
# pset 44
#sbatch fpit2_32cpu128TraiRhizo.sh 10 25 5 earlyDry --nodelist=node13   


# pset 5 fpit2_32cpu128TraiRhizo fpit4_32cpu128
# sbatch fpit2_32cpu128TraiRhizo.sh 10 25 5 baseline --nodelist=node13  
sbatch fpit2_32cpu128TraiRhizo.sh 10 25 5 earlyDry #--nodelist=node13  
sbatch fpit2_32cpu128TraiRhizo.sh 10 25 5 lateDry #--nodelist=node13  
# pset 44
sbatch fpit2_32cpu128TraiRhizo.sh 10 25 44 baseline #--nodelist=node16   
sbatch fpit2_32cpu128TraiRhizo.sh 10 25 44 earlyDry #--nodelist=node14    
sbatch fpit2_32cpu128TraiRhizo.sh 10 25 44 lateDry #--nodelist=node15   
# pset 49 
#sbatch fpit128cpu128.sh 10 25 49 baseline
#sbatch fpit128cpu256mem.sh 10 25 49 earlyDry #--nodelist=node06
#sbatch fpit128cpu256.sh 10 25 49 lateDry --nodelist=node10
# pset 61
sbatch fpit2_32cpu128TraiRhizo.sh 10 25 61 baseline #--nodelist=node16  
sbatch fpit2_32cpu128TraiRhizo.sh 10 25 61 earlyDry #--nodelist=node15   
sbatch fpit2_32cpu128TraiRhizo.sh 10 25 61 lateDry  #--nodelist=node15  


# sbatch fpit128cpu128.sh 10 25 17 none
# sbatch --nodelist=node02 fpit128cpu256.sh 10 25 49 earlyDry
# sbatch fpit128cpu128.sh 10 25 61 earlyDry 
# sbatch fpit128cpu128.sh 10 25 61 earlyDry

# sbatch fpit4_64cpu128.sh 10 25 5 lateDry