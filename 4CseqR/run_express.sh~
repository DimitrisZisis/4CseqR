#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l mem=8GB
#PBS -l walltime=080:00:00
#PBS -N 4CexpressresultPilot


DATA_PATH="/home/users/dzisis/NewPilotExp"
cd $DATA_PATH

#Alignment with Express
#chmod +x express


./express --no-bias-correct -o xprs_out1_5_100length fragment_endsAllChrAGATCT100length.fa INDEX_1_5_Perfect_2mis100length_0insite.sam 
./express --no-bias-correct -o xprs_out1_6_100length fragment_endsAllChrAGATCT100length.fa INDEX_1_6_Perfect_2mis100length_0insite.sam 

exit

