#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l mem=2GB
#PBS -l walltime=080:30:00
#PBS -N 4CSEQ100lenght


DATA_PATH="/home/users/dzisis/NewPilotExp"
cd $DATA_PATH

#Alignment with Bowtie 2
#chmod +x bowtie2-align

#max 2 mismatches, min score = .40*31 = 12.4
./bowtie2-align -q -x fragment_ends_AGATCT -a --score-min L,0,-0.4 -U INDEX_1_3_Perfect31.fq -S INDEX_1_3_Perfect_2mis.sam
./bowtie2-align -q -x fragment_ends_AGATCT -a --score-min L,0,-0.4 -U INDEX_1_4_Perfect29.fq -S INDEX_1_4_Perfect_2mis.sam
./bowtie2-align -q -x fragment_ends_AGATCT -a --score-min L,0,-0.4 -U INDEX_1_5_Perfect31.fq -S INDEX_1_5_Perfect_2mis.sam
./bowtie2-align -q -x fragment_ends_AGATCT -a --score-min L,0,-0.4 -U INDEX_1_6_Perfect29.fq -S INDEX_1_6_Perfect_2mis.sam

