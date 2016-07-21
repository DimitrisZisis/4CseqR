#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l mem=2GB
#PBS -l walltime=00:20:00
#PBS -N 4CSEQbuild

#BT=/home/users/dzisis/bedtools-2.17.0/bin/bedtools

DATA_PATH="/home/users/dzisis/NewPilotExp"
cd $DATA_PATH

./bowtie2-build --offrate 1 frag_ends_AGATCT.fa  fragment_ends_AGATCT
