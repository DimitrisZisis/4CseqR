#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l mem=4GB
#PBS -l walltime=00:10:00
#PBS -N 4CSEQfa

BT=/home/users/dzisis/reef/bedtools-2.17.0/bin/fastaFromBed


DATA_PATH="/home/users/dzisis/reef/IrisHovel/NewIrisData/Arabidopsis/Samples1-4"
cd $DATA_PATH

check()
{
	if diff $1 $2; then
    	echo ok
	else
    	echo fail
	fi
}

$BT -fi arabidopsis.fas -bed fragment_endsAllChrAGATCT100length.bed -fo fragment_endsAllChr_AGATCT100length.fa

exit



