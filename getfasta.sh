#!/bin/bash
#run bedtools function
BT=/pathto the bed tools installation /fastaFromBed
DATA_PATH="path to your data"
cd $DATA_PATH

check()
{
	if diff $1 $2; then
    	echo ok
	else
    	echo fail
	fi
}

$BT -fi arabidopsis.fas -bed fragment_endsAllChrAGATCT100.bed -fo fragment_endsAllChrAGATCT100.fa

exit



