#!/bin/bash
#run bedtools function
BT=/pathto the bed tools installation /fastaFromBed
DATA_PATH="path to your data"
cd $DATA_PATH
$BT -fi arabidopsis.fas -bed fragment_endsAllChrAGATCT100.bed -fo fragment_endsAllChrAGATCT100.fa

exit



