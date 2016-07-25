#!/bin/bash
DATA_PATH="/path/to/your/data"
cd $DATA_PATH

#Alignment with Bowtie 2
#chmod +x bowtie2-align
#max 2 mismatches, min score = .40*31 = 12.4
./bowtie2-align -q -x fragment_endsAllChrAGATCT100 -a --score-min L,0,-0.4 -U INDEX_1_5_treatment_results31.fq -S INDEX_1_5_treatment_2mis.sam
./bowtie2-align -q -x fragment_endsAllChrAGATCT100 -a --score-min L,0,-0.4 -U INDEX_1_6_treatment_results29.fq -S INDEX_1_6_treatment_2mis.sam

