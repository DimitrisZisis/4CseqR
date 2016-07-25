#!/bin/bash

DATA_PATH="/path/to/your/data"
cd $DATA_PATH

#Alignment with Express
#chmod +x express

./express --no-bias-correct -o xprs_INDEX_1_5_100length fragment_endsAllChrAGATCT100.fa INDEX_1_5_treatment_2mis.sam
./express --no-bias-correct -o xprs_INDEX_1_6_100length fragment_endsAllChrAGATCT100.fa INDEX_1_6_treatment_2mis.sam
exit

