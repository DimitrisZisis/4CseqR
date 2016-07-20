#Description of nescesarry libraries and packages: 
#library("Basic4Cseq", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.2")
#library(BSgenome)
#library(ShortRead)
#library(Biostrings)

#CREATE LIBRARIES OF FRAGMENTS AND FRAGMENT ENDS

# Display seed file for Athaliana TAIR10:

seed_files <- system.file("extdata", "GentlemanLab", package="BSgenome")
tail(list.files(seed_files, pattern="-seed$"))
Athaliana_seed <- list.files(seed_files, pattern="BSgenome.Athaliana.TAIR.TAIR9-seed$", full.names=TRUE)
cat(readLines(Athaliana_seed), sep="\n")

libraryName = "fragments_Athalianatest.csv"

#Create the fragments for all chromosomes. 
require(BSgenome.Athaliana.TAIR.TAIR9)
createVirtualFragmentLibrary(chosenGenome = Athaliana, firstCutter = "AGATCT", secondCutter = "CATG", readLength = 50, onlyNonBlind = FALSE, chromosomeName = Chr, libraryName = libraryName)
  
source("4CseqR_functions.R")

#keep only fragments >100 >50 etc
fragmentLenght = 100
newfragments=data.filtered(libraryName, fragmentLenght)
filteredLibraryName = "fragments_Athalianaall100test.csv"
write.table(newfragments, file = filteredLibraryName,quote = FALSE,row.names=FALSE , sep="\t")
#Create fragment ends and combine them in order to make bed files 
FragmentEnds=fragment.ends(filteredLibraryName, 7, 43, 6, 44)
write.table(FragmentEnds, file = "fragment_endsAllChrAGATCT100testt.bed",quote = FALSE,row.names=FALSE , col.names=FALSE)

#====================================================================================================================================#

#TREAMING THE DATA 
#Filter and tream the fastq files with reads in order to have the perfect.fq. Provide 2 bait files
#Create and read 2 txt file with the primer sequences and the restriction enzymes of each experiment 

# Provide the following parameters:
primer_filename = "Downloads/NewIrisData_08_2014/Clean/sample1/bait_seqTEST.txt"
restriction_filename = "Downloads/NewIrisData_08_2014/Clean/sample1/bait_seq2.txt"
fastq_filename = "Downloads/NewIrisData_08_2014/Clean/sample1/PilotExperimentJUL2015/INDEX_1_5.fq"
pattern= "AGATCT"
# Run the analysis of treaming the fastq files
PrimerSeq <- readLines(primer_filename)
RestrEnzyme <- readLines(restriction_filename)
SelectReads=select.reads(fastq_filename, PrimerSeq[1], RestrEnzyme[1], pattern)
writeFastq(SelectReads, paste("INDEX_1_5", sep="_", "treatment_results.fq"))

#select reads with the specific lenght(31)and write them in a file. 
select_reads <- grep("^.{31}$", sread(SelectReads))
new_reads <- SelectReads[select_reads]
writeFastq(new_reads, paste("INDEX_1_5", "_Perfect31.fq", sep=""))

#Extracting fasta from bed file. if you have bedtools installed (Mac/Unix) and it is in your path (export your path if not), 
#you can run the function in R to replicate the commandline task.

bedfile = "Downloads/NewIrisData_08_2014/Clean/new_fragment_ends/AGATCT_CATGfragmentends/fragment_endsAllChrAGATCT100length.bed"
reference = "Downloads/NewIrisData_08_2014/Clean/new_fragment_ends/AGATCT_CATGfragmentends/arabidopsis.fas"
Bedtools=bedTools(fstring="fastaFromBed", bedfile ,reference)

#Run bowtie2 ,run Express, Run Python Script  and the result of this is the input to normalization script 

#Normalization by Ranks 
#normalization for the 3 prime ends 
input_filename = "Downloads/NewIrisData_08_2014/Clean/new_fragment_ends/AGATCT_CATGfragmentends/finaltable100frag1_1new.csv"
categories_left = c("50", "100", "150", "200", "400","600", "800", ">800")
categories_right = c("50", "100", "150", "200","250", "300", "350", "400",">400")
categories_lenght = c("50", "100", "150", "200", "400","600", "800", ">800")
NormalizationRanks3= normalization_ranks_3prime_end(input_filename, categories_left, categories_right, categories_lenght)
#normalization for 5 prime ends 
NormalizationRanks5= normalization_ranks_5prime_end(input_filename, categories_left, categories_right, categories_lenght)
write.table(NormalizationRanks3, file = "rank3.txt", append=TRUE , quote = FALSE,  sep = "\t", row.names=FALSE , col.names=FALSE)
write.table(NormalizationRanks5, file = "rank5.txt", append=TRUE , quote = FALSE,  sep = "\t", row.names=FALSE , col.names=FALSE)
#sort the normalization rank results by chromosome and calculate the total normalization for both 3 and 5 prime ends 
rank_data5 <- read.table(file= "/home/dimitris/rank5.txt", header=FALSE, sep="\t")
rank_data3 <- read.table(file= "/home/dimitris/rank5.txt", header=FALSE, sep="\t")
sort_data5=rank_data5[with(rank_data5, order(rank_data5$V1, rank_data5$V2)), ]
sort_data3=rank_data3[with(rank_data3, order(rank_data3$V1, rank_data3$V2)), ]
totnormal = sort_data5$V5 + sort_data3$V5
total_normalized_coverage=cbind(as.character(sort_data3$V1),totnormal)
write.table(total_normalized_coverage, file = "total_normalized_coverage.txt", append=TRUE , quote = FALSE,  sep = "\t", row.names=FALSE , col.names=FALSE)

#Statistical analysis (ANOVA,pairwise test etc)


