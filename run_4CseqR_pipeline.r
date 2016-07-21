#Description of nescesarry libraries and packages: 
#library("Basic4Cseq", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.2")
#library(BSgenome)
#library(ShortRead)
#library(Biostrings)

#CREATE LIBRARIES OF FRAGMENTS AND FRAGMENT ENDS

# Display seed file for A. thaliana TAIR10:

seed_files <- system.file("extdata", "GentlemanLab", package="BSgenome")
tail(list.files(seed_files, pattern="-seed$"))
Athaliana_seed <- list.files(seed_files, pattern="BSgenome.Athaliana.TAIR.TAIR9-seed$", full.names=TRUE)
cat(readLines(Athaliana_seed), sep="\n")

libraryName = "fragments_Athalianatest.csv"

#Create the fragments for all chromosomes 
require(BSgenome.Athaliana.TAIR.TAIR9)
createVirtualFragmentLibrary(chosenGenome = Athaliana, firstCutter = "AGATCT", secondCutter = "CATG", readLength = 50, onlyNonBlind = FALSE, chromosomeName = Chr, libraryName = libraryName)
  
source("4CseqR_functions.R")

#keep only fragments longer than 100 nucleotides
fragmentLenght = 100
newfragments=data.filtered(libraryName, fragmentLenght)
filteredLibraryName = "fragments_Athalianaall100test.csv"
write.table(newfragments, file = filteredLibraryName,quote = FALSE,row.names=FALSE , sep="\t")
#Create fragment ends and combine them in order to make bed files 
FragmentEnds=fragment.ends(filteredLibraryName, 7, 43, 6, 44)
write.table(FragmentEnds, file = "fragment_endsAllChrAGATCT100testt.bed",quote = FALSE,row.names=FALSE , col.names=FALSE)

#====================================================================================================================================#

#TRIMING THE DATA 
#Filter and trim the fastq files with reads in order to have the file of "restriction site + contact" sequences
#Create and read 2 txt files with the primer sequences and the enzyme recognition sequences for each experiment 

# Provide the following parameters:
primer_filename = "Downloads/NewIrisData_08_2014/Clean/sample1/bait_seqTEST.txt"
restriction_filename = "Downloads/NewIrisData_08_2014/Clean/sample1/bait_seq2.txt"
fastq_filename = "Downloads/NewIrisData_08_2014/Clean/sample1/PilotExperimentJUL2015/INDEX_1_5.fq"
!!!!!!!!!!!!!!!! pattern= "AGATCT"
# Run the triming the reads in fastq files
PrimerSeq <- readLines(primer_filename)
RestrEnzyme <- readLines(restriction_filename)
SelectReads=select.reads(fastq_filename, PrimerSeq[1], RestrEnzyme[1], pattern)
writeFastq(SelectReads, paste("INDEX_1_5", sep="_", "treatment_results.fq"))

#select reads with the specific length(31)and write them in a file. 
select_reads <- grep("^.{31}$", sread(SelectReads))
SelectReads <- SelectReads[select_reads]
writeFastq(SelectReads, paste("INDEX_1_5", "_Perfect31.fq", sep=""))

#Run bowtie2, run Express, run Python script to extract results to make input for normalization script 

#Normalization by Ranks 

#normalization for the 3 and 5 prime ends 
input_filename = "Downloads/NewIrisData_08_2014/Clean/new_fragment_ends/AGATCT_CATGfragmentends/finaltable100frag1_1new.csv"
#create categories of distance and length (edit the limits to change intervals)
lefS_categories_name = c("50", "100", "150", "200", "400","600", "800", ">800")
lefS_categories_limits = c(-Inf, 51, 101, 151, 201, 251, 301, 351, 401, Inf)
rigS_categories_name = c("50", "100", "150", "200","250", "300", "350", "400",">400")
rigS_categories_limits = c(-Inf, 51, 101, 151, 201, 251, 301, 351, 401, Inf)
lenght_categories_names = c("50", "100", "150", "200", "400","600", "800", ">800")
lenght_categories_limits =  c(-Inf, 51, 101, 151, 201, 401, 601, 801,  Inf)

NormalizationRanks= normalization_ranks_prime_end(input_filename)
write.table(NormalizationRanks, file = "testrs.txt", append=TRUE , quote = FALSE,  sep = "\t", row.names=FALSE , col.names=FALSE)


#sort the normalized data by chromosome and calculate the total normalized coverage of fragments 
rank_data5 <- read.table(file= "/home/dimitris/rank5.txt", header=FALSE, sep="\t")
rank_data3 <- read.table(file= "/home/dimitris/rank5.txt", header=FALSE, sep="\t")
sort_data5=rank_data5[with(rank_data5, order(rank_data5$V1, rank_data5$V2)), ]
sort_data3=rank_data3[with(rank_data3, order(rank_data3$V1, rank_data3$V2)), ]
totnormal = sort_data5$V5 + sort_data3$V5
total_normalized_coverage=cbind(as.character(sort_data3$V1),totnormal)
write.table(total_normalized_coverage, file = "total_normalized_coverage.txt", append=TRUE , quote = FALSE,  sep = "\t", row.names=FALSE , col.names=FALSE)




