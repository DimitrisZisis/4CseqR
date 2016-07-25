# 4CSeqR

4CSeqR is a pipeline to investigate DNA contacts across the genome. A basic characteristics of this pipeline is that it is based on tools commonly used for NGS data analysis and the main algorithm developed anew in R. 
For reviews about the 4C-seq methodology refer to the following articles : 
<ul>
	<li> van de Werken HJ et al. (2012) 4C technology: protocols and data

analysis. Methods Enzymol. 513: 89-112.

	<li> Stadhouders R, Kolovos P, Brouwer R et al. (2013) Multiplexed chromosome conformation

capture sequencing for rapid genome-scale high-resolution detection of long-range chromatin

interactions. Nature Protocols  doi: 10.1038/nprot.2013.018
</ul>


The required R packages to run 4CSeqR are: Basic4Cseq, Bsgenome, ShortRead, Biostrings.
```
library(Basic4Cseq")
library(BSgenome)
library(ShortRead)
showMethods(readFastq)
showMethods(writeFastq)
library(Biostrings)
```
The required  tools to run 4CseqR are: bedtools, bowtie2, eXpress and an installation of Python. 

The input to 4CseqR is the data from a 4C-seq experiment in fastq format. The 4CseqR pipeline will take this fastq files, map reads to the genome, estimate coverage and normalize it providing tables at each step for further exploration. 

An example of typical use of 4CseqR pipeline can be found in :
A library of restriction fragments and fragment end sequences is part of the standard download (refer to step 1). In the standard available download there are also example fastq files for two data sets. 

**1. Creating the library of fragments and fragment ends**


To create the library of fragments for Arabidopsis thaliana genome, in 4CseqR we are using the function *createVirtualFragmentLibrary* from Basic4Cseq tool. Two restriction enzymes have to be specified that cut the DNA, the read length is needed to check the fragment ends of corresponding length for uniqueness. Filtering options (minimum and maximum size) are provided on fragment level and on fragment end level: 
```
# Display seed file for A. thaliana TAIR10:
seed_files <- system.file("extdata", "GentlemanLab", package="BSgenome")
tail(list.files(seed_files, pattern="-seed$"))
Athaliana_seed<-list.files(seed_files,pattern="Bsgenome.Athaliana.TAIR.TAIR9-seed$", full.names=TRUE)
cat(readLines(Athaliana_seed), sep="\n")
libraryName = "fragments_Athalianatest.csv"
#Create the fragments for all chromosomes 
require(BSgenome.Athaliana.TAIR.TAIR9)
createVirtualFragmentLibrary(chosenGenome = Athaliana, firstCutter = "AGATCT", secondCutter = "CATG", readLength = 50, onlyNonBlind = FALSE, chromosomeName = Chr, libraryName = libraryName)
```
In case of the example experiment only fragments longer than 100 nucleotides were used for mapping. For this reason there is  the *data.filtered* function to filter fragments:
```
fragmentLenght = 100
newfragments=data.filtered(libraryName, fragmentLenght)
filteredLibraryName = "fragments_Athalianaall100test.csv"
write.table(newfragments, file = filteredLibraryName,quote = FALSE,row.names=FALSE , sep="\t")  
```
After the creation of library of fragments 4CseqR pipeline provides another R function *fragment.ends* for the creation or the library of fragment ends. By analyzing the length (in this experiment is 50) and the sequence of the restriction enzyme (in this experiment is AGATCT) we calculate fragment ends around the start and end position of each fragment as bellow:

```

fragment__start_min <- filteredLibraryName[,2]-7
fragment_start_plus<- filteredLibraryName[,2]+43
fragment_end_plus<- filteredLibraryName[,3]+6
fragment_end_min <- filteredLibraryName[,3]-44
FragmentEnds=fragment.ends(filteredLibraryName,fragment__start_min,fragment_start_plus,fragment_end_plus,fragment_end_min)

```
The resulting table with fragment ends is saved in .bed format and is used in bedtools to create the fasta file which will be input to the mapping procedure. A script (*getfasta.sh*) has been provided for this reason. 

**2. Data Preprocessing**

4CseqR takes as input fastq files directly after sequencing. It filters the original reads to keep only those that start directly at a restriction enzyme cutting site. The NGS data obtained for all samples (INDEX1_5.fq, index1_6.fq) are at the beginning processed by filtering out truncated reads and the reads that did not contain the primary restriction site AGATCT (BglII). Main purpose of this preprocessing step is to find the reads from the original files which are legitimate 4C reads which means that they contain the sequence (primer + restr.site ). For filtering and trimming the primer sequences from the fastq reads the R function *select.reads* is provided.

The user have to create two text files (primer_filename, restriction_filename) with the primer sequences and the enzyme recognition sequences for each experiment  and also the fastq read files (fastq_filename). 
The restiction_filename which includes the restriction enzyme is also used as pattern in the trimming step of this function for the specific experiment.
```
primer_filename = "path to the directory where primer.txt is t"
restriction_filename = "patht to the directory where bait.txt is "
fastq_filename="path to the directory where fastq file is"
PrimerSeq <- readLines(primer_filename)
RestrEnzyme <- readLines(restriction_filename)
SelectReads=select.reads(fastq_filename, PrimerSeq[1], RestrEnzyme[1], pattern)
writeFastq(SelectReads, paste("INDEX_1_5", sep="_", "treatment_results.fq"))
```

Another additional process which can be done in preprocessing is the filtering of short reads in the output file which can be done by grep 
```
select_reads <- grep("^.{31}$", sread(SelectReads))
SelectReads <- SelectReads[select_reads]
writeFastq(SelectReads, paste("INDEX_1_5", "_Perfect31.fq", sep=""))
```

**3. Mapping of primer sequence and Mapping of reads**

The reads were mapped to the library of all AGATCT fragment 5' and 3' ends (of 50 bp length) from the TAIR10 genome used as the reference. Mapping was conducted by Bowtie2 in all positions with up to 2 mismatches allowed in the contact region only and without mismatches in the restriction sequence. The reads mapped to both fragment ends are combined for further analysis. A script (*mapping.sh*) is provided for this step. 



**4. Estimate Coverage**

4CseqR provides a script (*run_express.sh*) which is taking care about the non-uniquely mapped reads  especially in the pericentromeric region of each chromosome.  The pipeline applies the algorithm eXpress [Roberts and Pachter (2013)] for correction of estimated coverage with respect to non-uniquely mapped reads. The coverage is calculated as a value between „all” and „unique” mapping. By using this method we are able to estimate the coverage of a set of "target" reference sequences by NGS reads using a probabilistic model.

The result of express script is (*.xprs files) a table  with the 5’ end  and 3’ end mapping results based on the : 1) total number of reads 2) estimated number of reads 3) unique number of reads, and the coordinates of the fragments in which those reads are found. 
Command line,bedtools and a python script are used in order to  create the proper input for normalization. An R script (*Xprs_to_bed.r*) is provided to create bed files and the Terminal_commads_4CseqR file describes the steps which are necessary in order to have the proper .csv input in the normalization function. 

The finaltable1_5leng100.csv is the file in which there are informations about:
Chromosome ID, start and end position of the fragment, IsNonBlind, left and right length (lefS, rigS) and also the total, unique and estimate coverage as it is calculated from eXpress( total5, unique5, estim5, , total3, unique3, estim3)

Chromosome	| start	| end	| IsNonBlind|	lefS |	rigS|	total5|	unique5	|estim5|	total3|	unique3	|estim3
------------|-------|-----|-----------|------|------|-------|--------|-------|--------|--------|-------|
Chr1	|13176|	20064	|TRUE	|256	|152|	0	|0	|0	|24	|24	|24
Chr1	|20071|	28272|	TRUE|	230|	8|	1|	1|	1|	0|	0|	0|
Chr1	|28279|	29959	|TRUE|	66|	127|	6|	6|	6|	51|	51|	51|



**5. Normalization**

The counts of reads are normalized according to the procedure based on ranking within classes of fragments (blind/non blind, length of ends, total length). Ranks within each category are expressed in the scale from 0 to 1 and the total coverage of 5’ and 3’ ends is expressed as a value from 0 to 2. The input of the R function is the input_filename and the categories of distance and length (lefS_categories_name, lefS_categories_limits, rigS_categories_name, rigS_categories_limits, lenght_categories_names, lenght_categories_limits) based to the needs of each experiment. To normalize the example data:
```
input_filename = "/home/dimitris/4CseqR/Data/finaltable1_5leng100.csv"
#Create categories of distance and length (edit the limits to change intervals)
lefS_categories_limits = c(50, 100, 150, 200, 250, 300, 350, 400)
rigS_categories_limits = c(50, 100, 150, 200, 250, 300, 350, 400)
lenght_categories_limits = c(50, 100, 150, 200, 400, 600, 800)
NormalizationRanks= normalization_ranks(input_filename,lefS_categories_limits, rigS_categories_limits, lenght_categories_limits)```





