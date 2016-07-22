# 4CseqR

4CseqR is pipeline to  investigate DNA contacts across the genome in the  model plant Arabidopsis thaliana.4CseqR is a full data analysis procedure including comparison of different samples. Basic characteristics of this pipeline is that it is based on tools commonly used for NGS data analysis and the the main algorithm is developed anew in R. This pipeline is intended for processing data for different viewpoints, genotypes or conditions. More over 4CseqR can analyze experiments including replications that allow to test the significance of differences between treatments.

For reviews about 4C-seq refer to the following articles : 
<ul>
	<li> van de Werken HJ et al. (2012) 4C technology: protocols and data

analysis. Methods Enzymol. 513: 89-112.

	<li> Stadhouders R, Kolovos P, Brouwer R et al. (2013) Multiplexed chromosome conformation

capture sequencing for rapid genome-scale high-resolution detection of long-range chromatin

interactions. Nature Protocols doi:10.1038/nprot.2013.018.
</ul>

4CseqR is designed to identify through a simple way regions across the genome that interact with the chosen 4C bait. 4CseqR provides a general method for the analysis of 4C-seq data starting from preprocessing of the NGS reads and ending with statistical inference on different 4C samples obtained for different viewpoints, genotypes or conditions, possibly including replications that allow to test the significance of differences between treatments.

The required R packages to run 4CseqR are : Basic4Cseq, Bsgenome, ShortRead, Biostrings.
```
library("Basic4Cseq", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.2")
library(BSgenome)
library(ShortRead)
showMethods(readFastq)
showMethods(writeFastq)
library(Biostrings)
```
The required  tools to run 4CseqR are : bedtools, bowtie2, eXpress and an installation of Python. 

The input to 4CseqR is the output of a 4C-seq experiment in fastq format. The 4CseqR pipeline will take this fastq output will map it to the genome, will create the estimate coverage and will normalize it prividing tables at each step for further exploration. 

An example of typical use of 4CseqR pipeline can be found in :
https://github.com/DimitrisZisis/4CseqR
A library of restriction fragments and fragment end sequences is part of the standard download (refer to step 1). In the standard avaiklable download there are also example fastq files for two data sets. 

**1. Creating the library of fragments and fragment ends**

Unlike to other methods in 4CseqR there is the opportunity to use either all or reduced genome for mapping. The reduced genome which is prepared in our case is called library of fragments. Main purpose is to create an in-silico library of restriction fragments (5’ and 3’ ends) from the given genome restricted with the primary and the secondary restriction enzymes. In case of our experiment the library of fragments contains all AGATCT fragment 5' and 3' ends (of 100bp length) from the TAIR10 genome used as the reference (file frag_ends_100AGATCT.fa ).

In case of 4CseqR an R function  is produced for this reason and it is combined with the package Basic4Cseq (BSgenome package or DNAString object). For more information about the installation and use of this package visit :
https://www.bioconductor.org/packages/devel/bioc/html/Basic4Cseq.html

To create the library of fragments in 4CseqR for Arabidopsis thaliana genome we are using the function *createVirtualFragmentLibrary* from Basic4Cseq tool. Two restriction enzymes have to be specified to cut the DNA, the read length is needed to check the fragment ends of corresponding length for uniqueness. Filter options (minimum and maximum size) are provided on fragment level and on fragment end level: 
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
In case of the example experiment only fragments longer than 100 nucleotides where used for mapping. For this reason there is  the *data.filtered* function to filter fragments :
```
fragmentLenght = 100
newfragments=data.filtered(libraryName, fragmentLenght)
filteredLibraryName = "fragments_Athalianaall100test.csv"
write.table(newfragments, file = filteredLibraryName,quote = FALSE,row.names=FALSE , sep="\t")  
```
After the creation of library of fragments  4CseqR pipeline  provides another R function *fragment.ends* for the creation or the library of fragment ends. The user is able to set the limits of each end by specifying the length of each fragment end from the start and end of each fragment. In case of the example experiment fragment ends are created based to the rule of start fragment position – 7, start fragment position + 43 and end fragment position + 6 end fragment position – 44 for the 3 and 5 prime end at each time. 
```
FragmentEnds=fragment.ends(filteredLibraryName, 7, 43, 6, 44)
write.table(FragmentEnds, file = "fragment_endsname.bed",quote = FALSE,row.names=FALSE , col.names=FALSE)
```
The result table with fragment ends is saved in .bed format and is used from bedtools in order to create the .fasta file which will be input to the mapping procedure.A script (*getfasta.sh*) has been provided for this reason. 

**2. Data Preprocessing**

Unlike other 4C-seq pipelines and methods, 4CseqR takes as input fastq files directly after sequencing. The basic criteria for filtering the original reads is to keep only those reads that start directly at a restriction enzyme cutting site. The NGS data obtained for all samples (INDEX 5-6.fq) are at the beginning processed by filtering out truncated reads and the reads that did not contain the primary restriction site AGATCT (BglII). Main purpose of this preprocessing step is to find the reads from the original files which are legitimate 4C reads which means (primer + restr.site ). For filtering and trimming the primer sequences from the fastq reads the  R function *select.reads*  is provided.

The user have to create and read 2 txt files (primer_filename, restriction_filename) with the primer sequences and the enzyme recognition sequences for each experiment  and also the fastq read files (fastq_filename). 
The restiction_filename which includes the restriction enzyme is also used as pattern in the trimming step of this function for the specific experiment.
```
primer_filename = "path to the directory where primer.txt is t"
restriction_filename = "patht to the directory where bait.txt is "
fastq_filename="path to the directory where fastq file is"
pattern= "AGATCT"
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

The primer restriction site is mapped to the library of fragments which is used as reference genome in order to create library of that contains the primary sequences and use them for mapping. 
In order to do this a script (*bowt_buid.sh*) is using the tool bowtie2 and more specific the option bowtie2-build. This function builds a Bowtie index from a set of DNA sequences and outputs a set of 6 files with suffixes .1.bt2, .2.bt2, .3.bt2, .4.bt2, .rev.1.bt2, and .rev.2.bt2.These files together are all necessary to align reads to that library of fragments. The original sequence FASTA files are no longer used by Bowtie 2 once the index is built. 
Based to the preprocessing filtering process 4CseqR ensure us that the sequence reads generated from 4C-Seq experiment contain the primer sequence starting with the primer restriction enzyme and followed by the capture. Those reads were mapped to the library of all AGATCT fragment 5' and 3' ends (of 50bp length) from the TAIR10 genome used as the reference. Mapping was conducted by Bowtie2 in all positions with up to 2 mismatches allowed in the contact region only and without mismatches in the restriction sequence. The reads mapped to both fragment ends are combined for further analysis. A script (*mapping.sh*) is provided for this step. 


To check for coherence between replicates or samples, we can visualize the .bam files in a genome browser. The mapping results can be aslo processed to obtain data concerning the total number of reads mapped to all fragments (i.e., 5' mapped reads + 3' mapped reads). These data can be formatted as "wig" files for inspection in genome browser 

**4. Estimate Coverage**

4CseqR provides a script (*run_express.sh*) which is taking care about the non- uniquely mapped reads  especially in the pericentromeric region of each chromosome.  The pipeline applies the algorithm eXpress [Roberts and Pachter (2013)] for correction of estimated coverage with respect to non-uniquely mapped reads. The coverage is calculated as a value between „all” and „unique” mapping. By using this method we are able to estimate the coverage of a set of "target" reference sequences by NGS reads using a probabilistic model.

The result of express script is (*.xprs files) a table  with the 5’ end  and 3’ end mapping results based on the : 1) total number of reads 2) estimate number of reads 3) unique number of reads and the coordinated of the fragments in which those reads are found. In case of the example experiment both estimated and unique reads were used for further analysis. 

In this part it is quicker  to use some command line orders a python script and bedtools in order to create the proper input for normalization. For this reason in 4CseqR pipeline we run the following to sort intersect and merge the express result:
```
bedtools sort -i results_eXpress.bed > results_eXpress_sorted.bed
bedtools intersect -a fragments_Athaliana.bed -b results_eXpress_sorted.bed -wa -wb > intersect_results.txt
```
The above result from bedtools is the input to python script which merges and creates the input for normalization
```
python merge_4Cnew2.py iintersect_results.txt > finaltable_results.csv
```
The finaltable_results.csv is file in which there are information about :
Chromosome, start, end position , IsNonBlind,lefS, rigS, total5, unique5, estim5, , total3, unique3, estim3

Chromosome	| start	| end	| IsNonBlind|	lefS |	rigS|	total5|	unique5	|estim5|	total3|	unique3	|estim3
------------|-------|-----|-----------|------|------|-------|--------|-------|--------|--------|-------|
Chr1	|13176|	20064	|TRUE	|256	|152|	0	|0	|0	|24	|24	|24
Chr1	|20071|	28272|	TRUE|	230|	8|	1|	1|	1|	0|	0|	0|
Chr1	|28279|	29959	|TRUE|	66|	127|	6|	6|	6|	51|	51|	51|



**5. Normalization**

4CseqR pipeline is proposing a ranked based normalization method. The resulting coverage of fragments after some table modifications as are described in step 4  (.csv files) is the input to an R function *normalization_ranks* responsible for the normalization process. The mapped counts are normalized according to the procedure based on ranking within classes of fragments (blind/non blind, have short/long ends, short/long themselves). Ranks within each category are expressed in the scale from 0 to 1 and the total coverage of 5’ and 3’ ends is expressed as values from 0 to 2. The input of the R function is the input_filename and the categories of distance and length (lefS_categories_name, lefS_categories_limits, rigS_categories_name, rigS_categories_limits, lenght_categories_names, lenght_categories_limits) based to the needs of each experiment. To normalize the example data:
```
input_filename ="path_to_your_finaltable_results.csv"
lefS_categories_name = c("50", "100", "150", "200", "400","600", "800", ">800")
lefS_categories_limits = c(-Inf, 51, 101, 151, 201, 251, 301, 351, 401, Inf)
rigS_categories_name = c("50", "100", "150", "200","250", "300", "350", "400",">400")
rigS_categories_limits = c(-Inf, 51, 101, 151, 201, 251, 301, 351, 401, Inf)
lenght_categories_names = c("50", "100", "150", "200", "400","600", "800", ">800")
lenght_categories_limits =  c(-Inf, 51, 101, 151, 201, 401, 601, 801,  Inf)

NormalizationRanks= normalization_ranks(input_filename)
```






