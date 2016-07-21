# Display seed file for Athaliana TAIR10:
library("Basic4Cseq", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.2")
library(BSgenome)
library(ShortRead)
showMethods(readFastq)
showMethods(writeFastq)
library(Biostrings)
#Filter data bigger than 100
data.filtered <- function(filename,n) {
  Fragments <- read.table (filename,header = TRUE)
  FilterFragments<-subset(Fragments, fragmentLength>=n )
  return(FilterFragments)
}
#save("data.filtered", file="FilterFragmentsFunct.Rdata")

#Create fragment ends and combine them in order to make bed files 

fragment.ends <- function (filename2, start_min,start_plus, end_plus, end_min) {
  #read the file
  MyData <- read.table (filename2,header = TRUE)
  #create the fragment ends
  a <- MyData[,2]-start_min
  b <- MyData[,2]+start_plus
  d <- MyData[,3]+end_plus
  c <- MyData[,3]-end_min
  #combine the fragments ends 
  x <- cbind(a,b,c,d)
  vec <- matrix(t(x),nrow=1)
  final <- t(matrix((vec),nrow=2))
  names2 <- cbind(levels(MyData$chromosomeName)[MyData$chromosomeName], levels(MyData$chromosomeName)[MyData$chromosomeName])
  names1 <- matrix(t(names2))
  final_correct <- cbind(names1, final)
  return(final_correct)

}
#Filter and tream the fastq files with reads in order to have the perfect.fq. Provide 2 bait files
select.reads  <- function (filename,PrimerSeq,RestrEnzyme,pattern){
  # create and read the data file 
  fq <- readFastq(dirPath=filename)
  reads <- sread(fq)
  quality <- quality(fq)
  #count ocurrences of bait
  pos <- grepl(PrimerSeq, reads)
  # fuzzy match - 2 substitutions
  fpos1 <- agrep(RestrEnzyme, reads, max = list(sub =0, ins=0, del=0))
  # use writeFastq functions in ShortRead package to save "good" reads from fuzzy match in "*.good.fq" file and "bad" reads in *".bad.fq" file
  good.reads <- reads[fpos1]
  good.quality <- quality[fpos1]
  good.result <- ShortReadQ(good.reads, good.quality)
  #Create a new file to save only reads without primer(PrimerSeq) in the *.perfect.fq file 
  match <- regexpr(pattern=pattern, text=good.reads)
  if (match >= 0) {
    perfect.reads <- substr(good.reads, match, width(good.reads))
    perfect.quality <- substr(good.quality@quality, match, width(good.quality@quality))
     
  }
  final.reads <- DNAStringSet(perfect.reads)
  perfect.quality.final <- SFastqQuality(quality=BStringSet(perfect.quality))
  perfect.result <- ShortReadQ(final.reads, perfect.quality.final)

  return(perfect.result)
}

#Running bedtools and their functions in R
#fasta from bed 
bedTools<-function(fstring="fastaFromBed", bedframe, fasta_file, opt="-s -tab"){
  #create temp files
  bed.file= bedfile = "Downloads/NewIrisData_08_2014/Clean/new_fragment_ends/AGATCT_CATGfragmentends/fragment_endsAllChrAGATCT100length.bed"
  out   = tempfile()
  
  # writing temporary files
  write.table(bedframe,file=bed.file,quote=F,sep="\t",col.names=F,row.names=F)
  
  # call the function
  command=paste(fstring,opt,"-fi",fasta_file,"-bed",bed.file,"-fo",out,sep=" ")
  cat(command,"\n")
  try(system(command))
  
  res=read.table(out,header=F,stringsAsFactors =FALSE)
  unlink(bedfile);unlink(out)
  return(res)
}

#normalization by Ranks for the 3 prime ends 
normalization_ranks_3prime_end<- function (input_filename, categories_left, categories_right, categories_lenght ){
  MyData <- read.table (input_filename,header = TRUE)
  #check if any row of each column is 0
  any(MyData$rigS == 0)
  #calculate the lenght and the middle point
  lenght = MyData$end - MyData$start
  middle = (MyData$end + MyData$start)/2
  
  #replace 0 with 1 if there are
  MyData$IsNonBlind <- as.character(MyData$IsNonBlind)
  MyData$IsNonBlind[MyData$IsNonBlind %in% c("TRUE")] <- 0
  MyData$IsNonBlind[MyData$IsNonBlind %in% c("FALSE")] <- 1
  fblind <- MyData$IsNonBlind
  #calculate left and right distance
  ldist = MyData$lefS
  rdist =MyData$rigS
  #create categories of distance and lenght 
  fldist<-categories_left[
    findInterval(ldist , c(-Inf, 51, 101, 151, 201, 251, 301, 351, 401, Inf) ) ]
  frdist <- categories_right [
    findInterval(rdist , c(-Inf, 51, 101, 151, 201, 251, 301, 351, 401, Inf) ) ]
  flenght <- categories_lenght[
    findInterval(lenght , c(-Inf, 51, 101, 151, 201, 401, 601, 801,  Inf) ) ]
  
  kat3 <- cbind(fblind,frdist,flenght)
  #Normalization for 3 prime ends
  V3=paste(kat3[,1],kat3[,2],kat3[,3])
  at3=sort(unique(V3))
  nlev=length(at3)
  ach3=NULL
  arst3=NULL
  aren3=NULL
  ac3=NULL
  ac3h=NULL
  rank3=NULL

  for (kk in 1:nlev){
    ach3[[kk]]=MyData$Chromosome[which(V3==at3[kk])]
    arst3[[kk]]=MyData$start[which(V3==at3[kk])]
    aren3[[kk]]=MyData$end[which(V3==at3[kk])]
    ac3[[kk]]=MyData$estim3[which(V3==at3[kk])]
    v3=ac3[[kk]]
    v3[which(v3==0)]=NA
    rank3[[kk]]=rank(v3,na.last=FALSE)
    rank3[[kk]]=rank3[[kk]]/length(v3)
    rank3[[kk]][is.na(v3)] <- 0
    #U3=cbind(as.character(ach3[[kk]]),arst3[[kk]],aren3[[kk]],ac3[[kk]],rank3[[kk]],at3[[kk]])
    
  }
  total_results_3prime=NULL
  for (kk in 1:nlev) {
    rank_results_3prime=cbind(as.character(ach3[[kk]]),arst3[[kk]],aren3[[kk]],ac3[[kk]],rank3[[kk]],at3[[kk]])
    total_results_3prime <- rbind(total_results_3prime, rank_results_3prime)
  }
  return(total_results_3prime)
}

#normalization by Ranks for the 5 prime ends 
normalization_ranks_5prime_end <- function (input_filename, categories_left, categories_right, categories_lenght ){
  MyData <- read.table (input_filename,header = TRUE)
  #check if any row of each column is 0
  any(MyData$rigS == 0)
  #calculate the lenght and the middle point
  lenght = MyData$end - MyData$start
  middle = (MyData$end + MyData$start)/2
  
  #replace 0 with 1 if there are 
  MyData$IsNonBlind <- as.character(MyData$IsNonBlind)
  MyData$IsNonBlind[MyData$IsNonBlind %in% c("TRUE")] <- 0
  MyData$IsNonBlind[MyData$IsNonBlind %in% c("FALSE")] <- 1
  fblind <- MyData$IsNonBlind
  #calculate left and right distance
  ldist = MyData$lefS
  rdist = MyData$rigS
  #create categories of distance and lenght 
  fldist <- categories_left [
    findInterval(ldist , c(-Inf, 51, 101, 151, 201, 251, 301, 351, 401, Inf) ) ]
  frdist <- categories_right [
    findInterval(rdist , c(-Inf, 51, 101, 151, 201, 251, 301, 351, 401, Inf) ) ]
  flenght <- categories_lenght [
    findInterval(lenght , c(-Inf, 51, 101, 151, 201, 401, 601, 801,  Inf) ) ]
  
  kat5 <- cbind(fblind,frdist,flenght)
  #Normalization for 5 prime ends
  V5=paste(kat5[,1],kat5[,2],kat5[,3])
  at5=sort(unique(V5))
  nlev=length(at5)
  ach5=NULL
  arst5=NULL
  aren5=NULL
  ac5=NULL
  ac5h=NULL
  rank5=NULL
  
  for (kk in 1:nlev){
    ach5[[kk]]=MyData$Chromosome[which(V5==at5[kk])]
    arst5[[kk]]=MyData$start[which(V5==at5[kk])]
    aren5[[kk]]=MyData$end[which(V5==at5[kk])]
    ac5[[kk]]= MyData$estim5[which(V5==at5[kk])]
    v5=ac5[[kk]]
    v5[which(v5==0)]=NA
    rank5[[kk]]=rank(v5,na.last=FALSE)
    rank5[[kk]]=rank5[[kk]]/length(v5)
    rank5[[kk]][is.na(v5)] <- 0
    #U3=cbind(as.character(ach3[[kk]]),arst3[[kk]],aren3[[kk]],ac3[[kk]],rank3[[kk]],at3[[kk]])
    
  }
  total_results_5prime = NULL
  for (kk in 1:nlev) {
    rank_results_5prime=cbind(as.character(ach5[[kk]]),arst5[[kk]],aren5[[kk]],ac5[[kk]],rank5[[kk]],at5[[kk]])
    total_results_5prime <- rbind(total_results_5prime, rank_results_5prime)
  }
  return(total_results_5prime)
}
