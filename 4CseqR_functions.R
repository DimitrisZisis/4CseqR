# Display seed file for A. thaliana TAIR10:
library("Basic4Cseq", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.2")
library(BSgenome)
library(ShortRead)
showMethods(readFastq)
showMethods(writeFastq)
library(Biostrings)

#Filter out fragments shorter than n nucleotides(e.g. 100 nucleotides)
data.filtered <- function(filename,n) {
  Fragments <- read.table (filename,header = TRUE)
  FilterFragments<-subset(Fragments, fragmentLength>=n )
  return(FilterFragments)
}

#Create fragment ends and combine them in order to make bed files 
fragment.ends <- function (filename2, start_min, start_plus, end_plus, end_min) {
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
#Filter and trim the fastq files with reads in order to have the file of "restriction site + contact" sequences. 
select.reads  <- function (filename,PrimerSeq,RestrEnzyme,pattern){
  # create and read the data file 
  fq <- readFastq(dirPath=filename)
  reads <- sread(fq)
  quality <- quality(fq)
  #count ocurrences of bait
  pos <- grepl(PrimerSeq, reads)
  # fuzzy search - 2 mismatches
  fpos1 <- agrep(RestrEnzyme, reads, max = list(sub =0, ins=0, del=0))
  # use writeFastq functions in ShortRead package to separate "good" reads from "bad" reads
  good.reads <- reads[fpos1]
  good.quality <- quality[fpos1]
  good.result <- ShortReadQ(good.reads, good.quality)
  #Create a new file to save only reads without primer(PrimerSeq)
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

#normalization by Ranks for the 3 and 5 prime ends
normalization_ranks<- function (input_filename,lefS_categories_limits,rigS_categories_limits, lenght_categories_limits){
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
  #create categories of distance and length (edit the limits to change intervals)
  l = factor(lefS_categories_limits)
  maxlabel = paste(">", max(lefS_categories_limits), sep="")
  l[length(l) + 1] <- levels(l)[length(l)+1] <- maxlabel
  
  r = factor(rigS_categories_limits)
  maxlabel = paste(">", max(rigS_categories_limits), sep="")
  r[length(r) + 1] <- levels(r)[length(r)+1] <- maxlabel
  
  len = factor(lenght_categories_limits)
  maxlabel = paste(">", max(lenght_categories_limits), sep="")
  len[length(len) + 1] <- levels(len)[length(len)+1] <- maxlabel
  
  fldist = l[findInterval(ldist , lefS_categories_limits+1 )+1]
  frdist = r[findInterval(rdist , rigS_categories_limits+1 )+1]
  flenght = len [findInterval(lenght , lenght_categories_limits+1 )+1 ]
  
  kat3 <- cbind(fblind,frdist,flenght)
  kat5 <- cbind(fblind,frdist,flenght)
  #Normalization parameters for 3 prime ends
  V3=paste(kat3[,1],kat3[,2],kat3[,3])
  at3=sort(unique(V3))
  nlev3=length(at3)
  ach3=NULL
  arst3=NULL
  aren3=NULL
  ac3=NULL
  ac3h=NULL
  rank3=NULL
  #Normalization parameters for 5 prime ends
  V5=paste(kat5[,1],kat5[,2],kat5[,3])
  at5=sort(unique(V5))
  nlev5=length(at5)
  ach5=NULL
  arst5=NULL
  aren5=NULL
  ac5=NULL
  ac5h=NULL
  rank5=NULL
  
  for (kk in 1:nlev3){
    ach3[[kk]]=MyData$Chromosome[which(V3==at3[kk])]
    arst3[[kk]]=MyData$start[which(V3==at3[kk])]
    aren3[[kk]]=MyData$end[which(V3==at3[kk])]
    ac3[[kk]]=MyData$estim3[which(V3==at3[kk])]
    v3=ac3[[kk]]
    v3[which(v3==0)]=NA
    rank3[[kk]]=rank(v3,na.last=FALSE)
    rank3[[kk]]=rank3[[kk]]/length(v3)
    rank3[[kk]][is.na(v3)] <- 0
    
  }
  total_results_3prime=NULL
  for (kk in 1:nlev3) {
    rank_results_3prime=cbind(as.character(ach3[[kk]]),arst3[[kk]],aren3[[kk]],ac3[[kk]],rank3[[kk]],at3[[kk]])
    total_results_3prime <- rbind(total_results_3prime, rank_results_3prime)
  }
  for (kk in 1:nlev5){
    ach5[[kk]]=MyData$Chromosome[which(V5==at5[kk])]
    arst5[[kk]]=MyData$start[which(V5==at5[kk])]
    aren5[[kk]]=MyData$end[which(V5==at5[kk])]
    ac5[[kk]]= MyData$estim5[which(V5==at5[kk])]
    v5=ac5[[kk]]
    v5[which(v5==0)]=NA
    rank5[[kk]]=rank(v5,na.last=FALSE)
    rank5[[kk]]=rank5[[kk]]/length(v5)
    rank5[[kk]][is.na(v5)] <- 0
    
  }
  total_results_5prime = NULL
  for (kk in 1:nlev5) {
    rank_results_5prime=cbind(as.character(ach5[[kk]]),arst5[[kk]],aren5[[kk]],ac5[[kk]],rank5[[kk]],at5[[kk]])
    total_results_5prime <- rbind(total_results_5prime, rank_results_5prime)
  }
  #Choose if you want to return a list with  the normalized results for the 3  the 5 or both prime ends 
  result = list(total_results_3prime, total_results_5prime)
  return(list(result))

}


