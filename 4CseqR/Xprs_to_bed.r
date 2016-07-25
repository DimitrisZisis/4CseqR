library(stringr)

#read the data file
input_data<- read.table ("/home/dimitris/4CseqR/Data/results1_6_100length.xprs",header = TRUE)
names(input_data)
#select the columns for the analysis
split_col = as.data.frame(str_match(input_data$target_id, "^(.*):(.*)-(.*)$")[,-1])
select_data = subset(input_data,  select = c(tot_counts,uniq_counts, est_counts))
new_data = cbind(split_col,select_data)
#sort the data based on the Chromosome order 
sort_data= new_data[order(new_data$V1,sort_data$V2),]
#save file as csv/bed.Bed file format is an extension of csv but tab delimited
write.table(sort_data, file = "Xprs_reluts_INDEX_1_5.bed", sep = " ", quote = FALSE,col.names = FALSE)

#Use bedtools sort to sort the data in bed file formal.
#bedtools sort -i results1_4_100leng.bed > results1_4_100leng_sorted.bed

