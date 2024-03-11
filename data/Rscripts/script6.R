#Rscript6 - Filter


suppressPackageStartupMessages(require('stringr'))

settings <- read.delim("/settings.sh", header = F, stringsAsFactors = F, comment.char = "#")
folder_a <- substring(settings[1,1], 10, )
folder_b <- substring(settings[2,1], 10, )


#set A
files_a<-list.files(paste0(folder_a, "/counting/"), recursive=T, pattern=".count", full.names = TRUE)
files_a<-files_a[grep("filter_", files_a, invert = FALSE)]
files_a<-files_a[grep("summary", files_a, invert = TRUE)]
length(files_a)
count_data_a<-data.frame(matrix())
for (n in 1:length(files_a)) {
  count_a<-read.delim(files_a[n], header=TRUE, sep= "\t", stringsAsFactors = FALSE, comment.char= "#")
  count_a<-count_a[,c(1,7)]
  count_data_a<-cbind(count_data_a,count_a)
}
rownames(count_data_a)<-count_data_a[,2]
count_data_a<-count_data_a[,seq(from = 3, to = ncol(count_data_a), by = 2)]

write.csv(count_data_a,file=paste0(folder_a, "/R_outputs/set_a_counting.csv"),quote = FALSE, row.names = TRUE)


file_a<-read.csv(paste0(folder_a, "/R_outputs/set_a_counting.csv"), sep = ",", header = T)
rownames(file_a)<-file_a$X
file_a<-file_a[,2:ncol(file_a)]
file_a$count<-0
file_a[,1:ncol(file_a)] = lapply(file_a[,1:ncol(file_a)], FUN = function(y){as.numeric(y)})

for (i in 1:nrow(file_a)) {
  for (j in 1:ncol(file_a)){
    if (file_a[i,j] >=1) {
      file_a[i,(ncol(file_a))]=file_a[i,(ncol(file_a))]+1
    }
  }
}

file_a$maxcount<- rowSums(file_a[,1:ncol(file_a)-1])
transcripts_a_not_counted<-subset(file_a,file_a$count <= (0.10*ncol(file_a))) # +90% of samples without counting
transcripts_a_counted<-subset(file_a,file_a$count > (0.10*ncol(file_a)))
transcripts_a_not_counted_b<-as.data.frame(transcripts_a_not_counted[,"count"])
transcripts_a_not_counted_b$transcript<-rownames(transcripts_a_not_counted)
transcripts_a_not_counted_b<-transcripts_a_not_counted_b[,c(2,1)]
colnames(transcripts_a_not_counted_b)<- c("transcript", "count")
transcripts_a_not_counted_c<-transcripts_a_not_counted_b[,1]


gtf_a<-read.table(paste0(folder_a, "/comparisons/lncRNA/set_a_all_flags_u_clean_concatenated2.gtf"), sep = "\t")
gtf_a$gene_id<-word(gtf_a$V9, start = 2, end = 2, sep = ";")
gtf_a$gene_id<-gsub("gene_id ", "", gtf_a$gene_id)
gtf_a$gene_id<-gsub(" ", "", gtf_a$gene_id)
gtf_a$transcript<-word(gtf_a$V9, start = 1, end = 1, sep = ";")
gtf_a$transcript<-gsub("transcript_id ", "", gtf_a$transcript)
gtf_a$transcript<-gsub(" ", "", gtf_a$transcript)
gtf_a_transcripts<-gtf_a[gtf_a$V3 == "transcript",]


#Set B
files_b<-list.files(paste0(folder_b, "/counting/"), recursive=T, pattern=".count", full.names = TRUE)
files_b<-files_b[grep("filter_", files_b, invert = FALSE)]
files_b<-files_b[grep("summary", files_b, invert = TRUE)]
length(files_b)
count_data_b<-data.frame(matrix())
for (n in 1:length(files_b)) {
  count_b<-read.delim(files_b[n], header=TRUE, sep= "\t", stringsAsFactors = FALSE, comment.char= "#")
  count_b<-count_b[,c(1,7)]
  count_data_b<-cbind(count_data_b,count_b)
}
rownames(count_data_b)<-count_data_b[,2]
count_data_b<-count_data_b[,seq(from = 3, to = ncol(count_data_b), by = 2)]

write.csv(count_data_b,file=paste0(folder_b, "/R_outputs/set_b_counting.csv"),quote = FALSE, row.names = TRUE)


file_b<-read.csv(paste0(folder_b, "/R_outputs/set_b_counting.csv"), sep = ",", header = T)
rownames(file_b)<-file_b$X
file_b<-file_b[,2:ncol(file_b)]
file_b$count<-0
file_b[,1:ncol(file_b)] = lapply(file_b[,1:ncol(file_b)], FUN = function(y){as.numeric(y)})

for (i in 1:nrow(file_b)) {
  for (j in 1:ncol(file_b)){
    if (file_b[i,j] >=1) {
      file_b[i,(ncol(file_b))]=file_b[i,(ncol(file_b))]+1
    }
  }
}

file_b$maxcount<- rowSums(file_b[,1:ncol(file_b)-1])
transcripts_b_not_counted<-subset(file_b,file_b$count <= (0.10*ncol(file_b))) # +90% of samples without counting
transcripts_b_counted<-subset(file_b,file_b$count > (0.10*ncol(file_b)))
transcripts_b_not_counted_b<-as.data.frame(transcripts_b_not_counted[,"count"])
transcripts_b_not_counted_b$transcript<-rownames(transcripts_b_not_counted)
transcripts_b_not_counted_b<-transcripts_b_not_counted_b[,c(2,1)]
colnames(transcripts_b_not_counted_b)<- c("transcript", "count")
transcripts_b_not_counted_c<-transcripts_b_not_counted_b[,1]


gtf_b<-read.table(paste0(folder_b, "/comparisons/lncRNA/set_b_all_flags_u_clean_concatenated2.gtf"), sep = "\t")
gtf_b$gene_id<-word(gtf_b$V9, start = 2, end = 2, sep = ";")
gtf_b$gene_id<-gsub("gene_id ", "", gtf_b$gene_id)
gtf_b$gene_id<-gsub(" ", "", gtf_b$gene_id)
gtf_b$transcript<-word(gtf_b$V9, start = 1, end = 1, sep = ";")
gtf_b$transcript<-gsub("transcript_id ", "", gtf_b$transcript)
gtf_b$transcript<-gsub(" ", "", gtf_b$transcript)
gtf_b_transcripts<-gtf_b[gtf_b$V3 == "transcript",]

#Summary

#set A
transcripts_a_to_remove<-subset(gtf_a_transcripts, gtf_a_transcripts$gene_id %in% transcripts_a_not_counted_c)
print("Transcripts to remove from set A without counts:")
table(transcripts_a_to_remove$V3)
transcripts_a_to_remove_list<-transcripts_a_to_remove$transcript

#set B
transcripts_b_to_remove<-subset(gtf_b_transcripts, gtf_b_transcripts$gene_id %in% transcripts_b_not_counted_c)
print("Transcripts to remove from set B without counts:")
table(transcripts_b_to_remove$V3)
transcripts_b_to_remove_list<-transcripts_b_to_remove$transcript


write.table(transcripts_a_to_remove_list, file = paste0(folder_a, "/R_outputs/not_counted_transcripts_set_a.txt"), sep = "\t",append = FALSE, quote = FALSE, row.names = FALSE, col.names=FALSE)
write.table(transcripts_b_to_remove_list, file = paste0(folder_b, "/R_outputs/not_counted_transcripts_set_b.txt"), sep = "\t",append = FALSE, quote = FALSE, row.names = FALSE, col.names=FALSE)


save.image(file=paste0(folder_a, "/R_outputs/Rscript_6.RData"))
save.image(file=paste0(folder_b, "/R_outputs/Rscript_6.RData"))
