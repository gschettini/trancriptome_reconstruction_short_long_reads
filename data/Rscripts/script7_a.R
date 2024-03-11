#Rscript 7


suppressPackageStartupMessages(require('stringr'))

settings <- read.delim("/settings.sh", header = F, stringsAsFactors = F, comment.char = "#")
folder_a <- substring(settings[1,1], 10, )

#set A
files_a<-list.files(paste0(folder_a,"/functional_characterization"), recursive=F, pattern="_match.tsv", full.names = TRUE)
files_a<-files_a[grep("_no_match", files_a, invert = TRUE)]
info<-file.info(files_a)
empty<-rownames(info[info$size == 0L, ])
files_a<-files_a[!(files_a %in% empty)]


organism_a<-data.frame(data.frame())
for (i in 1:length(files_a))  {
  organism<-read.delim(files_a[i], header=FALSE, sep= "\t", stringsAsFactors = FALSE)
  organism<-organism[,c(1,2,15)]
  organism$V3<-files_a[i]
  organism$V3<-word(organism$V3, start = 5, end = 5, sep = "/")
  organism$V3<-gsub("set_a_nr_", "organism_", organism$V3)
  organism$V3<-gsub("_match.tsv", "", organism$V3)
  organism_a<-rbind(organism_a,organism)
}

organism_a$label<- NA
organism_a[is.na(organism_a)]<- "Protein"
colnames(organism_a) <- c("Unknown_gene_loci", "V2", "V3", "name", "label")

for (i in 1:nrow(organism_a)){
  if (str_detect(organism_a[i,4], 'ncRNA')){
    organism_a[i,5] <- "NC - Remove"
  }
  if (str_detect(organism_a[i,4], 'reverse transcriptase')){
    organism_a[i,5] <- "TE - Remove"
  }
  if (str_detect(organism_a[i,4], 'ranspos')){
    organism_a[i,5] <- "TE - Remove"
  }
    if (str_detect(organism_a[i,4], 'etrovirus')){
    organism_a[i,5] <- "TE - Remove"
  }
  if (str_detect(organism_a[i,4], 'ORF')){
    organism_a[i,5] <- "TE - Remove"
  }
  if (str_detect(organism_a[i,4], 'open reading frame')){
    organism_a[i,5] <- "TE - Remove"
  }
  if (str_detect(organism_a[i,4], 'NCBI Reference Sequence')){
    organism_a[i,5] <- "Unknown - Remove"
  }
}


organism_a2 <- organism_a[(organism_a$label != "Protein"),]

proteins_a<-subset(organism_a, organism_a$label == "Protein")

print("Set A have:")
table(proteins_a$label)

write.table(proteins_a, paste0(folder_a, "/functional_characterization/all_hits_set_a.txt"), row.names = F, col.names = F,  sep = "\t")
write.table(organism_a2, paste0(folder_a, "/functional_characterization/unknown_set_a.txt"), row.names = F, col.names = F)

save.image(file=paste0(folder_a,"/R_outputs/Rscript_7.RData"))