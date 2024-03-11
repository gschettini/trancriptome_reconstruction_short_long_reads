#Rscript2 - Filter

suppressPackageStartupMessages(require('stringr'))

settings <- read.delim("/settings.sh", header = F, stringsAsFactors = F, comment.char = "#")
folder_a <- substring(settings[1,1], 10, )

#set A
GFF_GMAP_a<-read.table(paste0(folder_a,"/annotation_GMAP/ncbi/set_a_ncbi.gff3"), sep = "\t")

GFF_GMAP_a$gene<-word(GFF_GMAP_a$V9, start = 1, end = 1, sep = ";")
GFF_GMAP_a$gene<-gsub("ID=", "", GFF_GMAP_a$gene)
GFF_GMAP_a$gene<-gsub(".path.*", "", GFF_GMAP_a$gene)
GFF_GMAP_a$gene<-gsub(".mrna.*", "", GFF_GMAP_a$gene)

GFF_filtered_a<-read.table(paste0(folder_a,"/annotation_GMAP/ncbi/set_a_ncbi_filtered.gff3"), sep = "\t")
GFF_filtered_a$gene<-word(GFF_filtered_a$V9, start = 1, end = 1, sep = ";")
GFF_filtered_a$gene<-gsub("ID=", "", GFF_filtered_a$gene)
GFF_filtered_a$gene<-gsub(".path.*", "", GFF_filtered_a$gene)
GFF_filtered_a$gene<-gsub(".mrna.*", "", GFF_filtered_a$gene)

GFF_GMAP2_a<-subset(GFF_GMAP_a, GFF_GMAP_a$gene %in% GFF_filtered_a$gene)
GFF_GMAP2_a<- GFF_GMAP2_a[,c(1:9)]

#   summary clean GFF file NCBI

print("Set A summary")
nrow(subset(GFF_GMAP2_a, GFF_GMAP2_a$V3 == "gene"))

write.table(GFF_GMAP2_a, paste0(folder_a, "/annotation_GMAP/ncbi/set_a_ncbi_filtered_2.gff3"), quote = FALSE, col.names = FALSE, sep = "\t", row.names = FALSE)

save.image(file=paste0(folder_a,"/R_outputs/Rscript_2.RData"))

print("2nd filtering finished")
