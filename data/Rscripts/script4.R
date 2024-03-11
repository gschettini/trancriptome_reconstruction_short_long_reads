#Rscript4 - Filter

suppressPackageStartupMessages(require('stringr'))

settings <- read.delim("/settings.sh", header = F, stringsAsFactors = F, comment.char = "#")
folder_a <- substring(settings[1,1], 10, )
folder_b <- substring(settings[2,1], 10, )


#Set A
lncRNA_a<-read.table(paste0(folder_a,"/comparisons/lncRNA/set_a_lncRNA.annotated.gtf"), sep ="\t")

lncRNA_flag_u_a<-read.table(paste0(folder_a,"/comparisons/lncRNA/unknown_set_a_lncRNA.annotated.gtf"), sep ="\t")

lncRNA_flag_u_a$gene_name <- word(lncRNA_flag_u_a$V9, start = 3, end = 4, sep = " ")
lncRNA_flag_u_a$gene_name <- paste(lncRNA_flag_u_a$V1, lncRNA_flag_u_a$gene_name)

lncRNA_transcript_a<-subset(lncRNA_a, lncRNA_a$V3 == "transcript")

lncRNA_transcript_a$gene_name <- word(lncRNA_transcript_a$V9, start = 3, end = 4, sep = " ")
lncRNA_transcript_a$gene_name <- paste(lncRNA_transcript_a$V1, lncRNA_transcript_a$gene_name)

lncRNA_transcript_a$keep<-lncRNA_transcript_a$gene_name %in% lncRNA_flag_u_a$gene_name

lncRNA_transcript_a2<-subset(lncRNA_transcript_a, lncRNA_transcript_a$keep == TRUE)
lncRNA_exons_a<-subset(lncRNA_a, lncRNA_a$V3 == "exon")

lncRNA_exons_a$gene_name <- word(lncRNA_exons_a$V9, start = 3, end = 4, sep = " ")
lncRNA_exons_a$gene_name <- paste(lncRNA_exons_a$V1, lncRNA_exons_a$gene_name)

lncRNA_exons_a$keep<-lncRNA_exons_a$gene_name %in% lncRNA_transcript_a2$gene_name
lncRNA_exons_a2<-subset(lncRNA_exons_a, lncRNA_exons_a$keep == TRUE)

lncRNA_a$keep<-lncRNA_a$V9 %in% lncRNA_transcript_a2$V9
lncRNA_a$keep2<-lncRNA_a$V9 %in% lncRNA_exons_a2$V9

lncRNA_a2<-subset(lncRNA_a, lncRNA_a$keep == TRUE | lncRNA_a$keep2 == TRUE)
lncRNA_a3<-lncRNA_a2[,-c(10:11)]


#Set B
lncRNA_b<-read.table(paste0(folder_b,"/comparisons/lncRNA/set_b_lncRNA.annotated.gtf"), sep ="\t")

lncRNA_flag_u_b<-read.table(paste0(folder_b,"/comparisons/lncRNA/unknown_set_b_lncRNA.annotated.gtf"), sep ="\t")

lncRNA_flag_u_b$gene_name <- word(lncRNA_flag_u_b$V9, start = 3, end = 4, sep = " ")
lncRNA_flag_u_b$gene_name <- paste(lncRNA_flag_u_b$V1, lncRNA_flag_u_b$gene_name)

lncRNA_transcript_b<-subset(lncRNA_b, lncRNA_b$V3 == "transcript")

lncRNA_transcript_b$gene_name <- word(lncRNA_transcript_b$V9, start = 3, end = 4, sep = " ")
lncRNA_transcript_b$gene_name <- paste(lncRNA_transcript_b$V1, lncRNA_transcript_b$gene_name)

lncRNA_transcript_b$keep<-lncRNA_transcript_b$gene_name %in% lncRNA_flag_u_b$gene_name

lncRNA_transcript_b2<-subset(lncRNA_transcript_b, lncRNA_transcript_b$keep == TRUE)
lncRNA_exons_b<-subset(lncRNA_b, lncRNA_b$V3 == "exon")

lncRNA_exons_b$gene_name <- word(lncRNA_exons_b$V9, start = 3, end = 4, sep = " ")
lncRNA_exons_b$gene_name <- paste(lncRNA_exons_b$V1, lncRNA_exons_b$gene_name)

lncRNA_exons_b$keep<-lncRNA_exons_b$gene_name %in% lncRNA_transcript_b2$gene_name
lncRNA_exons_b2<-subset(lncRNA_exons_b, lncRNA_exons_b$keep == TRUE)

lncRNA_b$keep<-lncRNA_b$V9 %in% lncRNA_transcript_b2$V9
lncRNA_b$keep2<-lncRNA_b$V9 %in% lncRNA_exons_b2$V9

lncRNA_b2<-subset(lncRNA_b, lncRNA_b$keep == TRUE | lncRNA_b$keep2 == TRUE)
lncRNA_b3<-lncRNA_b2[,-c(10:11)]


#summary clean GFF file lncRNA

print("Set A summary")
nrow(subset(lncRNA_a3, lncRNA_a3$V3 == "transcript"))
print("Set B summary")
nrow(subset(lncRNA_b3, lncRNA_b3$V3 == "transcript"))

write.table(lncRNA_a3, paste0(folder_a, "/comparisons/lncRNA/set_a_lncRNA_filtered.gff3"), quote = FALSE, col.names = FALSE, sep = "\t", row.names = FALSE)
write.table(lncRNA_b3, paste0(folder_b, "/comparisons/lncRNA/set_b_lncRNA_filtered.gff3"), quote = FALSE, col.names = FALSE, sep = "\t", row.names = FALSE)

save.image(file=paste0(folder_a,"/R_outputs/Rscript_4.RData"))
save.image(file=paste0(folder_b,"/R_outputs/Rscript_4.RData"))

print("4th filtering finished")
