#Rscript3 - Filter

suppressPackageStartupMessages(require('stringr'))

settings <- read.delim("/settings.sh", header = F, stringsAsFactors = F, comment.char = "#")
folder_a <- substring(settings[1,1], 10, )

#Set A
ncbi_a<-read.table(paste0(folder_a,"/comparisons/ncbi/set_a_ncbi.annotated.gtf"), sep ="\t")

ncbi_flag_u_a<-read.table(paste0(folder_a,"/comparisons/ncbi/unknown_set_a_ncbi.annotated.gtf"), sep ="\t")

un_chr_a<-ncbi_flag_u_a[startsWith(ncbi_flag_u_a$V1, "Leftover"),]
un_chr_a<-rbind(un_chr_a,ncbi_flag_u_a[startsWith(ncbi_flag_u_a$V1, "MT"),] )

ncbi_flag_u_a$remove<-ncbi_flag_u_a$V1 %in% un_chr_a$V1

ncbi_flag_u_a<-subset(ncbi_flag_u_a, ncbi_flag_u_a$remove == FALSE)

ncbi_flag_u_a$gene_name <- word(ncbi_flag_u_a$V9, start = 3, end = 4, sep = " ")
ncbi_flag_u_a$gene_name <- paste(ncbi_flag_u_a$V1, ncbi_flag_u_a$gene_name)

ncbi_exons_a<-subset(ncbi_a, ncbi_a$V3 == "exon")

ncbi_exons_a$gene_name <- word(ncbi_exons_a$V9, start = 3, end = 4, sep = " ")
ncbi_exons_a$gene_name <- paste(ncbi_exons_a$V1, ncbi_exons_a$gene_name)

ncbi_exons_a$keep<-ncbi_exons_a$gene_name %in% ncbi_flag_u_a$gene_name

ncbi_exons_a2<-subset(ncbi_exons_a, ncbi_exons_a$keep == TRUE)

ncbi_a$keep<-ncbi_a$V9 %in% ncbi_flag_u_a$V9

ncbi_a$keep2<-ncbi_a$V9 %in% ncbi_exons_a2$V9

ncbi_a2<-subset(ncbi_a, ncbi_a$keep == TRUE | ncbi_a$keep2 == TRUE)

ncbi_a3<-ncbi_a2[,-c(10:11)]

write.table(ncbi_a3, paste0(folder_a, "/comparisons/ncbi/set_a_ncbi_filtered_3.gff3"), sep = "\t" , col.names = FALSE, row.names =  FALSE, quote = FALSE)

#summary clean GFF file NCBI

print("Set A summary")
nrow(subset(ncbi_a3, ncbi_a3$V3 == "transcript"))

save.image(file=paste0(folder_a,"/R_outputs/Rscript_3.RData"))

print("3rd filtering finished")