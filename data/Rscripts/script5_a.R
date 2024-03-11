#Rscript5 - Filter

suppressPackageStartupMessages(require('stringr'))
suppressPackageStartupMessages(require('rtracklayer'))
suppressPackageStartupMessages(require('GenomicRanges'))

settings <- read.delim("/settings.sh", header = F, stringsAsFactors = F, comment.char = "#")
folder_a <- substring(settings[1,1], 10, )
gtf_ensembl <- substring(settings[9,1], 13, )

#Set A
lncRNA_4_a <- read.table(paste0(folder_a,"/comparisons/lncRNA/set_a_all_flags_u_3.gtf"), sep ="\t")
lncRNA_4_a_transcripts <- subset(lncRNA_4_a, lncRNA_4_a$V3 == "transcript")
lncRNA_4_a_transcripts$V1<-as.factor(lncRNA_4_a_transcripts$V1)
lncRNA_4_a_transcripts$V7<-as.factor(lncRNA_4_a_transcripts$V7)


lncRNA_4_a_transcripts_granges <- GRanges(
  seqnames = lncRNA_4_a_transcripts$V1,
  ranges = IRanges(
    start = lncRNA_4_a_transcripts$V4,
    end = lncRNA_4_a_transcripts$V5,
    names = lncRNA_4_a_transcripts$V9
  ),
  strand = lncRNA_4_a_transcripts$V7
)

dist_filter<-5000

official_annotation_ensembl<-as.data.frame(rtracklayer::import(gtf_ensembl))

official_annotation_ensembl_granges <- GRanges(
  seqnames = official_annotation_ensembl$seqnames,
  ranges = IRanges(
    start = official_annotation_ensembl$start,
    end = official_annotation_ensembl$end
  ),
  strand = official_annotation_ensembl$strand
)

no_overlaps_a<- !overlapsAny(lncRNA_4_a_transcripts_granges, official_annotation_ensembl_granges,
    maxgap=dist_filter,
    ignore.strand=TRUE)

no_overlapping_transcripts_a <- lncRNA_4_a_transcripts_granges[(no_overlaps_a)]

no_overlapping_transcripts_a_df <- data.frame(
  seqnames = seqnames(no_overlapping_transcripts_a),
  start = start(no_overlapping_transcripts_a),
  end = end(no_overlapping_transcripts_a),
  strand = strand(no_overlapping_transcripts_a),
  names = names(no_overlapping_transcripts_a)
)

lncRNA_4_a_transcripts$keep<- lncRNA_4_a_transcripts$V9 %in% no_overlapping_transcripts_a_df$names
lncRNA_4_a_transcripts<-subset(lncRNA_4_a_transcripts, lncRNA_4_a_transcripts$keep == TRUE)
lncRNA_4_a_exons <- subset(lncRNA_4_a, lncRNA_4_a$V3 == "exon")

lncRNA_4_a_transcripts$gene_name <- word(lncRNA_4_a_transcripts$V9, start = 1, end = 1, sep = ";")
lncRNA_4_a_transcripts$gene_name <- gsub('ID=','',lncRNA_4_a_transcripts$gene_name)
lncRNA_4_a_transcripts$gene_name <- gsub(';','',lncRNA_4_a_transcripts$gene_name)
lncRNA_4_a_exons$gene_name <- lncRNA_4_a_exons$V9
lncRNA_4_a_exons$gene_name <- word(lncRNA_4_a_exons$gene_name, start = 1, end = 1, sep = ";")
lncRNA_4_a_exons$gene_name <- gsub('Parent=','',lncRNA_4_a_exons$gene_name)
lncRNA_4_a_exons$V1<-as.factor(lncRNA_4_a_exons$V1)
lncRNA_4_a_transcripts_exons<-data.frame()

lncRNA_4_a_transcripts_exons <- merge(lncRNA_4_a_transcripts[,c(1, 3:5, 11)], lncRNA_4_a_exons[,c(1, 3:5, 10)], by= "gene_name",  all = TRUE)

lncRNA_4_a_transcripts_exons <- lncRNA_4_a_transcripts_exons[!is.na(lncRNA_4_a_transcripts_exons$V1.x), ]

lncRNA_4_a_transcripts_exons$remove = NA


for (i in 1:nrow(lncRNA_4_a_transcripts_exons)) {
  if (lncRNA_4_a_transcripts_exons[i,2] == lncRNA_4_a_transcripts_exons[i,6] &
      lncRNA_4_a_transcripts_exons[i,4] == lncRNA_4_a_transcripts_exons[i,8] &
      lncRNA_4_a_transcripts_exons[i,5] == lncRNA_4_a_transcripts_exons[i,9]) {
    lncRNA_4_a_transcripts_exons[i,"remove"]<- TRUE
  }
  else {
    lncRNA_4_a_transcripts_exons[i,"remove"]<- FALSE
  }
}

lncRNA_4_a_transcripts_exons <- subset(lncRNA_4_a_transcripts_exons, lncRNA_4_a_transcripts_exons$remove == FALSE)

lncRNA_4_a$gene_name <- word(lncRNA_4_a$V9, start = 1, end = 1, sep = ";")
lncRNA_4_a$gene_name <- gsub('ID=','',lncRNA_4_a$gene_name)
lncRNA_4_a$gene_name <- gsub('Parent=','',lncRNA_4_a$gene_name)

lncRNA_4_a$remove <- NA

for (i in 1:nrow(lncRNA_4_a)){
  if (lncRNA_4_a[i, 10] %in% lncRNA_4_a_transcripts_exons$gene_name){
    lncRNA_4_a[i, "remove"] <- FALSE
  }
  else {
    lncRNA_4_a[i, "remove"] <- TRUE
  }
}

lncRNA_4_a <- subset(lncRNA_4_a, lncRNA_4_a$remove == FALSE)

lncRNA_4_a <- lncRNA_4_a[,1:9]


##Summary
only_transcripts_a <- subset(lncRNA_4_a, lncRNA_4_a$V3 == "transcript")
print("Unknwon transcripts - Set A")
nrow(only_transcripts_a)

write.table(lncRNA_4_a, paste0(folder_a, "/comparisons/lncRNA/set_a_all_flags_u_clean.gtf"), col.names = F, row.names = F, quote = F, sep = "\t")

save.image(file=paste0(folder_a,"/R_outputs/Rscript_5.RData"))

print("5th filtering finished")
