#Rscript8 - Summary

suppressPackageStartupMessages(require('stringr'))
suppressPackageStartupMessages(require('rtracklayer'))
suppressPackageStartupMessages(require('GenomicRanges'))

settings <- read.delim("/settings.sh", header = F, stringsAsFactors = F, comment.char = "#")
folder_a <- substring(settings[1,1], 10, )
folder_b <- substring(settings[2,1], 10, )


#Set A

rnasamba_a <- read.table(paste0(folder_a, "/functional_characterization/output.tsv"), header = T)
rnasamba_a$sequence_name <- rownames(rnasamba_a)
rnasamba_a <- rnasamba_a[,c(1,3)]
proteins_a <- read.table(paste0(folder_a, "/functional_characterization/all_hits_set_a.txt"))
proteins_a <- proteins_a[,c(1,2,4,5)]
colnames(proteins_a) <- c("sequence_name", "NCBI_ID", "name", "classification")
summary_a <- merge(proteins_a, rnasamba_a, by="sequence_name", all = T)
summary_a$final_classification<-NA
for (i in 1:nrow(summary_a)){
  if (is.na(summary_a[i,"classification.x"]) &&
  summary_a[i,"classification.y"] == "coding"){
    summary_a[i,"final_classification"] <- "C_RS"
  }
  if (is.na(summary_a[i,"classification.x"]) &&
  summary_a[i,"classification.y"] == "noncoding"){
    summary_a[i,"final_classification"] <- "NC_RS"
  }
  if (summary_a[i,"classification.x"] == "Protein" &&
  is.na(summary_a[i,"classification.y"])) {
    summary_a[i,"final_classification"] <- "C_Or"
  }
}

summary_a <- summary_a[,c(1:3,6)]
write.csv(summary_a, paste0(folder_a, "/functional_characterization/summary_transcript_loci_set_a"), row.names = FALSE)

#Set B

rnasamba_b <- read.table(paste0(folder_b, "/functional_characterization/output.tsv"), header = T)
rnasamba_b$sequence_name <- rownames(rnasamba_b)
rnasamba_b <- rnasamba_b[,c(1,3)]
proteins_b <- read.table(paste0(folder_b, "/functional_characterization/all_hits_set_b.txt"))
proteins_b <- proteins_b[,c(1,2,4,5)]
colnames(proteins_b) <- c("sequence_name", "NCBI_ID", "name", "classification")
summary_b <- merge(proteins_b, rnasamba_b, by="sequence_name", all = T)
summary_b$final_classification<-NA
for (i in 1:nrow(summary_b)){
  if (is.na(summary_b[i,"classification.x"]) &&
  summary_b[i,"classification.y"] == "coding"){
    summary_b[i,"final_classification"] <- "C_RS"
  }
  if (is.na(summary_b[i,"classification.x"]) &&
  summary_b[i,"classification.y"] == "noncoding"){
    summary_b[i,"final_classification"] <- "NC_RS"
  }
  if (summary_b[i,"classification.x"] == "Protein" &&
  is.na(summary_b[i,"classification.y"])) {
    summary_b[i,"final_classification"] <- "C_Or"
  }
}

summary_b <- summary_b[,c(1:3,6)]
write.csv(summary_b, paste0(folder_b, "/functional_characterization/summary_transcript_loci_set_b"), row.names = FALSE)


#Set A

gtf_file_a <- read.table(paste0(folder_a, "/comparisons/lncRNA/set_a_all_flags_u_final3.gtf"), sep = "\t")
gtf_file_a <- subset(gtf_file_a, gtf_file_a$V3 == "transcript")

gtf_file_a$name <- word(sep = ";", start = 1 , end = 1, gtf_file_a$V9)
gtf_file_a$name <- gsub("transcript_id ", "", gtf_file_a$name)

gtf_file_a_granges <- GRanges(
  seqnames = gtf_file_a$V1,
  ranges = IRanges(
    start = gtf_file_a$V4,
    end = gtf_file_a$V5,
    names = gtf_file_a$V9
    ),
  strand = gtf_file_a$V7,
  names = gtf_file_a$V9
)


#Set B

gtf_file_b <- read.table(paste0(folder_b, "/comparisons/lncRNA/set_b_all_flags_u_final3.gtf"), sep = "\t")
gtf_file_b <- subset(gtf_file_b, gtf_file_b$V3 == "transcript")

gtf_file_b$name <- word(sep = ";", start = 1 , end = 1, gtf_file_b$V9)
gtf_file_b$name <- gsub("transcript_id ", "", gtf_file_b$name)

gtf_file_b_granges <- GRanges(
  seqnames = gtf_file_b$V1,
  ranges = IRanges(
    start = gtf_file_b$V4,
    end = gtf_file_b$V5,
    names = gtf_file_b$V9
    ),
  strand = gtf_file_b$V7,
  names = gtf_file_b$V9
)

#Overlapping Set A vs Set B

overlaps_a<- suppressWarnings(overlapsAny(gtf_file_a_granges, gtf_file_b_granges,
    maxgap=0,
    ignore.strand=FALSE))

overlapping_a_gtf <- gtf_file_a_granges[(overlaps_a)]

overlapping_a_gtf_df <- data.frame(
  seqnames = seqnames(overlapping_a_gtf),
  start = start(overlapping_a_gtf),
  end = end(overlapping_a_gtf),
  strand = strand(overlapping_a_gtf),
  names = names(overlapping_a_gtf)
)

b_overlapping_a_gtf_df <- overlapping_a_gtf_df
b_overlapping_a_gtf_df$names <- word(sep = ";", start = 1, end = 1, b_overlapping_a_gtf_df$names)
b_overlapping_a_gtf_df$names <- gsub("transcript_id ", "", b_overlapping_a_gtf_df$names)
b_overlapping_a_gtf_df$seqnames <- as.character(b_overlapping_a_gtf_df$seqnames)
b_overlapping_a_gtf_df <- as.data.frame(b_overlapping_a_gtf_df[order(b_overlapping_a_gtf_df$seqnames, b_overlapping_a_gtf_df$start),])

write.table(b_overlapping_a_gtf_df, paste0(folder_a, "/R_outputs/set_a_transcripts_in_common_with_set_b.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)

no_overlaps_a<- suppressWarnings(!overlapsAny(gtf_file_a_granges, gtf_file_b_granges,
    maxgap=0,
    ignore.strand=FALSE))

no_overlapping_a_gtf <- gtf_file_a_granges[(no_overlaps_a)]

no_overlapping_a_gtf_df <- data.frame(
  seqnames = seqnames(no_overlapping_a_gtf),
  start = start(no_overlapping_a_gtf),
  end = end(no_overlapping_a_gtf),
  strand = strand(no_overlapping_a_gtf),
  names = names(no_overlapping_a_gtf)
)

b_no_overlapping_a_gtf_df <- no_overlapping_a_gtf_df
b_no_overlapping_a_gtf_df$names <- word(sep = ";", start = 1, end = 1, b_no_overlapping_a_gtf_df$names)
b_no_overlapping_a_gtf_df$names <- gsub("transcript_id ", "", b_no_overlapping_a_gtf_df$names)
b_no_overlapping_a_gtf_df$seqnames <- as.character(b_no_overlapping_a_gtf_df$seqnames)
b_no_overlapping_a_gtf_df <- b_no_overlapping_a_gtf_df[order(b_no_overlapping_a_gtf_df$seqnames, b_no_overlapping_a_gtf_df$start),]


write.table(b_no_overlapping_a_gtf_df, paste0(folder_a, "/R_outputs/set_a_transcripts_exclusive.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)

#Overlapping Set B vs Set A


overlaps_b<- suppressWarnings(overlapsAny(gtf_file_b_granges, gtf_file_a_granges,
    maxgap=0,
    ignore.strand=FALSE))

overlapping_b_gtf <- gtf_file_b_granges[(overlaps_b)]

overlapping_b_gtf_df <- data.frame(
  seqnames = seqnames(overlapping_b_gtf),
  start = start(overlapping_b_gtf),
  end = end(overlapping_b_gtf),
  strand = strand(overlapping_b_gtf),
  names = names(overlapping_b_gtf)
)

b_overlapping_b_gtf_df <- overlapping_b_gtf_df
b_overlapping_b_gtf_df$names <- word(sep = ";", start = 1, end = 1, b_overlapping_b_gtf_df$names)
b_overlapping_b_gtf_df$names <- gsub("transcript_id ", "", b_overlapping_b_gtf_df$names)
b_overlapping_b_gtf_df$seqnames <- as.character(b_overlapping_b_gtf_df$seqnames)
b_overlapping_b_gtf_df <- b_overlapping_b_gtf_df[order(b_overlapping_b_gtf_df$seqnames, b_overlapping_b_gtf_df$start),]


write.table(b_overlapping_b_gtf_df, paste0(folder_b, "/R_outputs/set_b_transcripts_in_common_with_set_a.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)

no_overlaps_b<- suppressWarnings(!overlapsAny(gtf_file_b_granges, gtf_file_a_granges,
    maxgap=0,
    ignore.strand=FALSE))

no_overlapping_b_gtf <- gtf_file_b_granges[(no_overlaps_b)]

no_overlapping_b_gtf_df <- data.frame(
  seqnames = seqnames(no_overlapping_b_gtf),
  start = start(no_overlapping_b_gtf),
  end = end(no_overlapping_b_gtf),
  strand = strand(no_overlapping_b_gtf),
  names = names(no_overlapping_b_gtf)
)

b_no_overlapping_b_gtf_df <- no_overlapping_b_gtf_df
b_no_overlapping_b_gtf_df$names <- word(sep = ";", start = 1, end = 1, b_no_overlapping_b_gtf_df$names)
b_no_overlapping_b_gtf_df$names <- gsub("transcript_id ", "", b_no_overlapping_b_gtf_df$names)
b_no_overlapping_b_gtf_df$seqnames <- as.character(b_no_overlapping_b_gtf_df$seqnames)
b_no_overlapping_b_gtf_df <- b_no_overlapping_b_gtf_df[order(b_no_overlapping_b_gtf_df$seqnames, b_no_overlapping_b_gtf_df$start),]


write.table(b_no_overlapping_b_gtf_df, paste0(folder_b, "/R_outputs/set_b_transcripts_exclusive.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)


save.image(file=paste0(folder_a, "/R_outputs/Rscript_8.RData"))
save.image(file=paste0(folder_b, "/R_outputs/Rscript_8.RData"))

print("Total number of common trancript loci:")
nrow(b_overlapping_a_gtf_df)
print("Specific from set A:")
nrow(b_no_overlapping_a_gtf_df)
print("Specific from set B:")
nrow(b_no_overlapping_b_gtf_df)
