#Rscript8 - Summary

suppressPackageStartupMessages(require('stringr'))
suppressPackageStartupMessages(require('rtracklayer'))

suppressWarnings({
settings <- read.delim("/settings.sh", header = F, stringsAsFactors = F, comment.char = "#")
folder_a <- substring(settings[1,1], 10, )

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

save.image(file=paste0(folder_a, "/R_outputs/Rscript_8.RData"))
})
