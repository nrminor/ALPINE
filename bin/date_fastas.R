#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(tidyverse)

fasta_path = args[1]

fasta_defline <- read.delim(fasta_path, header = F, sep = "|", nrows = 1)
fasta_defline <- fasta_defline[,-which(grepl("EPI",t(fasta_defline)))]
fasta_defline$V1 <- str_remove(fasta_defline$V1, ">")
fasta_defline$V1 <- str_replace_all(fasta_defline$V1, "/", "_")
fasta_defline <- fasta_defline[,1:2]
fasta_defline <- cbind(fasta_path, fasta_defline)
colnames(fasta_defline) <- NULL

write.table(fasta_defline, "fasta_date_tuple.csv", 
          sep = ",", quote = F, row.names = F, col.names = F)
