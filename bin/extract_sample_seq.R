#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

fasta <- read.delim(as.character(args[1]), header = F)
seq_limits <- data.frame(deflines = fasta[which(grepl(">", fasta$V1, fixed = T)),],
						 deflines_indices = which(grepl(">", fasta$V1, fixed = T)),
						 seq_ends = c(which(grepl(">", fasta$V1, fixed = T))[2]-1,
									  nrow(fasta)))
to_keep <- seq_limits[seq_limits$deflines!=">BA.2",]

fasta <- fasta[to_keep$deflines_indices[1]:to_keep$seq_ends[1],]

write.table(fasta, paste(args[2], ".fasta", sep = ""), quote = F, row.names = F, col.names = F)
