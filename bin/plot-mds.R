#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# load packages
library(tidyverse)
library(ape)
library(ggplot2)
library(readr)

# Read command-line arguments
yearmonth <- args[1]
tsv_file <- args[2]
ref_id <- args[3]

# Read FASTA files
fasta_files <- list.files(".", "-cluster-seqs")
for (i in fasta_files){
  
  if (i == fasta_files[1]){
    
    fasta_data <- ape::read.dna(i, format = "fasta")
    
    # Do away with the reference
    fasta_data <- fasta_data[which(!grepl(ref_id, labels(fasta_data))),]
    
  } else {
    
    next_fasta <- ape::read.dna(i, format = "fasta")
    
    # Do away with the reference
    next_fasta <- next_fasta[which(!grepl(ref_id, labels(next_fasta))),]
    
    # merge with out cluster sequences
    fasta_data <- rbind(fasta_data, next_fasta)
    
  }
  
}

if (length(labels(fasta_data)) > 1){
  
  # Calculate number of nucleotide differences
  dist_matrix <- ape::dist.dna(fasta_data, model = "N")
  
  # Perform MDS
  mds_result <- stats::cmdscale(dist_matrix) %>% 
    as.data.frame()
  mds_result <- cbind(row.names(mds_result), mds_result) ; rownames(mds_result) <- NULL
  colnames(mds_result) <- c("Accession", "Dim1", "Dim2")
  
  # Read TSV table
  tsv_data <- readr::read_tsv(tsv_file,col_names=F)
  tsv_data <- tsv_data[tsv_data$X1=="S" | 
                         tsv_data$X1=="H",]
  no_seqs <- nrow(tsv_data)
  colnames(tsv_data)[2] <- "Cluster"
  tsv_data$Cluster <- as.factor(tsv_data$Cluster)
  colnames(tsv_data)[9] <- "Accession"
  
  # Merge MDS results with TSV table
  merged_data <- merge(mds_result, tsv_data, by = "Accession", all.x = TRUE)
  
  # Plot MDS with colored points based on the table
  ggplot(merged_data, aes(x = Dim1, y = Dim2, color = Cluster)) +
    geom_point() +
    labs(title = paste("MDS Plot for", no_seqs, "sequences in", yearmonth, sep = " "), 
         x = "Dimension 1", y = "Dimension 2")
  
  # Save the plot as a PDF
  ggsave(paste(yearmonth, "mds_plot.pdf", sep = "_"), width = 8.5, height = 10)
  
  # Save MDS results as a CSV
  write.csv(merged_data, paste(yearmonth, "mds_results.csv", sep = "_"),
            quote = F, row.names = F)
  
} else {
  
  # Read TSV table
  tsv_data <- readr::read_tsv(tsv_file,col_names=F)
  tsv_data <- tsv_data[tsv_data$X1=="S" | 
                         tsv_data$X1=="H",]
  no_seqs <- nrow(tsv_data)
  colnames(tsv_data)[2] <- "Cluster"
  tsv_data$Cluster <- as.factor(tsv_data$Cluster)
  colnames(tsv_data)[9] <- "Accession"
  
  mds_failed <- data.frame(`MDS Failure Report` = 
                             paste("Accession ", tsv_data$Accession, 
                                   " is alone in its own cluster in ", yearmonth,
                                   " and thus a distance matrix cannot be generated.", sep = ""))
  
  # Save MDS report as a CSV
  write.csv(mds_failed, paste(yearmonth, "mds_failure_report.csv", sep = "_"),
            quote = F, row.names = F)
  
}
