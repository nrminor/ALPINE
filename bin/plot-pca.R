#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# load packages
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

# Calculate number of nucleotide differences
dist_matrix <- ape::dist.dna(fasta_data, model = "N")

# Perform PCA
pca_result <- stats::prcomp(dist_matrix)
pca_result <- pca_result$x
pca_result <- cbind(row.names(pca_result), pca_result) ; rownames(pca_result) <- NULL
colnames(pca_result)[1] <- "Accession"

# Read TSV table
tsv_data <- readr::read_tsv(tsv_file,col_names=F)
tsv_data <- tsv_data[tsv_data$X1=="S" | 
                       tsv_data$X1=="H",]
no_seqs <- nrow(tsv_data)
colnames(tsv_data)[2] <- "Cluster"
tsv_data$Cluster <- as.factor(tsv_data$Cluster)
colnames(tsv_data)[9] <- "Accession"

# Merge PCA results with TSV table
merged_data <- merge(pca_result, tsv_data, by = "Accession", all.x = TRUE, sort = TRUE)

# Plot PCA with colored points based on the table
ggplot(merged_data, aes(PC1, PC2, color = Cluster)) +
  geom_point() +
  labs(title = paste("PCA Plot for", no_seqs, "sequences in", yearmonth, sep = " "), 
       x = "PC1", y = "PC2") +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank())

# Save the plot as a PDF
ggsave(paste(yearmonth, "pca_plot.pdf", sep = "_"), width = 8.5, height = 10)

# save the results as a csv
write.csv(merged_data, paste(yearmonth, "pca_results.csv", sep = "_"),
          quote = F, row.names = F)
