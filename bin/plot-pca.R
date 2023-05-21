#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# load packages
library(ape)
library(ggplot2)
library(readr)

# Read command-line arguments
fasta_file <- args[1]
tsv_file <- args[2]

# Read FASTA file
fasta_data <- ape::read.dna(fasta_file, format = "fasta")

# Calculate number of nucleotide differences
dist_matrix <- ape::dist.dna(fasta_data)

# Perform PCA
pca_result <- stats::prcomp(dist_matrix)

# Read TSV table
tsv_data <- readr::read_tsv(tsv_file)

# Merge PCA results with TSV table
merged_data <- merge(pca_result$x, tsv_data, by = 0, all.x = TRUE)

# Plot PCA with colored points based on the table
ggplot(merged_data, aes(PC1, PC2, color = Category)) +
  geom_point() +
  labs(title = "PCA Plot", x = "PC1", y = "PC2")

# Save the plot as a PDF
ggsave("pca_plot.pdf", width = 8.5, height = 10)
