#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# load necessary libraries
library(tidyverse)
library(Biostrings)
library(compiler)

# record input file path for cluster table
cluster_table_path <- args[1]
metadata_path <- args[2]

# make sure vsearch cluster FASTA files are present in the current working directory
stopifnot(TRUE %in% grepl("meta-cluster-seqs", list.files()))

# define main function
main <- function(cluster_table_path, metadata_path){
  
  # read the cluster table into a dataframe and filter to repeats only
  cluster_df <- read_tsv(cluster_table_path, 
                         show_col_types = FALSE, trim_ws = TRUE, col_names = FALSE) %>%
    filter(X1!="C") %>%
    group_by(X2) %>%
    filter(n() > 1) %>%
    ungroup() %>%
    arrange(X2)
  
  stopifnot(length(unique(cluster_df$X9))==nrow(cluster_df))
  stopifnot(nrow(cluster_df)>2)
  
  # compile FASTA sequences for the above clusters
  seq_tally <- vapply(1:length(unique(cluster_df$X2)), FUN = function(i){
    
    cluster <- unique(cluster_df$X2)[i]
    fasta <- readDNAStringSet(paste("meta-cluster-seqs", cluster, sep = ""))
    writeXStringSet(fasta, 
                    paste("repeat-lineage-", i, ".fasta", sep = ""))
    
    return(length(labels(fasta)))
    
  }, integer(1))
  seq_tally <- sum(seq_tally)
  stopifnot(seq_tally==length(cluster_df$X2))
  
  # create new table to store accessions and repeat cluster number
  repeat_clusts <- tibble("Accession" = cluster_df$X9,
                          "Repeat Cluster Number" = match(cluster_df$X2, unique(cluster_df$X2)))
  
  # collate metadata about these new clusters
  metadata <- read_tsv(metadata_path, show_col_types = FALSE)
  merged_metadata <- merge(metadata, repeat_clusts, by = "Accession", all.x = TRUE) %>%
    filter(!is.na(`Repeat Cluster Number`)) %>%
    arrange(`Repeat Cluster Number`) %>%
    select(!c("Sum_weighted_distances", "Cluster_Size"))
  stopifnot(nrow(merged_metadata)==nrow(cluster_df))
  
  # write repeat lineage metadata
  write_tsv(merged_metadata, "repeat-lineage-metadata.tsv", na = "")
  
}
main_cp <- cmpfun(main)

# run the main function
main_cp(cluster_table_path, metadata_path)
