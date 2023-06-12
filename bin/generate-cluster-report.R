#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# load necessary packages
library(tidyverse)
library(ape)
library(adephylo)

# Bringing in NCBI metadata
metadata <- read.delim(args[1])

# defining Reference Sequence Accession number
ref_id <- as.character(args[2])

# creating list of year-month combinations to loop through
yearmonths <- str_remove_all(list.files(path = ".", pattern = "*-dist-matrix.csv"), "-dist-matrix.csv")

# make an empty data frame to store high distance FASTA sequences
high_dist_seqs <- data.frame()

# make an empty data frame to store high distance metadata 
high_dist_meta <- matrix(ncol = ncol(metadata), nrow = 0)
colnames(high_dist_meta) <- colnames(metadata)
high_dist_meta <- as.data.frame(high_dist_meta)

# looping through each year-month and adding the highest distance cluster to 
# a final table of hits
first_yearmonth = yearmonths[1]
for (i in yearmonths){
  
  # read in distance matrix
  distmat = read.csv(paste(i, "-dist-matrix.csv", sep = ""))
  
  # define high-distance accession
  max_distance <- 0
  distant_acc <- ""
  for (j in 1:nrow(distmat)){
    
    accession <- distmat[j, "Sequence_Name"]
    distances <- as.numeric(distmat[j,2:ncol(distmat)])
    distances <- distances[distances!=0]
    new_mean <- mean(distances)
    
    # distmat <- distmat[!grepl(ref_id, distmat$Sequence_Name),!grepl(ref_id, colnames(distmat))]
    
    if (grepl(ref_id, accession)){
      next
    } else if (new_mean > max_distance){
      max_distance <- new_mean
      distant_acc <- accession
    } else {
      break
    }
    
  }
  
  # determining which cluster corresponds to this more evolved virus
  cluster_table <- read.delim(paste(i, "-clusters.uc", sep = ""), header = F)
  row_hits <- cluster_table[grepl(distant_acc, cluster_table$V9),]
  if (nrow(row_hits)==0){
    break
  }
  cluster_number <- as.numeric(unique(row_hits$V2))
  
  # creating table for that cluster and adding branch length to it
  high_dist_cluster <- cluster_table[cluster_table$V2==cluster_number,]
  high_dist_cluster$mean_nucleotide_distance <- max_distance
  high_dist_cluster$year_month <- i
  high_dist_cluster <- high_dist_cluster[high_dist_cluster$V1!="C",]
  rownames(high_dist_cluster) <- NULL
  
  # bring in FASTA
  cluster_seqs <- read.delim(paste(i, "-cluster-seqs", cluster_number, sep = ""), header = F)
  
  # add that FASTA to a growing FASTA object of all hits, which will be exported
  # after normalizing later
  high_dist_seqs <- rbind(high_dist_seqs, cluster_seqs) ; remove(cluster_seqs)
  
  # isolate the high distance cluster metadata from NCBI, and save it in its own
  # data frame that will be exported later
  first_accession <- high_dist_cluster$V9[1]
  for (j in high_dist_cluster$V9){
    
    ncbi_row <- metadata[metadata$Accession==j,]
    
    # Add in some useful columns from the vsearch --clust_fast results to the
    # final metadata table
    if (i == first_yearmonth & j == first_accession){
      
      ncbi_row$mean_nucleotide_distance <- high_dist_cluster[high_dist_cluster$V9==j,"mean_nucleotide_distance"]
      ncbi_row$year_month <- high_dist_cluster[high_dist_cluster$V9==j, "year_month"]
      high_dist_meta <- ncbi_row
      
    } else {
      
      ncbi_row$mean_nucleotide_distance <- high_dist_cluster[high_dist_cluster$V9==j, "mean_nucleotide_distance"]
      ncbi_row$year_month <- high_dist_cluster[high_dist_cluster$V9==j, "year_month"]
      high_dist_meta <- rbind(high_dist_meta, ncbi_row)
      
    }
    
  }
  
}

# normalizing year-month by down-prioritizing branches that are not exceptionally
# long compared to all other year-months. I do this here by retaining only sequences that
# are in clusters above the 90% quantile of distances, which pulls out true outliers.
retention_threshold = quantile(high_dist_meta$mean_nucleotide_distance, seq(0, 1, 0.1))[10]
high_dist_meta <- subset(high_dist_meta, mean_nucleotide_distance >= retention_threshold)
high_dist_seqs$keep <- NA
high_dist_seqs$keep <- vapply(high_dist_seqs$V1, function(i){
  
  if (grepl(">", i) && !(str_remove(i, ">") %in% high_dist_meta$Accession)){
    keep_value <- FALSE
  } else if (grepl(">", i) && str_remove(i, ">") %in% high_dist_meta$Accession){
    keep_value <- TRUE
  } else {
    keep_value <- NA
  }
  
  return(keep_value)
  
}, logical(1))
for (i in 1:nrow(high_dist_seqs)){
  
  if (grepl(">", high_dist_seqs$V1[i])){
    next
  } else {
    high_dist_seqs$keep[i] <- high_dist_seqs$keep[i-1]
  }
  
}
high_dist_seqs <- as.data.frame(high_dist_seqs[high_dist_seqs$keep==T,1]) ; colnames(high_dist_seqs) <- NULL

# check here to prevent silent errors
stopifnot(as.numeric(table(grepl(">",high_dist_seqs[,1]))[2]) == nrow(high_dist_meta))

# writing candidate FASTA of high-distance sequences that passed normalization and 
# filtering
write.table(high_dist_seqs, "high_distance_candidates.fasta", row.names = F, col.names = F, quote = F, sep = "\t")

# writing candidate metadata, which combines vsearch --clust_fast results with 
# NCBI metadata and iqTree branch lengths
high_dist_meta <- high_dist_meta[order(high_dist_meta$mean_nucleotide_distance, decreasing = T),]
write.table(high_dist_meta, "high_distance_candidates.tsv", row.names = F, quote = F, sep = "\t")
