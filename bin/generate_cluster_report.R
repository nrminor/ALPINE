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
yearmonths <- str_remove_all(list.files(path = ".", pattern = "*.treefile"), ".treefile")

# make an empty data frame to store high distance FASTA sequences
high_dist_seqs <- data.frame("V1" = as.character(NA))
colnames(high_dist_seqs) <- NULL
high_dist_seqs <- high_dist_seqs[NULL,]

# make an empty data frame to store high distance metadata 
high_dist_meta <- matrix(ncol = ncol(metadata), nrow = 0)
colnames(high_dist_meta) <- colnames(metadata)
high_dist_meta <- as.data.frame(high_dist_meta)

# looping through each year-month and adding the highest distance cluster to 
# a final table of hits
first_yearmonth = yearmonths[1]
for (i in yearmonths){
  
  # read in the tree file
  tree <- read.tree(paste(i, ".treefile", sep = ""))
  
  # root the tree, if it isn't already
  if (is.rooted(tree)==F){
    tree <- root(tree, ref_id, resolve.root = T)
  }
  
  # fixing non-unique node labels
  node_support <- tree$node.label
  tree$node.label <- 1:length(node_support)
  
  # compute total branch lengths from root to each tip
  roottotip <- distRoot(tree, tips = tree$tip.label, method = "patristic")
  roottotip <- roottotip[names(roottotip)!=ref_id]
  
  # finding max branch length sample
  max_distance <- roottotip[roottotip==max(roottotip)]
  
  # restoring original isolate name formatting
  if (grepl("hCoV-19", names(max_distance))){
    split_name <- unlist(strsplit(names(max_distance), "_EPI"))[1]
    name <- str_replace_all(split_name, "_", "/")
  } else {
    name <- names(max_distance)
  }
  
  # determining which cluster corresponds to this more evolved virus
  cluster_table <- read.delim(paste(i, "-clusters.uc", sep = ""), header = F)
  row_hits <- cluster_table[grepl(name, cluster_table$V9),]
  cluster_number <- as.numeric(unique(row_hits$V2))
  
  # creating table for that cluster and adding branch length to it
  high_dist_cluster <- cluster_table[cluster_table$V2==cluster_number,]
  high_dist_cluster$root_to_tip_distance <- max(roottotip)
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
      
      high_dist_meta <- rbind(high_dist_meta, ncbi_row)
      high_dist_meta <- cbind(high_dist_meta, high_dist_cluster[high_dist_cluster$V9==j, "root_to_tip_distance"])
      colnames(high_dist_meta)[ncol(high_dist_meta)] <- "root_to_tip_distance"
      high_dist_meta <- cbind(high_dist_meta, high_dist_cluster[high_dist_cluster$V9==j, "year_month"])
      colnames(high_dist_meta)[ncol(high_dist_meta)] <- "year_month"
      
    } else {
      
      ncbi_row$root_to_tip_distance <- high_dist_cluster[high_dist_cluster$V9==j, "root_to_tip_distance"]
      ncbi_row$year_month <- high_dist_cluster[high_dist_cluster$V9==j, "year_month"]
      high_dist_meta <- rbind(high_dist_meta, ncbi_row)
      
    }
    
  }
  
}

# normalizing year-month by down-prioritizing branches that are not exceptionally
# long compared to all other year-months

# NOTE TO SELF: read in distance matrices here

# writing candidate FASTA of high-distance sequences that passed normalization and 
# filtering
write.table(high_dist_seqs, "high_distance_candidates.fasta", row.names = F, col.names = F, quote = F, sep = "\t")

# writing candidate metadata, which combines vsearch --clust_fast results with 
# NCBI metadata and iqTree branch lengths
high_dist_meta <- high_dist_meta[order(high_dist_meta$root_to_tip_distance, decreasing = T),]
write.table(high_dist_meta, "high_distance_candidates.tsv", row.names = F, quote = F, sep = "\t")
