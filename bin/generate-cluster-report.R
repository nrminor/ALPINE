#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# load necessary packages
library(tidyverse)
library(compiler)

# set just-in-time compilation setting so that top-level loops are JIT compiled
compiler::enableJIT(3)

# Bringing in NCBI metadata
metadata <- read.delim(args[1])

# creating list of year-month combinations to loop through
yearmonths <- str_remove_all(list.files(path = ".", pattern = "*-dist-matrix.csv"), "-dist-matrix.csv")

# create function to process metadata
process_meta <- function(yearmonths, metadata){
  
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
    min_distance <- 100.0
    distant_accs <- ""
    least_distant_acc <- ""
    distmat$Sum_weighted_distances <- 0.0
    distmat$Sum_weighted_distances <- vapply(1:nrow(distmat), FUN = 
                                               function(j){
                                                 
                                                 accession <- distmat[j, "Sequence_Name"]
                                                 distances <- as.numeric(distmat[j,2:ncol(distmat)])
                                                 distances <- distances[distances!=0]
                                                 new_sum <- sum(distances)
                                                 
                                                 return(new_sum)
                                                 
                                               }, numeric(1))
    distmat <- distmat[order(distmat$Sequence_Name),] ; rownames(distmat) <- NULL
    
    # collating cluster metadata with distances
    cluster_table <- read.delim(paste(i, "-clusters.uc", sep = ""), header = F)
    cluster_table <- if (nrow(distmat)==1){
      cluster_table <- cluster_table[cluster_table$V9==distmat$Sequence_Name[1],]
    }
    centroids <- cluster_table[cluster_table$V1=="C",]
    centroids <- centroids[order(centroids$V9),] ; rownames(centroids) <- NULL
    stopifnot(nrow(distmat) == nrow(centroids))
    stopifnot(unique(distmat$Sequence_Name == centroids$V9) == TRUE)
    distmat$Cluster_Size <- 0
    distmat$Cluster_Size <- centroids$V3
    distmat <- distmat[order(distmat$Sum_weighted_distances, decreasing = T),] ; rownames(distmat) <- NULL
    distmat$Month <- i 
    distmat$Cluster <- centroids$V2
    
    # new cluster table
    all_seqs <- distmat[, 
                        c("Sequence_Name", "Sum_weighted_distances", 
                          "Cluster_Size", "Month", "Cluster")]
    
    # adding all cluster sequences
    if ("H" %in% cluster_table$V1){
      hits <- cluster_table[cluster_table$V1 == "H",]
      for (j in 1:nrow(hits)){
        
        # record accession and the centroid for that accession
        accession <- hits$V9[j]
        centroid <- hits$V10[j]
        distance <- all_seqs[all_seqs$Sequence_Name==centroid, "Sum_weighted_distances"]
        cluster_size <- all_seqs[all_seqs$Sequence_Name==centroid, "Cluster_Size"]
        month <- i
        cluster <- all_seqs[all_seqs$Sequence_Name==centroid, "Cluster"]
        
        # construct a new row for the all_seqs dataframe
        new_row <- c(accession, distance, cluster_size, month, cluster)
        
        # bind to table
        all_seqs <- rbind(all_seqs, new_row)
        
      }
      
    } else {
      
      # record accession and the centroid for that accession
      accession <- centroids$V9[1]
      centroid <- centroids$V9[1]
      distance <- as.numeric(all_seqs[all_seqs$Sequence_Name==centroid, "Sum_weighted_distances"])
      cluster_size <- as.numeric(all_seqs[all_seqs$Sequence_Name==centroid, "Cluster_Size"])
      month <- i
      cluster <- all_seqs[all_seqs$Sequence_Name==centroid, "Cluster"]
      
      # construct a new row for the all_seqs dataframe
      new_row <- c(accession, distance, cluster_size, month, cluster)
      
      # bind to table
      all_seqs <- rbind(all_seqs, new_row)
      
    }
    
    all_seqs <- all_seqs[order(all_seqs$Sum_weighted_distances),]
    rownames(all_seqs) <- NULL
    
    
    # bring these new data into the NCBI metadata
    for (j in 1:nrow(all_seqs)){
      
      accession <- all_seqs$Sequence_Name[j]
      
      if (!("Sum_weighted_distances" %in% colnames(metadata))){
        
        # create new columns to fill
        metadata$Sum_weighted_distances <- 0.0
        metadata$Cluster_Size <- 0
        metadata$Month_Cluster <- NA
        
      }
      
      # fill in distmat data in the NCBI metadata
      metadata[metadata$Accession==accession, 
               "Sum_weighted_distances"] <- as.numeric(all_seqs$Sum_weighted_distances[j])
      metadata[metadata$Accession==accession, 
               "Cluster_Size"] <- as.numeric(all_seqs$Cluster_Size[j])
      metadata[metadata$Accession==accession, 
               "Month_Cluster"] <- paste(all_seqs$Month[j], all_seqs$Cluster[j], sep = "_")
      
      return(metadata)
      
    }
    
  }
  
}
# compile the metadata function
process_meta_cp <- cmpfun(process_meta)

# run the metadata function in a parallelized fashion
metadata <- process_meta_cp(yearmonths, metadata)

# check to make sure metadata has clustering information in it
stopifnot(nrow(metadata[!is.na(metadata$Month_Cluster),])>0)
metadata <- metadata[metadata$Sum_weighted_distances>0,] ; rownames(metadata) <- NULL

# make an empty data frame to store high distance FASTA sequences
high_dist_seqs <- data.frame()

# bring in FASTAs
fastas <- list.files(path = ".", pattern = "*-cluster-seqs*")
for (fa in fastas){
  
  # bring in FASTA
  cluster_seqs <- read.delim(fa, header = F)
  
  # add that FASTA to a growing FASTA object of all hits, which will be exported
  # after normalizing later
  high_dist_seqs <- rbind(high_dist_seqs, cluster_seqs) ; remove(cluster_seqs)
  
}

# define retention threshold based on the data
retention_threshold = quantile(metadata$Sum_weighted_distances, seq(0, 1, 0.001))[995]

# plot distribution of summed weighted distances
if nrow(metadata > 1){
  pdf("distance_distribution.pdf", width = 7, height = 5.5)
  hist(metadata$Sum_weighted_distances, freq = FALSE, col = "lightblue", ylim = c(0, 0.15), 
       xlab = "Cluster Size Weighted Distances", ylab = "Frequency", 
       main = "Frequency Distribution of Nucleotide Distances")
  lines(density(metadata$Sum_weighted_distances), col = "darkblue", lwd = 2)
  abline(v = retention_threshold, col = "red", lwd = 3)
  text(x = retention_threshold+5, y = (0.15/2), adj = 0,
       labels = paste("Retention Threshold:\n", retention_threshold ))
  dev.off()
}

# normalizing year-month by down-prioritizing branches that are not exceptionally
# long compared to all other year-months. I do this here by retaining only sequences that
# are in clusters above the 90% quantile of distances, which pulls out true outliers.
high_dist_meta <- subset(metadata, Sum_weighted_distances >= retention_threshold)
normalize_meta <- function(high_dist_seqs, high_dist_meta){
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
  
  return(high_dist_seqs)
}

# compile normalize function
normalize_meta_cp <- cmpfun(normalize_meta)

# run the function
high_dist_seqs <- normalize_meta_cp(high_dist_seqs, high_dist_meta)

# check here to prevent silent errors
stopifnot(as.numeric(table(grepl(">",high_dist_seqs[,1]))[2]) == nrow(high_dist_meta))

# writing candidate FASTA of high-distance sequences that passed normalization and 
# filtering
write.table(high_dist_seqs, "high_distance_candidates.fasta", row.names = F, col.names = F, quote = F, sep = "\t")

# writing candidate metadata, which combines vsearch --clust_fast results with 
# NCBI metadata and iqTree branch lengths
high_dist_meta <- high_dist_meta[order(high_dist_meta$Sum_weighted_distances, decreasing = T),]
write.table(high_dist_meta, "high_distance_candidates.tsv", row.names = F, quote = F, sep = "\t")
