#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# load necessary packages
require(tidyverse)
require(arrow)
require(Biostrings)

# set just-in-time compilation setting so that top-level loops are JIT compiled
compiler::enableJIT(3)

# Bringing in NCBI metadata
metadata <- read_ipc_file(args[1])

# Defining strictness level
if (args[2] == "strict") {
  stringency <- 995
} else if (args[2] == "intermediate") {
  stringency <- 990
} else if (args[2] == "lenient") {
  stringency <- 980
} else {
  stringency <- 995
}

# create function to process metadata
process_meta <- compiler::cmpfun(function(yearmonths, metadata){
  
  # make an empty data frame to store high distance metadata 
  high_dist_meta <- matrix(ncol = ncol(metadata), nrow = 0)
  colnames(high_dist_meta) <- colnames(metadata)
  high_dist_meta <- as.data.frame(high_dist_meta)
  
  # looping through each year-month and adding the highest distance cluster to 
  # a final table of hits
  first_yearmonth = yearmonths[1]
  for (i in yearmonths){
    
    # read in distance matrix
    distmat = read_csv(paste(i, "-dist-matrix.csv", sep = ""),
                       trim_ws = TRUE, show_col_types = FALSE)
    
    # define high-distance accession
    distmat$Distance_Score <- vapply(2:ncol(distmat), FUN = 
                                       function(j){
                                         
                                         distances <- as.numeric(distmat[,j])
                                         distances <- distances[distances!=0]
                                         new_sum <- sum(distances)
                                         
                                         return(new_sum)
                                         
                                       }, numeric(1))
    distmat <- distmat[order(distmat$Sequence_Name),] ; rownames(distmat) <- NULL
    
    # collating cluster metadata with distances
    cluster_table <- read_tsv(paste(i, "-clusters.uc", sep = ""), 
                              col_names = FALSE, trim_ws = TRUE, show_col_types = FALSE)
    if (nrow(distmat)==1){
      cluster_table <- cluster_table[cluster_table$V9==distmat$Sequence_Name[1],]
    }
    centroids <- cluster_table[cluster_table$V1=="C",]
    centroids <- centroids[order(centroids$V9),] ; rownames(centroids) <- NULL
    stopifnot(nrow(distmat) == nrow(centroids))
    stopifnot(unique(distmat$Sequence_Name == centroids$V9) == TRUE)
    distmat$Cluster <- 0
    distmat$Cluster_Size <- centroids$V3
    distmat <- distmat[order(distmat$Distance_Score, decreasing = T),] ; rownames(distmat) <- NULL
    distmat$Month <- i 
    distmat$Cluster <- centroids$V2
    
    # new cluster table
    all_seqs <- distmat[, 
                        c("Sequence_Name", "Distance_Score", 
                          "Cluster_Size", "Month", "Cluster")]
    
    # adding all cluster sequences
    if ("H" %in% cluster_table$V1){
      hits <- cluster_table[cluster_table$V1 == "H",]
      for (j in 1:nrow(hits)){
        
        # record accession and the centroid for that accession
        accession <- hits$V9[j]
        centroid <- hits$V10[j]
        distance <- all_seqs[all_seqs$Sequence_Name==centroid, "Distance_Score"]
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
      distance <- as.numeric(all_seqs[all_seqs$Sequence_Name==centroid, "Distance_Score"])
      cluster_size <- as.numeric(all_seqs[all_seqs$Sequence_Name==centroid, "Cluster_Size"])
      month <- i
      cluster <- all_seqs[all_seqs$Sequence_Name==centroid, "Cluster"]
      
      # construct a new row for the all_seqs dataframe
      new_row <- c(accession, distance, cluster_size, month, cluster)
      
      # bind to table
      all_seqs <- rbind(all_seqs, new_row)
      
    }
    
    all_seqs <- all_seqs[order(all_seqs$Distance_Score),]
    rownames(all_seqs) <- NULL
    
    
    # bring these new data into the NCBI metadata
    for (j in 1:nrow(all_seqs)){
      
      accession <- all_seqs$Sequence_Name[j]
      
      if (!("Distance_Score" %in% colnames(metadata))){
        
        # create new columns to fill
        metadata$Distance_Score <- 0.0
        metadata$Cluster_Size <- 0
        metadata$Month_Cluster <- NA
        
      }
      
      # fill in distmat data in the NCBI metadata
      metadata[metadata$Accession==accession, 
               "Distance_Score"] <- as.numeric(all_seqs$Distance_Score[j])
      metadata[metadata$Accession==accession, 
               "Cluster_Size"] <- as.numeric(all_seqs$Cluster_Size[j])
      metadata[metadata$Accession==accession, 
               "Month_Cluster"] <- paste(all_seqs$Month[j], all_seqs$Cluster[j], sep = "_")
      
    }
    
  }
  
  return(metadata)
  
})

# create a function to normalize metadata
normalize_meta <- compiler::cmpfun(function(high_dist_seqs, high_dist_meta){
  
  high_dist_seqs <- high_dist_seqs[names(high_dist_seqs) %in% high_dist_meta$Accession]
  return(high_dist_seqs)
})

# candidate-identifying function to tie it all together
find_candidates <- compiler::cmpfun(function(metadata, stringency){
  
  # creating list of year-month combinations to loop through
  yearmonths <- str_remove_all(list.files(path = ".", pattern = "*-dist-matrix.csv"), "-dist-matrix.csv")
  
  # run the metadata function in a parallelized fashion
  metadata <- process_meta(yearmonths, metadata)
  
  # check to make sure metadata has clustering information in it
  stopifnot(nrow(metadata[!is.na(metadata$Month_Cluster),])>0)
  metadata <- metadata[metadata$Distance_Score>0,] ; rownames(metadata) <- NULL
  
  # make an empty data frame to store high distance FASTA sequences
  high_dist_seqs <- DNAStringSet()
  
  # bring in FASTAs
  fastas <- list.files(path = ".", pattern = "*-cluster-seqs*")
  for (fa in fastas){
    
    # bring in FASTA
    cluster_seqs <- readDNAStringSet(fa)
    
    # add that FASTA to a growing FASTA object of all hits, which will be exported
    # after normalizing later
    high_dist_seqs <- c(high_dist_seqs, cluster_seqs) ; remove(cluster_seqs)
    
  }
  
  # define retention threshold based on the data
  retention_threshold = quantile(metadata$Distance_Score, seq(0, 1, 0.001))[stringency]
  
  # plot distribution of summed weighted distances
  pdf("distance_distribution.pdf", width = 7, height = 5.5)
  hist_vector <-hist(metadata$Distance_Score, freq = FALSE, col = "lightblue",
                     xlab = "Distance Score", ylab = "Frequency", 
                     main = "Frequency Distribution of Nucleotide Distances")
  if (nrow(metadata) > 1){
    lines(density(metadata$Distance_Score), col = "darkblue", lwd = 2)
  }
  abline(v = retention_threshold, col = "red", lwd = 3)
  text(x = retention_threshold+5, y = (max(hist_vector$counts)/2), adj = 0,
       labels = paste("Retention Threshold:\n", retention_threshold ))
  dev.off()
  
  # normalizing year-month by down-prioritizing branches that are not exceptionally
  # long compared to all other year-months. I do this here by retaining only sequences that
  # are in clusters above the 90% quantile of distances, which pulls out true outliers.
  high_dist_meta <- subset(metadata, Distance_Score >= retention_threshold)
  
  # run the function
  high_dist_seqs <- normalize_meta(high_dist_seqs, high_dist_meta)
  
  # check here to prevent silent errors
  stopifnot(as.numeric(length(high_dist_seqs)) == nrow(high_dist_meta))
  
  # writing candidate FASTA of high-distance sequences that passed normalization and 
  # filtering
  writeXStringSet(high_dist_seqs, "high_distance_candidates.fasta")
  
  # writing candidate metadata, which combines vsearch --clust_fast results with 
  # NCBI metadata and iqTree branch lengths
  high_dist_meta <- high_dist_meta[order(high_dist_meta$Distance_Score, decreasing = T),]
  write.table(high_dist_meta, "high_distance_candidates.tsv", row.names = F, quote = F, sep = "\t")
  
})

# run the candidate-finder
find_candidates(metadata, stringency)

