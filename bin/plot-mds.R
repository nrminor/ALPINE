#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# load packages
library(tidyverse)
library(ape)
library(ggplot2)
library(readr)

# Read command-line arguments, which will serve as the script's only globals
yearmonth <- args[1]
tsv_path <- args[2]
fasta_path <- args[3]

# define a function that runs and plots the MDS
run_mds <- compiler::cmpfun(function(dist_mat, tsv_path, yearmonth, num_dims){
  
  # Perform MDS
  mds_result <- stats::cmdscale(dist_mat, k = num_dims) %>% 
    as.data.frame()
  mds_result <- cbind(row.names(mds_result), mds_result) ; rownames(mds_result) <- NULL
  if(num_dims == 3){
    colnames(mds_result) <- c("Accession", "Dim1", "Dim2", "Dim3")
  } else {
    colnames(mds_result) <- c("Accession", "Dim1", "Dim2")
  }                
  
  # Read TSV table
  tsv_data <- readr::read_tsv(tsv_path,col_names=F,show_col_types = FALSE, trim_ws = TRUE)
  tsv_data <- tsv_data[tsv_data$X1=="C",]
  no_seqs <- nrow(tsv_data)
  colnames(tsv_data)[2] <- "Cluster"
  tsv_data$Cluster <- as.factor(tsv_data$Cluster)
  colnames(tsv_data)[3] <- "Size"
  tsv_data$Size <- as.numeric(tsv_data$Size)
  colnames(tsv_data)[9] <- "Accession"
  
  # Merge MDS results with TSV table
  merged_data <- merge(mds_result, tsv_data, by = "Accession", all.x = TRUE)
  
  # Plot MDS with colored points based on the table
  ggplot(merged_data, aes(x = Dim1, y = Dim2, color = Cluster, size = Size)) +
    geom_point(alpha=0.6) +
    geom_text(aes(label = Size), vjust = 2.5, size = 3, show.legend = FALSE) +
    labs(title = paste("Multidimensional Scaling of", no_seqs, "clusters in", yearmonth, sep = " "), 
         x = "Dimension 1", y = "Dimension 2") +
    theme(text = element_text(size = 12),
          legend.position = "none") +
    scale_size(range = c(1, 35))
  
  # Save the plot as a PDF
  ggsave(paste(yearmonth, "mds_plot.pdf", sep = "_"), width = 7, height = 5)
  
  if (num_dims == 3){
    
    # plotly code for interactive 3D-plotting
    library(plotly)
    library(htmlwidgets)
    p <- plot_ly(data = merged_data, x = ~Dim1, y = ~Dim2, z = ~Dim3,
                 color = ~Cluster, size = ~Size,
                 type = "scatter3d", mode = "markers") %>%
      add_markers(text = ~Accession, hoverinfo = "text",
                  marker = list(sizemode = "diameter", sizeref = 1)) %>%
      layout(title = paste("Multidimensional Scaling of", no_seqs, "clusters in", yearmonth, sep = " "),
             scene = list(
               xaxis = list(title = "Dimension 1"),
               yaxis = list(title = "Dimension 2"),
               zaxis = list(title = "Dimension 3")
             ),
             showlegend = FALSE,
             hoverlabel = list(bgcolor = "white",
                               font = list(size = 12)))
    saveWidget(p, file = paste(yearmonth, "3D_interactive.html", sep = "_"))
    
    # Save MDS results as a CSV
    write.csv(merged_data, paste(yearmonth, "mds_results.csv", sep = "_"),
              quote = F, row.names = F)
    
  }
  
})

main <- compiler::cmpfun(function(fasta_path, tsv_path, yearmonth){
  
  # Read FASTA file
  fasta_data <- ape::read.dna(fasta_path, format = "fasta")
  
  if (length(labels(fasta_data)) > 1){
    
    # Calculate number of nucleotide differences
    dist_matrix <- ape::dist.dna(fasta_data, model = "N")
    
    # Perform MDS
    if (length(labels(fasta_data)) > 3){
      run_mds(dist_matrix, tsv_path, yearmonth, 3)
    } else {
      run_mds(dist_matrix, tsv_path, yearmonth, 2)
    }
    
    
  } else {
    
    # Read TSV table
    tsv_data <- readr::read_tsv(tsv_path,col_names=F)
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
  
})

# run the code
main(fasta_path, tsv_path, yearmonth)
