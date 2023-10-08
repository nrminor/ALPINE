#!/usr/bin/env Rscript

# retrieve positional command line arguments
args <- commandArgs(TRUE)

options(warn = 1)
suppressPackageStartupMessages({

  require(compiler)
  require(tidyr)
  require(stringr)
  require(dplyr)
  require(ggplot2)
  require(readr)

})

construct_file_paths <- compiler::cmpfun(function(search_dir) {
  
  paths <- list.dirs(path = search_dir,
                     full.names = TRUE, recursive = FALSE)
  
  stopifnot(length(paths) > 0)
  
  geographies <- list.dirs(path = search_dir,
                           full.names = FALSE, recursive = FALSE) %>%
    str_remove("LocalDataset_") %>%
    str_remove("GenBank_") %>%
    str_replace_all("_", " ")
  
  names(paths) <- geographies
  
  return(paths)
  
})

define_search_tree <- compiler::cmpfun(function(path) {
  
  subdirs <- list.dirs(path)
  query_pattern <- "double_candidates|high_distance_clusters|metadata_candidates"
  subdirs_to_search <- subdirs[grepl(query_pattern, subdirs)]
  
  if (length(subdirs_to_search) == 0 ) {
    return(list())
  }
  
  search_tree <- list(subdirs[1])
  names(search_tree) <- "parent"
  doubles <- grepl("double_candidates", subdirs_to_search)
  if (TRUE %in% doubles){
    search_tree$double_candidates <- subdirs_to_search[doubles]
  }
  anachrons <- grepl("metadata_candidates", subdirs_to_search)
  if (TRUE %in% anachrons){
    search_tree$anachronistics <- subdirs_to_search[anachrons]
  }
  high_dists <- grepl("high_distance_clusters", subdirs_to_search)
  if (TRUE %in% high_dists){
    search_tree$high_distance <- subdirs_to_search[high_dists]
  }
  
  return(search_tree)
  
})

get_early_counts <- compiler::cmpfun(function(path) {
  
  names(path) <- NULL
  
  if (!("early_stats.tsv" %in% list.files(path))) {
    return(0)
  }
  
  early_stats <- read_tsv(paste(path, "early_stats.tsv", sep = "/"),
                          show_col_types = FALSE,  trim_ws = TRUE)
  count <- early_stats$num_seqs[1]
  return(count)
  
})

get_late_counts <- compiler::cmpfun(function(path) {
  
  names(path) <- NULL
  
  if (!("late_stats.tsv" %in% list.files(path))){
    return(0)
  }
  
  late_stats <- read_tsv(paste(path, "late_stats.tsv", sep = "/"),
                         show_col_types = FALSE,  trim_ws = TRUE)
  count <- late_stats$num_seqs[1]
  return(count)
  
})

summarize_doubles <- compiler::cmpfun(function(template_df) {
  
  with_double_stats <- template_df %>%
    mutate(`Double Candidate Prevalence (%)` = 
             (`Double Candidate Count` / `Input Sequence Count`) * 100) %>%
    mutate(`Double Candidate Rate` = if_else(
      `Double Candidate Prevalence (%)` == 0, 
      "None",
      paste("1 in ", floor(1 / (`Double Candidate Prevalence (%)` / 100)))
    ))
  
  return(with_double_stats)
  
})

summarize_anachrons <- compiler::cmpfun(function(template_df, search_tree) {
  
  template_df$`Anachronistic Count` <- as.integer(NA)
  template_df$`Anachronistic Prevalence (%)` <- as.numeric(NA)
  template_df$`Anachronistic Rate` <- as.character(NA)
  
  for (geo in names(search_tree)) {
    
    index <- which(grepl(geo, template_df$Geography))
    
    if (length(search_tree[[geo]]) == 0){
      template_df$`Anachronistic Count`[index] <- NA
    }
    if ( !("anachronistics" %in% names(search_tree[[geo]])) ) {
      template_df$`Anachronistic Count`[index] <- NA
    }
    
    anachron_dir <- search_tree[[geo]][["anachronistics"]]
    anachron_df <- read_tsv(paste(anachron_dir, 
                                  "anachronistic_metadata_only_candidates.tsv",
                                  sep = "/"),
                            show_col_types = FALSE,  trim_ws = TRUE)
    
    template_df$`Anachronistic Count`[index] <- nrow(anachron_df)
    template_df$`Anachronistic Prevalence (%)`[index] <- (
      nrow(anachron_df) / template_df$`Input Sequence Count`[index]
    ) * 100
    template_df$`Anachronistic Rate`[index] <- paste(
      "1 in ",
      floor(
        1 / (template_df$`Anachronistic Prevalence (%)`[index] / 100)
      ),
      sep = "")
    
  }
  
  return(template_df)
  
})

summarize_highdist <- compiler::cmpfun(function(template_df, search_tree) {
  
  template_df$`High Distance Count` <- as.integer(NA)
  template_df$`High Distance Prevalence (%)` <- as.numeric(NA)
  template_df$`High Distance Rate` <- as.character(NA)
  
  for (geo in names(search_tree)) {
    
    index <- which(grepl(geo, template_df$Geography))
    
    if (length(search_tree[[geo]]) == 0){
      template_df$`Anachronistic Count`[index] <- NA
    }
    if ( !("high_distance_clusters" %in% names(search_tree[[geo]])) ) {
      template_df$`Anachronistic Count`[index] <- NA
    }
    
    high_dist_dir <- search_tree[[geo]][["high_distance"]][1]
    high_dist_df <- read_tsv(paste(high_dist_dir, 
                                   "high_distance_candidates.tsv",
                                   sep = "/"),
                             show_col_types = FALSE,  trim_ws = TRUE)
    
    template_df$`High Distance Count`[index] <- nrow(high_dist_df)
    template_df$`High Distance Prevalence (%)`[index] <- (
      nrow(high_dist_df) / template_df$`Input Sequence Count`[index]
    ) * 100
    template_df$`High Distance Rate`[index] <- paste(
      "1 in ",
      floor(
        1 / (template_df$`High Distance Prevalence (%)`[index] / 100)
      ),
      sep = "")
    
  }
  
  return(template_df)
  
})

main <- compiler::cmpfun(function() {

  # determine the file tree to be traversed in search of results
  paths <- construct_file_paths(args[1])

  # search the tree and construct a nested named list (like an R version of a
  # dictionary of dictionaries) that will guide downstream steps
  search_tree <- lapply(paths, define_search_tree)
  names(search_tree) <- names(paths)

  # use the search tree to construct dataframes of anachronistic, high-distance,
  # and double candidate results, and then joirn them
  parentdirs <- sapply(search_tree, function(x) x[["parent"]])
  run_results <- tibble("Geography" = names(search_tree),
                        "Input Sequence Count" = vapply(
                          paths,
                          get_early_counts,
                          double(1)),
                        "Double Candidate Count" = vapply(
                          paths,
                          get_late_counts,
                          double(1))) %>%
    summarize_doubles() %>%
    summarize_anachrons(search_tree) %>%
    summarize_highdist(search_tree)

  # Write the new results and create a few visualizations
  write_tsv(run_results, "alpine_run_summary.tsv", na = "")

})
main()
