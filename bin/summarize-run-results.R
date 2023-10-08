#!/usr/bin/env Rscript

options(warn = 1)
suppressPackageStartupMessages({

  require(compiler)
  require(argparse)
  require(tidyr)
  require(stringr)
  require(dplyr)
  require(ggplot2)
  require(readr)

})

parse_args <- compiler::cmpfun(function() {
  parser <- ArgumentParser(description = "R script that generates a table of \
                           ALPINE summary statistics.")
  parser$add_argument("--results", default = "./results",
                      help = "Parent directory for results.")
  args <- parser$parse_args()

  return(args)
})

construct_file_paths <- compiler::cmpfun(function(search_dir) {
  # Construct named vector where each item is a results parent directory and the
  # name is the geography

  paths <- list.dirs(path = search_dir,
                     full.names = TRUE, recursive = FALSE)

  stopifnot(length(paths) > 0)

  geographies <- list.dirs(path = search_dir,
                           full.names = FALSE, recursive = FALSE) %>%
    str_split_i("_", i = -1) %>%
    str_replace_all("_", " ")

  names(paths) <- geographies

  return(paths)

})

define_search_tree <- compiler::cmpfun(function(path) {

  subdirs <- list.dirs(path)
  query_pattern <- "double_candidates|high_distance_clusters|metadata_candidates"
  subdirs_to_search <- subdirs[grepl(query_pattern, subdirs)]

  if (length(subdirs_to_search) == 0) {
    return(list())
  }

  search_tree <- list(subdirs[1])
  names(search_tree) <- "parent"
  doubles <- grepl("double_candidates", subdirs_to_search)
  if (TRUE %in% doubles) {
    search_tree$double_candidates <- subdirs_to_search[doubles]
  }
  anachrons <- grepl("metadata_candidates", subdirs_to_search)
  if (TRUE %in% anachrons) {
    search_tree$anachronistics <- subdirs_to_search[anachrons]
  }
  high_dists <- grepl("high_distance_clusters", subdirs_to_search)
  if (TRUE %in% high_dists) {
    search_tree$high_distance <- subdirs_to_search[high_dists]
  }

  return(search_tree)

})

get_early_counts <- compiler::cmpfun(function(path) {

  if (!(paste(path, "early_stats.tsv", sep = "") %in% list.files(path))){
    return(0)
  }

  early_stats <- read_tsv(paste(path, "early_stats.tsv", sep = ""))
  count <- early_stats$num_seqs[1]
  return(count)

})

get_late_counts <- compiler::cmpfun(function(path) {

  if (!(paste(path, "late_stats.tsv", sep = "") %in% list.files(path))){
    return(0)
  }

  late_stats <- read_tsv(paste(path, "late_stats.tsv", sep = ""))
  count <- late_stats$num_seqs[1]
  return(count)

})

summarize_doubles <- compiler::cmpfun(function(template_df) {

  with_double_stats <- template_df %>%
    mutate(`Double Candidate Prevalence (%)` =
             `Double Candidate Count` / `Input Sequence Count`) %>%
    mutate(`Double Candidate Rate` =
             paste("1 in ", floor(1 / `Double Candidate Prevalence (%)`)),
           sep = "")

  return(with_double_stats)

})

summarize_anachrons <- compiler::cmpfun(function(template_df, search_tree) {

  template_df$`Anachronistic Count` <- as.integer(NA)
  template_df$`Anachronistic Prevalence (%)` <- as.numeric(NA)
  template_df$`Anachronistic Rate` <- as.character(NA)

  for (geo in names(search_tree)) {

    index <- which(template_df$Geography == geo)

    if (length(search_tree[[geo]]) == 0) {
      template_df$`Anachronistic Count`[index] <- NA
    }
    if (!("anachronistics" %in% names(search_tree[[geo]]))) {
      template_df$`Anachronistic Count`[index] <- NA
    }

    anachron_dir <- search_tree[[geo]][["anachronistics"]]
    anachron_df <- read_tsv(paste(anachron_dir, 
                                  "anachronistic_metadata_only_candidates.tsv",
                                  sep = "/"))

    template_df$`Anachronistic Count`[index] <- nrow(anachron_df)
    template_df$`Anachronistic Prevalence (%)`[index] <- (
      nrow(anachron_df) / template_df$`Input Sequence Count`[index]
    ) * 100
    template_df$`Anachronistic Rate`[index] <- paste(
      "1 in ",
      floor(
        1 / template_df$`Anachronistic Prevalence (%)`[index]
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

    index <- which(template_df$Geography == geo)

    if (length(search_tree[[geo]]) == 0){
      template_df$`Anachronistic Count`[index] <- NA
    }
    if (!("high_distance_clusters" %in% names(search_tree[[geo]]))) {
      template_df$`Anachronistic Count`[index] <- NA
    }

    high_dist_dir <- search_tree[[geo]][["high_distance_clusters"]]
    high_dist_df <- read_tsv(paste(anachron_dir, 
                                   "high_distance_candidates.tsv",
                                   sep = "/"))

    template_df$`High Distance Count`[index] <- nrow(high_dist_df)
    template_df$`High Distance Prevalence (%)`[index] <- (
      nrow(high_dist_df) / template_df$`Input Sequence Count`[index]
    ) * 100
    template_df$`High Distance Rate`[index] <- paste(
      "1 in ",
      floor(
        1 / template_df$`High Distance Prevalence (%)`[index]
      ),
      sep = "")

  }

  return(template_df)

})

main <- compiler::cmpfun(function() {

  # determine the file tree to be traversed in search of results
  paths <- construct_file_paths(parse_args()$results)

  # search the tree and construct a nested named list (like an R version of a
  # dictionary of dictionaries) that will guide downstream steps
  search_tree <- lapply(paths, define_search_tree)
  names(search_tree) <- names(paths)

  # use the search tree to construct dataframes of anachronistic, high-distance,
  # and double candidate results, and then joirn them
  parentdirs <- sapply(search_tree, function(x) x[["parent"]])
  run_results <- tibble("Geography" = names(search_tree),
                        "Input Sequence Count" = lapply(
                          parentdirs,
                          get_early_counts),
                        "Double Candidate Count" = lapply(
                          parentdirs,
                          get_late_counts
                        )) %>%
    summarize_doubles() %>%
    summarize_anachrons(search_tree) %>%
    summarize_highdist(search_tree)

  # Write the new results and create a few visualizations
  write_tsv(run_results, "alpine_run_summary.tsv", na = "")

})
main()
