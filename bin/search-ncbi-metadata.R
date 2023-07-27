#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(tidyverse)
library(parallel)

# set just-in-time compilation setting so that top-level loops are JIT compiled
compiler::enableJIT(3)

# reading in filtered metadata and pango lineage dates
metadata <- read_tsv(args[1], show_col_types = FALSE, trim_ws = TRUE)
dates <- read.csv(args[2])
infection_cutoff <- as.numeric(args[3])
cores <- as.numeric(args[4])

# cutting out any rows that don't have pango classifications from NCBI
metadata <- metadata %>%
  filter(!is.na(`Virus Pangolin Classification`))

# ensuring dates are properly formatted
metadata$`Isolate Collection date` <- as.Date(metadata$`Isolate Collection date`)
dates$designation_date <- as.Date(dates$designation_date)
dates[is.na(dates$designation_date),2] <- as.Date("2021-02-18")

# parallelize looping through filtered metadata to merge designation dates
metadata$lineage_designation <- as.Date(NA)
metadata$infection_duration <- 0
cluster <- makeCluster(cores)
pango <- metadata$`Virus Pangolin Classification`
clusterExport(cluster, varlist=c("pango", "dates"))
metadata$lineage_designation <- do.call(c, 
                                        clusterApplyLB(cluster, pango, function(i){
                                          return(as.Date(dates[dates$lineage==i, "designation_date"]))
                                        }))
ncbi_dates <- metadata$`Isolate Collection date`
desig_dates <- metadata$lineage_designation
clusterExport(cluster, varlist=c("ncbi_dates", "desig_dates"))
metadata$infection_duration <- do.call(c, 
                                       clusterApplyLB(cluster, 1:length(ncbi_dates), function(i){
                                         return(as.numeric(ncbi_dates[i] - desig_dates[i]))
                                       }))
metadata <- filter(metadata,!is.na(lineage_designation))

# filtering down to long infection candidates
long_infections <- metadata %>%
  filter(infection_duration >= infection_cutoff &
         lineage_designation > as.Date("2021-02-18")) %>%
  arrange(desc(infection_duration)) %>%
  filter(!is.na(Accession))

# exporting CSV of long infection candidates
write.table(long_infections,
            paste("anachronistic_metadata_only_candidates.tsv", sep = ""),
            row.names = F, quote = F, na = "", sep = "\t")
