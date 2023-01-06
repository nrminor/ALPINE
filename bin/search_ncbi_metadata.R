#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(parallel)

# reading in filtered metadata and pango lineage dates
metadata <- read.csv(args[1], na.strings = c("NA", ""))
dates <- read.csv(args[2])
infection_cutoff <- as.numeric(args[3])
cores <- as.numeric(args[4])

# cutting out any rows that don't have pango classifications from NCBI
metadata <- metadata[!is.na(metadata$pango),]

# ensuring dates are properly formatted
metadata$date <- as.Date(metadata$date)
dates$designation_date <- as.Date(dates$designation_date)
dates[is.na(dates$designation_date),2] <- as.Date("2021-02-18")

# parallelize looping through filtered metadata to merge designation dates
metadata$lineage_designation <- as.Date(NA)
metadata$infection_duration <- 0
cluster <- makeCluster(cores)
pango <- metadata$pango
clusterExport(cluster, varlist=c("pango", "dates"))
metadata$lineage_designation <- do.call(c, 
                                        clusterApplyLB(cluster, pango, function(i){
                                          return(as.Date(dates[dates$lineage==i, "designation_date"]))
                                        }))
ncbi_dates <- metadata$date
desig_dates <- metadata$lineage_designation
clusterExport(cluster, varlist=c("ncbi_dates", "desig_dates"))
metadata$infection_duration <- do.call(c, 
                                       clusterApplyLB(cluster, 1:length(ncbi_dates), function(i){
                                         return(as.numeric(ncbi_dates[i] - desig_dates[i]))
                                       }))
metadata <- metadata[!is.na(metadata$lineage_designation),]

# filtering down to long infection candidates
long_infections <- metadata[metadata$infection_duration >= infection_cutoff &
                              metadata$date > as.Date("2021-02-18"),]
rownames(long_infections) <- NULL
long_infections <- long_infections[!is.na(long_infections$accession),]

# exporting CSV of long infection candidates
write.csv(long_infections,
          paste("long_infections_ncbi_metadata_",
                Sys.Date(), 
                ".csv", sep = ""),
          row.names = F, quote = F, na = "")
