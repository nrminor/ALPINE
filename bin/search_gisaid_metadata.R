#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(parallel)

# reading in filtered metadata and pango lineage dates
metadata <- read.delim(args[1], na.strings = c("NA", ""))
dates <- read.csv(args[2])
infection_cutoff <- as.numeric(args[3])
cores <- as.numeric(args[4])

# cutting out any rows that don't have pango classifications from GISAID  
metadata <- metadata[!is.na(metadata$pango),]
metadata <- metadata[metadata$pango!="Unassigned",]

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
gisaid_dates <- metadata$date
desig_dates <- metadata$lineage_designation
clusterExport(cluster, varlist=c("gisaid_dates", "desig_dates"))
metadata$infection_duration <- do.call(c, 
                                        clusterApplyLB(cluster, 1:length(gisaid_dates), function(i){
                                          return(as.numeric(gisaid_dates[i] - desig_dates[i]))
                                        }))
metadata <- metadata[!is.na(metadata$lineage_designation),]
metadata <- metadata[!is.na(metadata$infection_duration),]

# filtering down to long infection candidates
long_infections <- metadata[metadata$infection_duration >= infection_cutoff &
                              metadata$date > as.Date("2021-02-18"),]
long_infections <- long_infections[order(long_infections$infection_duration, decreasing = T),]
long_infections <- long_infections[!is.na(long_infections$accession),]
rownames(long_infections) <- NULL

# exporting TSV of long infection candidates
write.table(long_infections,
          paste("long_infections_gisaid_metadata_",
                Sys.Date(), 
                ".tsv", sep = ""),
          row.names = F, quote = F, na = "", sep = "\t")
