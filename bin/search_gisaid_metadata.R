#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# reading in filtered metadata and pango lineage dates
metadata <- read.csv(args[1], na.strings = c("NA", ""))
dates <- read.csv(args[2])
infection_cutoff <- as.numeric(args[3])

# cutting out any rows that don't have pango classifications from NCBI
metadata <- metadata[!is.na(metadata$pango),]
metadata <- metadata[metadata$pango!="",]

# ensuring dates are properly formatted
metadata$date <- as.Date(metadata$date)
dates$designation_date <- as.Date(dates$designation_date)
dates[is.na(dates$designation_date),2] <- as.Date("2021-02-18")

# looping through filtered metadata to merge designation dates
metadata$lineage_designation <- as.Date(NA)
metadata$infection_duration <- 0
for (i in 1:nrow(metadata)){
  
  metadata$lineage_designation[i] <- as.Date(dates[dates$lineage==metadata$pango[i], "designation_date"])
  metadata$infection_duration[i] <- as.numeric(metadata$date[i] - metadata$lineage_designation[i])
  
}

# filtering down to long infection candidates
long_infections <- metadata[metadata$infection_duration >= infection_cutoff,]
rownames(long_infections) <- NULL

# exporting CSV of long infection candidates
write.csv(long_infections,
          paste("long_infections_gisaid_metadata_",
                Sys.Date(), 
                ".csv", sep = ""),
          row.names = F, quote = F, na = "")
