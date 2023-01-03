#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# reading in pango lineages, filtered GISAID metadata, and designation dates
lineage_csv <- read.csv(args[1])
metadata <- read.csv(args[2])
dates <- read.csv(args[3])

# defining number of days past designation to consider an infection prolonged
days_infected = as.numeric(args[4])

# preparing a data frame to hold long infection data, if detected
long_infections <- data.frame(sample = NA,
                              lineage = NA,
                              lineage_designation_date = NA,
                              collection_date = NA,
                              repository = NA,
                              pango_version = NA)