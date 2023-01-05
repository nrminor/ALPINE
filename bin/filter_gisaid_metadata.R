#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# collecting inputs/parameters
full_list <- read.delim(args[1], quote = "", 
                        row.names = NULL, 
                        stringsAsFactors = FALSE,
                        na.strings = c("NA", ""))
min_date <- as.Date(args[2])
max_date <- as.Date(args[3])
geography <- args[4]

# filtering
full_list$Collection.date <- as.Date(full_list$Collection.date)
full_list <- full_list[!is.na(full_list$Collection.date),]
full_list <- full_list[full_list$Collection.date >= min_date &
                         full_list$Collection.date <= max_date,]
full_list <- full_list[!is.na(full_list$Location),]
full_list <- full_list[grepl(geography,full_list$Location),]
# full_list <- full_list[grepl("usa",tolower(full_list$Location)),]

# cleaning up pango lineages
full_list <- full_list[!is.na(full_list$Pango.lineage),]
full_list <- full_list[full_list$Pango.lineage!="Unassigned",]

# formatting and exporting include list
include_list <- full_list[, c("Virus.name", "Accession.ID", "Collection.date", "Location", "Pango.lineage")]
rownames(include_list) <- NULL
colnames(include_list) <- c("strain", "accession", "date", "location", "pango")

write.csv(include_list, file = "include_list.csv", 
          quote = F, row.names = F)