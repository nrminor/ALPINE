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

# fixing improperly formatted dates
full_list <- full_list[!is.na(full_list$Isolate.Collection.date),]
full_list$Isolate.Collection.date <- as.Date(full_list$Isolate.Collection.date)
full_list <- full_list[!is.na(full_list$Isolate.Collection.date),]

# filtering
full_list <- full_list[full_list$Isolate.Collection.date >= min_date &
                         full_list$Isolate.Collection.date <= max_date,]
if (geography == "USA"){
  usa <- paste(geography, ": ", sep = "")
  full_list <- full_list[grepl(usa,full_list$Geographic.Location),]
} else {
  abbreviated <- paste("USA: ", toupper(substr(geography,1,2)), sep = "")
  full_list <- full_list[grepl(geography,full_list$Geographic.Location) | 
                           grepl(abbreviated,full_list$Geographic.Location),]
}

# cleaning up pango lineages
full_list <- full_list[!is.na(full_list$Virus.Pangolin.Classification),]
full_list <- full_list[full_list$Virus.Pangolin.Classification!="unclassifiable",]

# formatting and exporting include list
include_list <- full_list[, c("Accession", "Isolate.Collection.date", "Geographic.Location", "Virus.Pangolin.Classification")]
colnames(include_list) <- c("accession", "date", "location", "pango")

write.table(include_list, file = "include_list.tsv", 
          sep = "\t", quote = F, quote = F, row.names = F)