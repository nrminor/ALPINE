#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# collecting inputs/parameters
full_list <- read.delim(args[1], quote = "", 
                        row.names = NULL, 
                        stringsAsFactors = FALSE,
                        na.strings = c("NA", ""))
if (!is.na(args[2])){
  min_date <- as.Date(args[2])
} else {
  min_date <- as.Date("2019-10-31")
}
if (!is.na(args[3])){
  max_date <- as.Date(args[3])
} else {
  max_date <- as.Date(Sys.Date())
}
if (!is.na(args[4])){
  geography <- args[4]
}

# bringing in US State names/territories & abbreviations
states = c("Alabama", "Alaska", "American Samoa", "Arizona", "Arkansas",
"California", "Colorado", "Connecticut", "Delaware", "District of Columbia",
"Florida", "Georgia", "Guam", "Hawaii", "Idaho", "Illinois", "Indiana",
"Iowa", "Kansas", "Kentucky", "Louisiana", "Maine", "Maryland",
"Massachusetts", "Michigan", "Minnesota", "Minor Outlying Islands",
"Mississippi", "Missouri", "Montana", "Nebraska", "Nevada", "New Hampshire",
"New Jersey", "New Mexico", "New York", "North Carolina", "North Dakota",
"Northern Mariana Islands", "Ohio", "Oklahoma", "Oregon", "Pennsylvania",
"Puerto Rico", "Rhode Island", "South Carolina", "South Dakota", "Tennessee",
"Texas", "U.S. Virgin Islands", "Utah", "Vermont", "Virginia", "Washington",
"West Virginia", "Wisconsin", "Wyoming") 
states_abbrev <- c("AK", "AL", "AR", "AS", "AZ", "CA", "CO", "CT", "DC", "DE",
"FL", "GA", "GU", "HI", "IA", "ID", "IL", "IN", "KS", "KY", "LA", "MA", "MD",
"ME", "MI", "MN", "MO", "MP", "MS", "MT", "NC", "ND", "NE", "NH", "NJ", "NM",
"NV", "NY", "OH", "OK", "OR", "PA", "PR", "RI", "SC", "SD", "TN", "TX", "UM",
"UT", "VA", "VI", "VT", "WA", "WI", "WV", "WY")

# fixing improperly formatted dates
full_list <- full_list[!is.na(full_list$Isolate.Collection.date),]
full_list$Isolate.Collection.date <- as.Date(full_list$Isolate.Collection.date)
full_list <- full_list[!is.na(full_list$Isolate.Collection.date),]

# filtering
if (!is.na(args[2])){
  full_list <- full_list[full_list$Isolate.Collection.date >= min_date,]
}
if (!is.na(args[3])){
  full_list <- full_list[full_list$Isolate.Collection.date <= max_date,]
}
if(!is.na(args[4])){
  if (geography == "USA"){
    usa <- paste(geography, ": ", sep = "")
    full_list <- full_list[grepl(usa,full_list$Geographic.Location),]
  } else if ( geography %in% states | geography %in% states_abbrev ) {
    abbreviated <- paste("USA: ", toupper(substr(geography,1,2)), sep = "")
    full_list <- full_list[grepl(geography,full_list$Geographic.Location) | 
                             grepl(abbreviated,full_list$Geographic.Location),]
  } else {
    full_list <- full_list[grepl(tolower(geography),tolower(full_list$Geographic.Location)),]
  }
}

# cleaning up pango lineages
full_list <- full_list[!is.na(full_list$Virus.Pangolin.Classification),]
full_list <- full_list[full_list$Virus.Pangolin.Classification!="unclassifiable",]
rownames(full_list) <- NULL

# formatting and exporting include list
include_list <- full_list[, c("Accession", "Isolate.Collection.date", "Geographic.Location", "Virus.Pangolin.Classification")]
rownames(include_list) <- NULL ; remove(full_list)
colnames(include_list) <- c("accession", "date", "location", "pango")

write.table(include_list, file = "include_list.tsv", 
          sep = "\t", quote = F, row.names = F)