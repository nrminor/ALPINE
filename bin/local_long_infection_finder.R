#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# bringing in sequencing run metadata
fasta_basename <- args[1]
experiment_name <- args[2]
experiment_date <- as.Date(args[3]) ; experiment_date

# reading in pango lineages and designation dates
lineage_csv <- read.csv(args[4])
dates <- read.csv(args[5])

# defining number of days past designation to consider an infection prolonged
days_infected = as.numeric(args[6])

# preparing a data frame to hold long infection data, if detected
long_infections <- data.frame(sample = NA,
                              lineage = NA,
                              lineage_designation_date = NA,
                              experiment_date = NA,
                              experiment_name = NA,
                              pango_version = NA)

# Making sure the lineage csv isn't empty
if (nrow(lineage_csv) != 0){
  
  # ensuring dates are properly formatted
  dates$designation_date <- as.Date(dates$designation_date)
  
  # diff <- c()
  for (i in 1:nrow(lineage_csv)){
    
    lineage <- lineage_csv$lineage[i]
    designation_date <- dates[dates$lineage==lineage, 2]
    
    if ( is.na(designation_date) ) {
      designation_date <- as.Date("2021-02-18")
    }
    
    # diff <- c(diff, as.numeric(experiment_date - designation_date))
    # mean(diff) ; hist(diff, breaks = 20)
    
    if ( as.numeric(experiment_date - designation_date) >= days_infected ){
      
      new_row <- c(lineage_csv$taxon[i],
                   lineage,
                   as.character(designation_date),
                   as.character(experiment_date),
                   experiment_name,
                   lineage_csv$pangolin_version[i])
      long_infections <- rbind(long_infections, new_row)
      
    } else {
      next
    }
    
    
  }
  if ( nrow(long_infections) > 1){
    long_infections <- long_infections[2:nrow(long_infections),]
    rownames(long_infections) <- NULL
  }
  
}

if (nrow(long_infections)==1 &&
    is.na(long_infections[1, "sample"])){
  
  long_infections[1,1] <- paste("No putative long infections were identified in experiment",
                       experiment_name, "on", Sys.Date(), sep = " ")
  
}

write.csv(long_infections,
          paste(fasta_basename,
                "_putative_long_infections_",
                Sys.Date(), 
                ".csv", sep = ""),
          row.names = F, quote = F, na = "")
