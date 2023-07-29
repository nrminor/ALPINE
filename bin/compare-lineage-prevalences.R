#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# loading necessary libraries
require(parallel)
require(compiler)
require(tidyverse)
require(Biostrings)
require(cdcfluview)
require(outbreakinfo)

# set just-in-time compilation setting so that top-level loops are JIT compiled
compiler::enableJIT(3)

# Record pangolin report csv path, fasta path and metadata path
report_path <- args[1]
fasta_path <- args[2]
metadata_path <- args[3]

# launch authentication prompt to get GISAID token
# authenticateUser()

# printing GISAID Disclaimer
print("This data was obtained from GISAID via the outbreak.info API.")

# function that uses dplyr to arrange metadata
date_pango_calls <- cmpfun(function(report_path, metadata){
  
  # Read in pangolin report columns
  pango_report <- read_csv(report_path, col_select = c(`taxon`, lineage),
                           show_col_types = FALSE, trim_ws = TRUE) %>%
    filter(!grepl("taxon", `taxon`))
  
  pango_report$taxon <- str_split_i(pango_report$taxon,
                                    " ",
                                    i = 1)
  
  # joining the desired columns
  pango_with_dates <- pango_report %>% 
    left_join(metadata %>% select(Accession, `Isolate Collection date`), 
              by = c("taxon" = "Accession"), keep = TRUE) %>%
    select(-Accession) %>%
    arrange(`taxon`)
  
  return(pango_with_dates)
  
})

# function to create lookup of "rare dates" for the pango lineages in the 
# current pango report
create_rarity_lookup <- cmpfun(function(pango_report){
  
  # create vector of unique lineages to iterate through
  lineages <- unique(pango_report$lineage)
  
  # create rarity lookup vector with vapply
  pdf("lineage_prevalences.pdf", width = 8, height = 6)
  rarity_lookup <- vapply(lineages, function(lineage){
    
    # creating an empty table where we will fill in
    # prevalences
    prevalence_table <- getPrevalence(pangolin_lineage = lineage, 
                                      location = "United States")
    
    # Then, we get rid of zeroes, which will
    # skew our methods below. We also get
    # rid of 1.0's, which are likely
    # erroneous.
    prevalence_table <- prevalence_table[prevalence_table$proportion > 0 &
                                           prevalence_table$proportion < 1,]
    
    if (nrow(prevalence_table) > 0){
      
      # Next, we use a rather crude method to identify
      # when the wave for a lineage in question has
      # "crested". This is the date of peak prevalence
      # for the lineage in the HHS region.
      crest_date <- max(prevalence_table[
        prevalence_table$proportion==max(prevalence_table$proportion), 
        "date"])
      
      # We use the same logic to compute when a given
      # lineage has become rare, which we consider to
      # be in the lowest 10% quantile of prevalence
      # throughout its spread in the HHS, which will
      # be out cutoff for identifying prolonged
      # infection candidates. Importantly, it is
      # possible for a lineage to experience ups and
      # downs in its prevalence. To guarantee that we
      # don't erroneously call a sample anachronistic
      # based on a temporary prevalence low point, we
      # use the *last* date at which a lineage is in
      # its lowest 10% of prevalence in the HHS region
      # in question.
      q <- quantile(prevalence_table$proportion, 0.05)
      rare_prev <- prevalence_table$proportion[
        which.min(
          abs(prevalence_table$proportion - q))]
      rare_date <- max(prevalence_table[
        prevalence_table$proportion==rare_prev, "date"])
      
      if (rare_date <= crest_date) {
        
        post_crest <- prevalence_table %>%
          filter(date > crest_date)
        
        rare_prev <- post_crest$proportion[
          which.min(
            abs(post_crest$proportion - q))]
        rare_date <- max(post_crest[
          post_crest$proportion==rare_prev, "date"])
        
      }
      
      # plotting code for testing purposes
      plot(prevalence_table$proportion~prevalence_table$date, 
           pch = 20, col = "darkgray", xpd = T,
           xlim=c(min(prevalence_table$date)-60,
                  max(prevalence_table$date)+60),
           ylim = c(0, ( max(prevalence_table$proportion)) + (max(prevalence_table$proportion)) * 0.2),
           xlab = "Collection Date", ylab = "Prevalence Estimate", main = lineage)
      abline(v = crest_date, col = "green", lwd=2)
      abline(v = rare_date, col = "red", lwd=2)
      legend("topleft", c(
        paste("Crest Date: ", crest_date, " Crest Prevalence: ", round(max(prevalence_table$proportion), digits = 4)),
        paste("Rarity Date: ", rare_date, " Rarity Prevalence: ", round(rare_prev, digits = 8))
      ), col = c("green", "red"), lwd = c(2,2), bty = 'n')
      
      # if the lineage crested in the past 30 days, assume it cannot be
      # anachronistic (yet).
      if ( (Sys.Date() - crest_date) < 30 ){
        
        rare_date <- as.Date(NA)
        
      }
      
      # if the lineage has not yet crested, by definition it cannot be
      # anachronistic
      if (nrow(filter(prevalence_table, date > crest_date))==0){
        
        rare_date <- as.Date(NA)
        
      }
      
    } else {
      
      # assume the lineage has not crested if the length of the prevalence table
      # is zero
      rare_date <- as.Date(NA)
      
    }
    
    if (!is.na(rare_date)){
      stopifnot(rare_date > crest_date)
    }
    
    return(rare_date)
    
  }, Date(1))
  dev.off()
  
  # coerce the lookup to be properly date-formatted
  rarity_lookup <- as.Date(rarity_lookup)
  # make sure each item can be accessed by lineage name
  names(rarity_lookup) <- lineages
  
  return(rarity_lookup)
  
})

# define metadata function
subset_metadata <- cmpfun(function(metadata_path, report_path){
  
  # gathering long table of U.S. HHS regions
  # data(hhs_regions)
  # hhs_regions <- data.frame(hhs_regions)
  
  # Read in metadata
  metadata <- read_tsv(metadata_path, show_col_types = FALSE, trim_ws = TRUE)
  
  # collate pangolin metadata: new lineage calls, accession IDs, and collection dates
  pango_report <- date_pango_calls(report_path, metadata)
  
  # sort metadata so that accessions is in the same order as the pango report
  metadata <- arrange(metadata, Accession)
  stopifnot(!(FALSE %in% (metadata$Accession == pango_report$`taxon`)))
  
  # create rare-date lookup with the outbreak.info API
  rarity_lookup <- create_rarity_lookup(pango_report)
  
  # create vectorized function for identifying the HHS regions and determining
  # whether the probability that a given lineage has persisted is marginal
  metadata$`Anachronicity (days)` <- vapply(1:nrow(pango_report), 
                                            function(i, lineages, dates, rarity_lookup){
                                              
                                              # define anachronicity in days
                                              anachronicity <- as.numeric(dates[i] - rarity_lookup[lineages[i]])
                                              
                                              # Here an ifelse statement checks if the collection date for a the lineage
                                              # sample in question comes after the above-computed "final date of rarity"
                                              ifelse(!is.na(anachronicity) && anachronicity > 0, 
                                                     return(anachronicity), 
                                                     return(0))
                                              
                                            }, lineages = pango_report$lineage, 
                                            dates = pango_report$`Isolate Collection date`,
                                            rarity_lookup = rarity_lookup,
                                            numeric(1))
  
  metadata <- metadata %>%
    filter(`Anachronicity (days)` > 0)
  
  return(metadata)
  
})

# define FASTA function
subset_fasta <- cmpfun(function(metadata, fasta_path){
  
  # bring in FASTA
  all_seqs <- readDNAStringSet(fasta_path)
  
  # parse out names so that they only contain accession IDs
  names(all_seqs) <- str_split_i(names(all_seqs),
                                 " ",
                                 i = 1)
  
  # parse out anachronistic seqs
  anachron_seqs <- all_seqs[names(all_seqs) %in% metadata$Accession]
  
  stopifnot(length(anachron_seqs) == nrow(metadata))
  
  return(anachron_seqs)
  
})

# define main function that will bring it all together
main <- cmpfun(function(report_path, metadata_path, fasta_path){
  
  # filter the metadata
  new_metadata <- subset_metadata(metadata_path, report_path)
  
  # filter the sequences
  new_fasta <- subset_fasta(new_metadata, fasta_path)
  
  # write both
  writeXStringSet(new_fasta, "anachronistic_seq_candidates.fasta")
  new_metadata <- new_metadata[order(new_metadata$`Anachronicity (days)`, decreasing = T),]
  write.table(new_metadata, "anachronistic_seq_candidates.tsv", 
              row.names = F, quote = F, sep = "\t")
  
})

# run the hierarchy of functions
main(report_path, metadata_path, fasta_path)
