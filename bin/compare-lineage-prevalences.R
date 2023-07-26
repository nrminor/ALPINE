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
  pango_report <- read_csv(report_path, col_select = c(`Sequence name`, Lineage),
                           show_col_types = FALSE, trim_ws = TRUE) %>%
    filter(!grepl("taxon", `Sequence name`))
  
  # joining the desired columns
  pango_with_dates <- pango_report %>% 
    left_join(metadata %>% select(Accession, `Isolate Collection date`), 
              by = c("Sequence name" = "Accession"), keep = TRUE) %>%
    select(-Accession) %>%
    arrange(`Sequence name`)
  
  return(pango_with_dates)
  
})

# define metadata function
subset_metadata <- cmpfun(function(metadata_path, report_path){
  
  # gathering long table of U.S. HHS regions
  data(hhs_regions)
  hhs_regions <- data.frame(hhs_regions)
  
  # Read in metadata
  metadata <- read_tsv(metadata_path, show_col_types = FALSE, trim_ws = TRUE)
  
  # collate pangolin metadata: new lineage calls, accession IDs, and collection dates
  pango_report <- date_pango_calls(report_path, metadata)
  
  # sort metadata so that accessions is in the same order as the pango report
  metadata <- arrange(metadata, Accession)
  stopifnot(!(FALSE %in% (metadata$Accession == pango_report$`Sequence name`)))
  
  # create vectorized function for identifying the HHS regions and determining
  # whether the probability that a given lineage has persisted is marginal
  metadata$`Anachronicity (days)` <- vapply(1:nrow(pango_report), 
                                            function(i, lineages, dates){
    
    # identifying all states in each HHS region. All
    # prevalences will be based on U.S. HHS region
    # rather than individual states. This is because
    # data reporting and availability is variable
    # across states.
    # hhs <- hhs_regions[hhs_regions$state_or_territory==us_state[i], "region"]
    # region_states <- hhs_regions[hhs_regions$region==hhs,
    # "state_or_territory"]
    
    # creating an empty table where we will fill in
    # prevalences
    prevalence_table <- getPrevalence(pangolin_lineage = lineages[i], 
                                      location = "United States")
    # prevalence_table <- prevalence_table[NULL,]
    
    # looping through each state in the HHS region
    # in question and gathering prevalence estimates
    # from outbreak.info/GISAID
    # for (j in region_states){
    #   
    #   new_rows <- getPrevalence(pangolin_lineage = lineages[i], 
    #                             location = j)
    #   prevalence_table <- rbind(prevalence_table, new_rows)
    #   
    # }
    
    # Then, we get rid of zeroes, which will
    # skew our methods below. We also get
    # rid of 1.0's, which are likely
    # erroneous.
    prevalence_table <- prevalence_table[prevalence_table$proportion > 0 &
                                           prevalence_table$proportion < 1,]
    
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
    q <- quantile(prevalence_table$proportion, 0.01)
    rare_prev <- prevalence_table$proportion[
      which.min(
        abs(prevalence_table$proportion - q))]
    rare_date <- max(prevalence_table[
      prevalence_table$proportion==rare_prev, "date"]) # + 30
    
    # define anachronicity in days
    anachronicity <- as.numeric(dates[i] - rare_date)
    
    # Finally, an ifelse statement checks if the
    # collection date for a the lineage sample in
    # question comes after the above-computed "final date
    # of rarity"
    ifelse(anachronicity > 0, 
           return(anachronicity), 
           return(0))
    
  }, lineages = pango_report$Lineage, 
  dates = pango_report$`Isolate Collection date`,
  numeric(1))
  
  metadata <- metadata %>%
    filter(`Anachronicity (days)` > 0)
  
  return(metadata)
  
})

# define FASTA function
subset_fasta <- cmpfun(function(metadata, fasta_path){
  
  # bring in FASTA
  all_seqs <- readDNAStringSet(fasta_path)
  
  # parse out anachronistic seqs
  anachron_seqs <- all_seqs[names(all_seqs) %in% metadata$Accession]
  
  return(anachron_seqs)
  
})

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
