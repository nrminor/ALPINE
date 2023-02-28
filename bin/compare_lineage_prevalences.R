#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# loading necessary libraries
library(parallel)
library(tidyverse)
library(magrittr)
library(knitr)
library(lubridate)
library(cdcfluview)
library(outbreakinfo)

# Bringing in pangolin report csv 
pango_report <- read.csv(args[1])

# bringing in available cores supplied by Nextflow
available_cores <- as.numeric(args[2])

# printing GISAID Disclaimer
print("This data was obtained from GISAID via the outbreak.info API")

# gathering long table of U.S. HHS regions
data(hhs_regions)
hhs_regions <- data.frame(hhs_regions)

# Preparing a cluster of cores provided by nextflow and exporting the necessary
# data to each core
cluster <- makeCluster(available_cores-1)
taxa <- pango_report$taxon
# accessions <- pango_report$accession
us_state <- str_remove_all(pango_report$location, "US: ")
lineages <- pango_report$lineage
dates <- pango_report$collection_date
clusterExport(cluster, 
              varlist=c("hhs_regions", "taxa", "us_state", "lineages", "dates"))

# create parallelized function for identifying the HHS regions and determining
# whether the probability that a given lineage has persisted is marginal
pango_report$anachronistic <- do.call(c, 
                                      clusterApplyLB(cluster, 1:length(dates), function(i){
                                        
                                        # identifying all states in each HHS region. All
                                        # prevalences will be based on U.S. HHS region
                                        # rather than individual states. This is because
                                        # data reporting and availability is variable
                                        # across states.
                                        hhs <- hhs_regions[hhs_regions$state_or_territory==us_state[i], "region"]
                                        region_states <- hhs_regions[hhs_regions$region==hhs,
                                                                     "state_or_territory"]
                                        
                                        # creating an empty table where we will fill in
                                        # prevalences
                                        prevalence_table <- getPrevalence(pangolin_lineage = lineages[i], 
                                                                          location = "Wisconsin")
                                        prevalence_table <- prevalence_table[NULL,]
                                        
                                        # looping through each state in the HHS region
                                        # in question and gathering prevalence estimates
                                        # from outbreak.info/GISAID
                                        for (j in region_states){
                                          
                                          new_rows <- getPrevalence(pangolin_lineage = lineages[i], 
                                                                    location = j)
                                          prevalence_table <- rbind(prevalence_table, new_rows)
                                          
                                        }
                                        
                                        # Then, we get rid of zeroes, which will skew
                                        # our methods below
                                        prevalence_table <- prevalence_table[prevalence_table$proportion>0,]
                                        rownames(prevalence_table) <- NULL
                                        
                                        # Next, we use a rather crude method to identify
                                        # when the wave for a lineage in question has
                                        # "crested". This is the date of peak prevalence
                                        # for the lineage in the HHS region.
                                        crest <- max(prevalence_table[
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
                                        rare <- max(prevalence_table[
                                          prevalence_table$proportion<=quantile(prevalence_table$proportion, 
                                                                                0.1), "date"])
                                        
                                        # Finally, an ifelse statement checks if the
                                        # collection date for a the lineage sample in
                                        # question comes after the above-computed "final date
                                        # of rarity"
                                        ifelse(dates[i] > rare, return(TRUE), return(FALSE))
                                        
                                      }))

