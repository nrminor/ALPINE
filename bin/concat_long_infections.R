#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

csv_list <- list.files(getwd(), "*.csv")

days_of_infection <- as.numeric(args[1])

for (i in 1:length(csv_list)){
	
	csv <- read.csv(csv_list[i], header = T)
	
	if (i == 1 & !grepl("No putative long infections", csv[1,1]) ){
		
		long_inf_table <- csv
		
	} else if (i == 1 & grepl("No putative long infections", csv[1,1]) ){
		
		long_inf_table <- csv
		long_inf_table <- long_inf_table[FALSE, ]
		
	}
	
	if (nrow(csv) > 0 && !grepl("No putative long infections", csv[1,1])) {
		
		long_inf_table <- rbind(long_inf_table, csv)
		remove(csv)
		
	} else {
		
		remove(csv)
		
	}
}

if (nrow(long_inf_table)==0){
	
	null_row <- rep(NA, times = ncol(long_inf_table))
	null_row[1] <- "No putative long infections were identified in this dataset."
	null_row <- t(as.data.frame(null_row)) ; rownames(null_row) <- NULL
	long_inf_table <- rbind(long_inf_table, null_row)
	
}

if (nrow(long_inf_table)>0){
	
	long_inf_table$days_to_define_prolonged_infections <- days_of_infection
	
}

write.csv(long_inf_table, paste("putative_long_infections_", Sys.Date(), ".csv", sep = ""), row.names = F, quote = F, na = "")
