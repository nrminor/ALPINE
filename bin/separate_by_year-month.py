#!/usr/bin/env python3

import sys
from Bio import SeqIO
from datetime import datetime

# Get the input filename from the command line arguments
input_file = sys.argv[1]

# Define the date format expected in the FASTA headers
date_format = "%Y-%m-%d"

# Create a dictionary to hold the output file handles
output_files = {}

# Open the input file and iterate over the records
for record in SeqIO.parse(input_file, "fasta"):
	# Extract the date from the FASTA header
	header = record.description
	date_str = header.split("|")[-1].strip()
	
	try:
		date = datetime.strptime(date_str, "%Y-%m-%d")
	except ValueError:
		# Skip this sequence record if the date isn't parseable
		continue

	# Create a filename for the output file based on the year and month of the sequence
	year_month = f"{date.year}-{date.month:02}"
	output_file = f"{year_month}.fasta"

	# If the output file hasn't been opened yet, open it in write mode
	if output_file not in output_files:
		output_files[output_file] = open(output_file, "w")

	# Write the record to the appropriate output file
	SeqIO.write(record, output_files[output_file], "fasta")

# Close all the output files
for output_file in output_files.values():
	output_file.close()
