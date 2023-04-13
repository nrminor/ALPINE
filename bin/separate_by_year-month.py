#!/usr/bin/env python3

import sys
from Bio import SeqIO
from datetime import datetime
from multiprocessing import Pool, cpu_count

# Define a function to process a single sequence record
def process_seq_record(record):
    # Extract the date from the FASTA header
    header = record.description
    date_str = header.split("|")[-1].strip()
    
    try:
        date = datetime.strptime(date_str, "%Y-%m-%d")
    except ValueError:
        # Skip this sequence record if the date isn't parseable
        return None

    # Create a filename for the output file based on the year and month of the sequence
    year_month = f"{date.year}-{date.month:02}"
    output_file = f"{year_month}.fasta"

    # Write the record to the appropriate output file
    return (output_file, record)

# Get the input filename from the command line arguments
input_file = sys.argv[1]

# Define the number of processes to use (default to number of CPU cores)
num_processes = int(sys.argv[2]) if len(sys.argv) > 2 else cpu_count()

# Read in the input FASTA file
with open(input_file, "r") as in_handle:
    seq_records = list(SeqIO.parse(in_handle, "fasta"))

# Process the sequence records in parallel
with Pool(processes=num_processes) as pool:
    results = pool.map(process_seq_record, seq_records)

# Create a dictionary to hold the output file handles
output_files = {}

# Loop over the results and write the records to the appropriate output files
for result in results:
    if result is None:
        continue
    
    output_file, record = result
    
    # If the output file hasn't been opened yet, open it in write mode
    if output_file not in output_files:
        output_files[output_file] = open(output_file, "w")

    # Write the record to the appropriate output file
    SeqIO.write(record, output_files[output_file], "fasta")

# Close all the output files
for output_file in output_files.values():
    output_file.close()
