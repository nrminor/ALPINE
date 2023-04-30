#!/usr/bin/env python3

import sys
import pandas as pd
from Bio import SeqIO

# parse command line arguments
if len(sys.argv) != 3:
    sys.exit("Usage: python3 append_dates.py <input_metadata.tsv> <input.fasta> <desired_output_filename.fasta>")

# read in arguments
metadata = sys.argv[1]
fasta = sys.argv[2]
output_handle = sys.argv[3]

# Read in the metadata TSV file as a pandas dataframe
metadata_df = pd.read_csv(metadata, sep="\t")

# Create a dictionary that maps accession numbers to collection dates
accession_to_date = dict(zip(metadata_df["Accession"], metadata_df["Isolate Collection date"]))

# Open the FASTA file and iterate over each record
with open(fasta, "r") as fasta_file:
    for record in SeqIO.parse(fasta_file, "fasta"):
        # Parse out the accession number from the defline
        accession = record.id.split()[0]
        # Look up the collection date for the accession number
        collection_date = accession_to_date.get(accession, "")
        # Add the collection date to the defline after a pipe symbol
        record.id = f"{record.id}|{collection_date}"
        # Write the modified record to a new FASTA file
        with open(output_handle, "a") as output_file:
            SeqIO.write(record, output_file, "fasta")