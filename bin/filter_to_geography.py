#!/usr/bin/env python3

import sys
import pandas as pd
from Bio import SeqIO

# parse information supplied by nextflow as arguments
metadata = sys.argv[1]
fasta = sys.argv[2]
geography = str(sys.argv[3])

# Read the spreadsheet
df = pd.read_csv(metadata, sep='\t')

# Filter the spreadsheet for rows containing 'USA:'
selected_rows = df[df['Geographic location'].str.contains(geography, na=False)]

# Write the filtered spreadsheet to a new file
selected_rows.to_csv('filtered_to_geography.tsv', sep='\t', index=False)

# Read the FASTA file
fasta_records = list(SeqIO.parse(fasta, 'fasta'))

# Filter the FASTA records for those corresponding to USA accessions
selected_accessions = selected_rows['Accession'].tolist()
selected_fastas = [record for record in fasta_records if record.id.split()[0] in selected_accessions]

# Write the filtered FASTA records to a new file
with open('filtered_to_geography.fasta', 'w') as output_handle:
    SeqIO.write(selected_fastas, output_handle, 'fasta')
