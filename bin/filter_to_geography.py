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

# Read the FASTA file and write matching sequences to output file in real time
with open('filtered_to_geography.fasta', 'w') as output_handle:
    for record in SeqIO.parse(fasta, 'fasta'):
        if record.id.split()[0] in selected_rows['Accession'].tolist():
            SeqIO.write(record, output_handle, 'fasta')

# Close the output FASTA file
output_handle.close()
