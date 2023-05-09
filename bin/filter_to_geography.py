#!/usr/bin/env python3

import sys
import polars as pl

# parse information supplied by nextflow as arguments
metadata = sys.argv[1]
geography = str(sys.argv[2])

# Load the TSV file into a DataFrame
df = pl.read_csv(metadata, separator='\t')

# Filter for rows containing geography in the "Geographic location" column
selected_rows = df.filter(pl.col('Geographic location').str.contains(geography).alias('literal'))

# Write the matching "Accession" numbers to a separate file
accessions = selected_rows['Accession'].to_list()
with open('accessions.txt', 'w') as f:
    f.write('\n'.join(accessions))

# Close the output accession list file
f.close()

# Save the filtered DataFrame to a new TSV file
selected_rows.write_csv('filtered_to_geography.tsv', separator='\t')
