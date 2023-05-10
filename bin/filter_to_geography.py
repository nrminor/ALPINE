#!/usr/bin/env python3

import sys
import os
import polars as pl

# parse information supplied by nextflow as arguments
if os.path.islink(sys.argv[1]):
    metadata = os.readlink(sys.argv[1])
else:
    metadata = sys.argv[1]
geography = str(sys.argv[2])

# Load the TSV file into a DataFrame
df = pl.read_csv(metadata, separator='\t', has_header=True)

# check if any column names are "Geographic location" or "Geographic Location"
if "Geographic Location" in df.columns:
    # if the column name is already "Geographic Location", do nothing
    pass
elif "Geographic location" in df.columns:
    # if the column name is "Geographic location", rename it to "Geographic Location"
    df = df.rename({"Geographic location": "Geographic Location"})
else:
    # if neither column name is found, do nothing
    pass

# Filter for rows containing geography in the "Geographic location" column
selected_rows = df.filter(
    pl.col('Geographic Location').str.contains(geography).alias('literal')
)

# Write the matching "Accession" numbers to a separate file
accessions = selected_rows['Accession'].to_list()
with open('accessions.txt', 'w') as f:
    f.write('\n'.join(accessions))

# Close the output accession list file
f.close()

# Save the filtered DataFrame to a new TSV file
selected_rows.write_csv('filtered_to_geography.tsv', separator='\t')
