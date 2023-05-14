#!/usr/bin/env julia

# loading packages
using DelimitedFiles, DataFrames, CSV

# saving command line arguments supplied by nextflow
metadata_path = ARGS[1]
geography = string(ARGS[2])

# Read in the TSV file with metadata
metadata_df = CSV.read(metadata_path, DataFrame, delim="\t")

# Double check the column name for geographic locations
if "Geographic location" in names(metadata_df)
    idx = findfirst(isequal("Geographic location"), names(metadata_df))
    rename!(metadata_df, names(metadata_df)[idx] => "Geographic Location")
end

# filter metadata based on desired geography
filtered = metadata_df[[contains(string(value), geography) for value in metadata_df[!,"Geographic Location"]], :]

# Writing filtered metadata
CSV.write("filtered_to_geography.tsv", filtered, delim="\t")

# separating out accessions
accessions = filtered[!,"Accession"]

# Writing accessions to a text file for use by seqtk
writedlm("accessions.txt", accessions, "\n")