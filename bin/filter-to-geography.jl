#!/usr/bin/env julia

# loading packages
using DataFrames, CSV

# saving command line arguments supplied by nextflow
metadata_path = ARGS[1]
geography = ARGS[2]

# Read in the TSV file with metadata
metadata = CSV.read(metadata_path, DataFrame, delim="\t", follow_symlinks=true)

# Double check the column name for geographic locations
if "Geographic location" in names(metadata)
    idx = findfirst(isequal("Geographic location"), names(metadata))
    rename!(metadata, names(metadata)[idx] => "Geographic Location")
end

# filter metadata based on desired geography
filtered = metadata[[contains(geography, string(value)) for value in metadata[:"Geographic Location"]], :]

# separating out accessions
accessions = filtered[:"Accession"]

# Writing accessions to a text file for use by seqtk
CSV.write("accessions.txt", accessions, delim="\t")