#!/usr/bin/env julia

# loading necessary packages
using CSV, DataFrames, FastaIO, FileIO
import Base.Threads.@spawn

# locate metadata
metadata_path = ARGS[1]

# Check if the input metadata file is a symlink, and if it is, follow the symlink
if islink(metadata_path)
    metadata_path = readlink(metadata_path)
end

# locate FASTA
fasta_path = ARGS[2]

# Check if the input FASTA file is a symlink, and if it is, follow the symlink
if islink(fasta_path)
    fasta_path = readlink(fasta_path)
end

# define output file name
output_handle = ARGS[3]
touch(output_handle)

# read in metadata
metadata_df = CSV.read(metadata_path, DataFrame, delim="\t")

# create lookup of dates and accessions
accession_to_date = Dict(zip(metadata_df[!,"Accession"], metadata_df[!,"Isolate Collection date"]))

# add date for each FASTA record
FastaWriter(output_handle, "a") do fa
    FastaReader(fasta_path) do reader 
        for (name, seq) in reader
            # Parse out the accession number from the defline
            accession = split(name, ' ')[1]
            # Look up the collection date for the accession number
            collection_date = get(accession_to_date, accession, "")
            # Add the collection date to the defline after a pipe symbol
            new_name = string(name, "|", collection_date)
            # Write the modified record to a new FASTA file
            writeentry(fa, new_name, seq)
        end
    end
end