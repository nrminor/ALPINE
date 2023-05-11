#!/usr/bin/env julia

# loading necessary packages
using CSV, DataFrames, BioSequences

# locate metadata
metadata_path = ARGS[1]

# locate FASTA
fasta_path = ARGS[2]

# define output file name
output_handle = ARGS[3]

# read in metadata
metadata = CSV.read(metadata_path, DataFrame, delim="\t", follow_symlinks=true)

# create lookup of dates and accessions
accession_to_date = Dict(zip(metadata[!,"Accession"], metadata[!,"Isolate Collection date"]))

# add date for each FASTA record
open(FASTA.Reader, fasta_path) do reader 
    for record in reader
        # Parse out the accession number from the defline
        accession = record.identifier
        # Look up the collection date for the accession number
        collection_date = get(accession_to_date, accession, "")
        # Add the collection date to the defline after a pipe symbol
        record.description = string(record.description, "|", collection_date)
        # Write the modified record to a new FASTA file
        open(FASTA.Writer, output_handle, append=true) do w
            write(w, record)
        end

    end
end