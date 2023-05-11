#!/usr/bin/env julia

# loading necessary packages
using BioSequences

# saving command line arguments supplied by nextflow
fasta_path = ARGS[1]

# creating a loop that goes through sequence records and writes out
# any sequences that have less than the minimum N count
open(FASTA.Reader, fasta_path) do reader 
    for record in reader
        max_n_count = floor(length(record.sequence) * 0.1)
        n_count = count("N", convert(String, FASTA.sequence(record)))
        if n_count < max_n_count
            open(FASTA.Writer, "filtered-by-n.fasta", append=true) do w
                write(w, record)
            end
        end
    end
end
