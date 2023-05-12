#!/usr/bin/env julia

# loading necessary packages
using FastaIO

# saving command line arguments supplied by nextflow
fasta_path = ARGS[1]

# creating a loop that goes through sequence records and writes out
# any sequences that have less than the minimum N count
FastaReader(fasta_path) do fr
    for (name, seq) in fr
        max_n_count = floor(length(seq) * 0.1)
        n_count = count("N", convert(String, seq))
        if n_count < max_n_count
            FastaWriter("filtered-by-n.fasta" , "a") do fa
                writeentry(fa, name, seq)
            end
        end
    end
end
