#!/usr/bin/env julia

# loading necessary packages
using FastaIO, FileIO

# saving command line arguments supplied by nextflow
fasta_path = ARGS[1]

# Check if the input FASTA file is a symlink, and if it is, follow the symlink
if islink(fasta_path)
    fasta_path = readlink(fasta_path)
end

# create the output file
touch("filtered-by-n.fasta")

# creating a loop that goes through sequence records and writes out
# any sequences that have less than the minimum N count
FastaWriter("filtered-by-n.fasta" , "a") do fa
    FastaReader(fasta_path) do fr
        for (name, seq) in fr
            max_n_count = floor(length(seq) * 0.1)
            n_count = count("N", convert(String, seq))
            if n_count < max_n_count
                writeentry(fa, name, seq)
            end
        end
    end
end
