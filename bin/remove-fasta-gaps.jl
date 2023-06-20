#!/usr/bin/env julia

using LongInfectionFinder

# saving command line arguments supplied by nextflow
fasta_path = ARGS[1]
output_filename = "no-gaps.fasta"

# Check if the input FASTA file is a symlink, and if it is, follow the symlink
if islink(fasta_path)
    fasta_path = readlink(fasta_path)
end

# run the function
replace_gaps(fasta_path, output_filename)