#!/usr/bin/env julia

using LongInfectionFinder

# saving command line arguments supplied by nextflow
fasta_path =  islink(ARGS[1]) ? readlink(ARGS[1]) : ARGS[1]
output_filename = "no-gaps.fasta.gz"

# run the function
replace_gaps(fasta_path, output_filename)