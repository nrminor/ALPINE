#!/usr/bin/env julia

using ALPINE

# saving command line arguments supplied by nextflow
const FASTA_PATH::String =  islink(ARGS[1]) ? readlink(ARGS[1]) : ARGS[1]
const OUT_NAME::String = "no-gaps.fasta.gz"

# run the function
replace_gaps(FASTA_PATH, OUT_NAME)
