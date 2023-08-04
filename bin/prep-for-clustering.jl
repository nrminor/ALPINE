#!/usr/bin/env julia

# saving command line arguments supplied by nextflow as immutable
# constants for more performant static typing and memory usage
const metadata_path = islink(ARGS[1]) ? readlink(ARGS[1]) : ARGS[1]
const geo = ARGS[2]
const fasta_path = islink(ARGS[3]) ? readlink(ARGS[3]) : ARGS[3]

# run wrapper function for all sequence clustering prep
prep_for_clustering(metadata_path,geo,fasta_path)
