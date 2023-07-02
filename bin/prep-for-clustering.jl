#!/usr/bin/env julia

# saving command line arguments supplied by nextflow
const metadata_path = islink(ARGS[1]) ? readlink(ARGS[1]) : ARGS[1]
const geo = ARGS[2]
const fasta_path = islink(ARGS[3]) ? readlink(ARGS[3]) : ARGS[3]

prep_for_clustering(metadata_path,geo,fasta_path)