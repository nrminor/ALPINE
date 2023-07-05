#!/usr/bin/env julia

# loading necessary packages
using FastaIO, FileIO, LongInfectionFinder
import Base.Threads

# Check if the input FASTA file is a symlink, and if it is, follow the symlink
const fasta_path = islink(ARGS[1]) ? readlink(ARGS[1]) : ARGS[1]

# create the output file
const out = "filtered-by-n.fasta.gz"

# run the function
filter_by_n(fasta_path, out)