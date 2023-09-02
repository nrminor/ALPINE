#!/usr/bin/env -S julia -t auto

# loading necessary packages
using FastaIO, FileIO, ALPINE
import Base.Threads

# Check if the input FASTA file is a symlink, and if it is, follow the symlink,
# and specify the output filename
const FASTA_PATH::String = islink(ARGS[1]) ? readlink(ARGS[1]) : ARGS[1]
const MAX_AMBIGUITY::Float64 = parse(Float64, ARGS[2])
const OUT_NAME::String = "filtered-by-n.fasta.gz"

# run the function
filter_by_n(FASTA_PATH, MAX_AMBIGUITY, OUT_NAME)
