#!/usr/bin/env -S julia -t auto

# loading necessary packages
using ALPINE, FastaIO, FileIO, Dates, DataFrames, Arrow

# bring in command line arguments as static-typed constants, while
# checking for and following any symlinks
const FASTA_PATH::String = islink(ARGS[1]) ? readlink(ARGS[1]) : ARGS[1]
const TABLE_PATH::String = islink(ARGS[2]) ? readlink(ARGS[2]) : ARGS[2]

# create lookup of dates and accessions and save as a typed constant
const ACCESSION_TO_DATE::Dict = date_accessions(TABLE_PATH)

# run the function
separate_by_month(FASTA_PATH, ACCESSION_TO_DATE)
