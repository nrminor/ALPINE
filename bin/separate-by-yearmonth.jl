#!/usr/bin/env julia

# loading necessary packages
using LongInfectionFinder, FastaIO, FileIO, Dates, CSV, DataFrames, Arrow

# Check if the input FASTA file is a symlink, and if it is, follow the symlink
const fasta_path = islink(ARGS[1]) ? readlink(ARGS[1]) : ARGS[1]

# locate metadata
const input_table_path = islink(ARGS[2]) ? readlink(ARGS[2]) : ARGS[2]

# read in metadata
metadata_df = DateFrame(Arrow.Table(input_table_path, DataFrame))

# create lookup of dates and accessions
accession_to_date = Dict(zip(metadata_df[!,"Accession"], metadata_df[!,"Isolate Collection date"]))

# run the function
separate_by_month(fasta_path, accession_to_date)
