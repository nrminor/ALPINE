#!/usr/bin/env julia

# loading necessary packages
using LongInfectionFinder, FastaIO, BioSequences, Distances, Statistics, DataFrames, CSV

# parse supplied command line arguments to locate files and set parameters
const fasta_file = islink(ARGS[1]) ? readlink(ARGS[1]) : ARGS[1]
const table_path = islink(ARGS[2]) ? readlink(ARGS[2]) : ARGS[2] # this will be used in the future
const yearmonth = ARGS[3]
const count = parse(Int, ARGS[4])
const majority_centroid = ARGS[5]

# replace lowercase n symbols with uppercase Ns
const tmp = "tmp.fasta"
touch(tmp)
set_to_uppercase(fasta_file,tmp)

# read in the cluster table
cluster_table = CSV.read(table_path, DataFrame, delim="\t", header=false)

# run the function
distance_matrix(tmp, cluster_table, yearmonth)
