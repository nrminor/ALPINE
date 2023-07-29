#!/usr/bin/env julia

# loading necessary packages
using LongInfectionFinder, FastaIO, BioSequences, Distances, Statistics, DataFrames, CSV

# parse supplied command line arguments to locate files and set parameters
const fasta_file = islink(ARGS[1]) ? readlink(ARGS[1]) : ARGS[1]
const table_path = islink(ARGS[2]) ? readlink(ARGS[2]) : ARGS[2] # this will be used in the future
const yearmonth = ARGS[3]
const stringency = ARGS[4]

# replace lowercase n symbols with uppercase Ns
const tmp = "tmp.fasta"
touch(tmp)
set_to_uppercase(fasta_file,tmp)

# read in the cluster table
cluster_table = CSV.read(table_path, DataFrame, delim="\t", header=false)

# define majority centroid
centroids = filter(1 => ==("C"), cluster_table)
const majority_centroid = convert(String, centroids[argmax(centroids[:,3]), :Column9])

# count the number of clusters
const count = nrow(centroids)

# run the function
distance_matrix(tmp, cluster_table, count, majority_centroid, yearmonth, stringency)
