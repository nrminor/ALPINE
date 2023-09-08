#!/usr/bin/env -S julia -t auto

# loading necessary packages
using ALPINE, FastaIO, BioSequences, Distances, Statistics, DataFrames, CSV

# parse supplied command line arguments to locate files and set parameters
const FASTA_PATH::String = islink(ARGS[1]) ? readlink(ARGS[1]) : ARGS[1]
const TABLE_PATH::String = islink(ARGS[2]) ? readlink(ARGS[2]) : ARGS[2] # this will be used in the future
const YEARMONTH::String = ARGS[3]
const STRINGENCY::String = ARGS[4]

function main(fasta_file::String, table_path::String, yearmonth::String, stringency::String)

    # replace lowercase n symbols with uppercase Ns
    tmp = "tmp.fasta"
    touch(tmp)
    set_to_uppercase(fasta_file,tmp)

    # read in the cluster table
    cluster_table = CSV.read(table_path, DataFrame, delim="\t", header=false)

    # run the function
    distance_matrix(tmp, cluster_table, yearmonth, stringency)

end

# precompile main
precompile(main, (String,String,String,String,));

# run main
main(FASTA_PATH, TABLE_PATH, YEARMONTH, STRINGENCY)
