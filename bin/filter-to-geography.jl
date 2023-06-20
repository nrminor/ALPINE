#!/usr/bin/env julia

# loading packages
using LongInfectionFinder, DelimitedFiles, DataFrames, CSV

# saving command line arguments supplied by nextflow
const metadata_path = islink(ARGS[1]) ? readlink(ARGS[1]) : ARGS[1]
const geo = ARGS[2]

# run the function
filter_by_geo(metadata_path,geo)