#!/usr/bin/env -S julia -t auto

using ALPINE, CSV, DataFrames

# define input paths from the command line as typed constants
const STATS_FILE1::String = islink(ARGS[1]) ? readlink(ARGS[1]) : ARGS[1]
const STATS_FILE2::String = islink(ARGS[2]) ? readlink(ARGS[2]) : ARGS[2]

function main(file1::String, file2::String)

    # handle file I/O here
    early_stats = CSV.read(file1, DataFrame)
    late_stats = CSV.read(file2, DataFrame)

    # estimate prevalence
    prevalence, sample_size = estimate_prevalence(early_stats, late_stats)

    println("$prevalence% of $sample_size sequences were flagged.")
    
end
precompile(main, (String, String))

main(STATS_FILE1, STATS_FILE2)
