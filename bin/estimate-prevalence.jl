#!/usr/bin/env -S julia -t auto

using CSV, DataFrames

# define input paths from the command line as typed constants
const STATS_FILE1::String = islink(ARGS[1]) ? readlink(ARGS[1]) : ARGS[1]
const STATS_FILE2::String = islink(ARGS[2]) ? readlink(ARGS[2]) : ARGS[2]

"""
Function `estimate_prevalence` uses the outputs from `seqkit stats` to
estimate the proportion of sequences flagged at the end of ALPINE in the
sequences used as the workflow's input.
"""
function estimate_prevalence(early_stats::DataFrame, late_stats::DataFrame)

    # use an assertion to make sure the necessary column names are present
    @assert "num_seqs" in names(early_stats)
    @assert "num_seqs" in names(late_stats)

    # make sure that candidates were found at all before pulling out
    # the final sample size
    candidate_count = nrow(late_stats) == 0 ? 0 : unique(late_stats.num_seqs)

    # pull out the input sample size
    sample_size = unique(early_stats.num_seqs)[1]

    # compute the prevalence estimate
    prevalence = join(round((candidate_count / sample_size); digits=5) * 100, ";")

    return (prevalence, sample_size)

end

function main(file1::String, file2::String)

    # handle file I/O here
    early_stats = CSV.read(file1, DataFrame)
    late_stats = CSV.read(file2, DataFrame)

    # estimate prevalence
    prevalence, sample_size = estimate_prevalence(early_stats, late_stats)

    println("$prevalence% of $sample_size sequences were flagged as highly evolved and anachronistic.")
    
end

main(STATS_FILE1, STATS_FILE2)
