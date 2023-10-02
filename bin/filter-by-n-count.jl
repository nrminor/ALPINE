#!/usr/bin/env -S julia -t auto

using ALPINE, Pipe

# saving command line arguments supplied by nextflow
const FASTA_PATH::String = islink(ARGS[1]) ? readlink(ARGS[1]) : ARGS[1]
const MAX_AMBIGUITY::Float64 = parse(Float64, ARGS[2])
const REF_PATH::String = islink(ARGS[3]) ? readlink(ARGS[3]) : ARGS[3]

function main(fasta::String, max_ambiguity::Float64, ref_path::String)

    @assert sum(occursin(">", line) for line in eachline(fasta)) > 0 "FASTA is empty or corrupted and thus cannot be filtered."

    @pipe replace_gaps(fasta) |>
    filter_by_n(_, max_ambiguity, ref_path)

end
precompile(main, (String, Float64, String,));

main(FASTA_PATH, MAX_AMBIGUITY, REF_PATH)
