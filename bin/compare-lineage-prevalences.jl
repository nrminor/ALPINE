#!/usr/bin/env julia

using ALPINE, RCall, CSV, DataFrames

# define typed global constants based on command line arguments
const REPORT_PATH::String = islink(ARGS[1]) ? readlink(ARGS[1]) : ARGS[1]
const METADATA_PATH::String = islink(ARGS[2]) ? readlink(ARGS[2]) : ARGS[2]
const FASTA_PATH::String = islink(ARGS[3]) ? readlink(ARGS[3]) : ARGS[3]

function main(report_path::String, metadata_path::String, fasta_path::String)

    # printing GISAID Disclaimer; authenticate with R"authenticateUser()"
    println("This data was obtained from GISAID via the outbreak.info API.")

    # filter the metadata
    new_metadata = assign_anachron(metadata_path, report_path)
    sort!(new_metadata, Symbol("Anachronicity (days)"), rev = true)

    # write metadata (the FASTA was already written incrementally in subset_fasta)
    CSV.write("anachronistic_seq_candidates.tsv", new_metadata, delim='\t')

    # filter the sequences
    find_anachron_seqs(new_metadata, fasta_path)

end
precompile(main, (string,string,string))

main(REPORT_PATH, METADATA_PATH, FASTA_PATH)
