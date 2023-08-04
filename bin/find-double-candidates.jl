#!/usr/bin/env julia

# load pipeline module
using LongInfectionFinder, CSV, Dataframes

# load I/O
const metadata_files = filter(x->occursin(".tsv",x), readdir(abspath("."), join=true))
const fasta_files = filter(x->occursin(".fasta",x), readdir(abspath("."), join=true))

# use multiple dispatch to collect double candidates based on the 
# number of available files
function main(metadata_files::Vector{String}, fasta_files::Vector{String})

    num_meta = length(metadata_files)
    num_fasta = length(fasta_files)

    if num_meta == 2 && num_fasta == 1

        metadata1 = CSV.read(metadata_files[1], DataFrame, delim='\t', header=1)
        metadata2 = CSV.read(metadata_files[2], DataFrame, delim='\t', header=1)
        find_double_candidates(metadata1, metadata2, fasta_files[1])

    elseif num_meta == 3 && num_fasta == 1

        metadata1 = CSV.read(metadata_files[1], DataFrame, delim='\t', header=1)
        metadata2 = CSV.read(metadata_files[2], DataFrame, delim='\t', header=1)
        metadata3 = CSV.read(metadata_files[3], DataFrame, delim='\t', header=1)
        find_double_candidates(metadata1, metadata2, metadata3, fasta_files[1])
        
    elseif num_meta == 3 && num_fasta == 0

        metadata1 = CSV.read(metadata_files[1], DataFrame, delim='\t', header=1)
        metadata2 = CSV.read(metadata_files[2], DataFrame, delim='\t', header=1)
        metadata3 = CSV.read(metadata_files[3], DataFrame, delim='\t', header=1)
        find_double_candidates(metadata1, metadata2, metadata3)

    elseif num_meta == 2 && num_fasta == 0

        metadata1 = CSV.read(metadata_files[1], DataFrame, delim='\t', header=1)
        metadata2 = CSV.read(metadata_files[2], DataFrame, delim='\t', header=1)
        find_double_candidates(metadata1, metadata2)

    else
        throw(ErrorException("No comparable candidate files provided."))
    end
end

# run main
main(metadata_files, fasta_files)
