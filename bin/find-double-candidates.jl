#!/usr/bin/env julia

# load pipeline module
using LongInfectionFinder

# load I/O
const metadata_files = filter(x->occursin(".tsv",x), readdir(abspath("."), join=true))
const fasta_files = filter(x->occursin(".fasta",x), readdir(abspath("."), join=true))

# use multiple dispatch to collect double candidates based on the 
# number of available files
if length(metadata_files) == 2 && length(fasta_files) == 1

    find_double_candidates(metadata_files[1], metadata_files[2], fasta_files[1])

elseif length(metadata_files) == 3 && length(fasta_files) == 1

    find_double_candidates(metadata_files[1], metadata_files[2], metadata_files[3], fasta_files[1])
    
elseif length(metadata_files) == 3 && length(fasta_files) == 0

    find_double_candidates(metadata_files[1], metadata_files[2], metadata_files[3])

elseif length(metadata_files) == 2 && length(fasta_files) == 0

    find_double_candidates(metadata_files[1], metadata_files[2])

else
    throw(ErrorException("No comparable candidate files provided."))
end
