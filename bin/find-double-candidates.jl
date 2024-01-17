#!/usr/bin/env -S julia -t auto --sysimage alpine.so

using CSV, DataFrames

# load I/O
const METADATA_FILES::Vector{String} = filter(x->occursin(".tsv",x), readdir(abspath("."), join=true))
const FASTA_FILES::Vector{String} = filter(x->occursin(".fasta",x), readdir(abspath("."), join=true))

"""
Here, "double candidates" are defined as sequences that are both highly evolved,
with a large number of nucleotide substitutions relative to co-occurring sequences,
and are anachronistic, occurring long after a lineage is circulating and prevalent.
This method compares between the output of a highly evolved detector and a single
anachronistic reporter, and then filters a FASTA to any identifiable double candidates.
"""
function find_double_candidates(metadata1::DataFrame, metadata2::DataFrame, seqs::String)

    # use assertions to catch runtime errors more informatively and make
    # assumptions about the data explicit
    @assert ncol(metadata1) > 0
    @assert ncol(metadata2) > 0
    @assert "Accession" in names(metadata1)
    @assert "Accession" in names(metadata2)

    # Get the intersecting accessions
    common_accessions = intersect(metadata1[!, :Accession], metadata2[!, :Accession])

    # Filter metadata to only common accessions
    metadata1_filtered = filter(row -> row[1] in common_accessions, metadata1)
    metadata2_filtered = filter(row -> row[1] in common_accessions, metadata2)

    # Sort both dataframes by Accession
    sort!(metadata1_filtered, :Accession)
    sort!(metadata2_filtered, :Accession)

    if occursin("Anachronicity", names(metadata2_filtered)[end])

        # Add the last column from metadata2 to metadata1
        metadata1_filtered[!, :Anachronicity] = metadata2_filtered[!, end]

        # Write the combined metadata to a TSV file
        CSV.write("double_candidate_metadata.tsv", metadata1_filtered, delim = '\t')

    else

        # Add the last column from metadata2 to metadata1
        metadata2_filtered[!, :Anachronicity] = metadata1_filtered[!, end]

        # Write the combined metadata to a TSV file
        CSV.write("double_candidate_metadata.tsv", metadata2_filtered, delim = '\t')

    end

    # Open FASTA files and perform the filtering operation
    FastaReader(seqs) do fr
        FastaWriter("double_candidates.fasta", "w") do fw
            for (name, seq) in fr
                if split(name, " ")[1] in common_accessions
                    writeentry(fw, name, seq)
                end
            end
        end

    end

end

"""
Here, "double candidates" are defined as sequences that are both highly evolved,
with a large number of nucleotide substitutions relative to co-occurring sequences,
and are anachronistic, occurring long after a lineage is circulating and prevalent.
This method looks for double candidates among a single highly evolved sequence
detector and two anachronistic sequence detectors, and then filters a FASTA to any
identifiable double candidates.
"""
function find_double_candidates(
    metadata1::DataFrame,
    metadata2::DataFrame,
    metadata3::DataFrame,
    seqs::String,
)

    # use assertions to catch runtime errors more informatively and make
    # assumptions about the data explicit
    @assert ncol(metadata1) > 0
    @assert ncol(metadata2) > 0
    @assert ncol(metadata3) > 0
    @assert "Accession" in names(metadata1)
    @assert "Accession" in names(metadata2)
    @assert "Accession" in names(metadata3)

    # join the metadata based on anachronistic accessions
    if nrow(metadata2) >= nrow(metadata3)
        metadata_joined = leftjoin(metadata2, metadata3, on = :Accession)
    else
        metadata_joined = leftjoin(metadata3, metadata2, on = :Accession)
    end

    # Get the intersecting accessions
    common_accessions = intersect(metadata1[!, :Accession], metadata_joined[!, :Accession])

    # Filter metadata to only common accessions
    metadata1_filtered = filter(row -> row[1] in common_accessions, metadata1)
    metadata_joined_filtered = filter(row -> row[1] in common_accessions, metadata_joined)

    # Sort both dataframes by Accession
    sort!(metadata1_filtered, :Accession)
    sort!(metadata_joined_filtered, :Accession)

    if occursin("Anachronicity", names(metadata_joined_filtered)[end])

        # Add the last column from metadata2 to metadata1
        metadata1_filtered[!, :Anachronicity] = metadata_joined_filtered[!, end]

        # Write the combined metadata to a TSV file
        CSV.write("double_candidate_metadata.tsv", metadata1_filtered, delim = '\t')

    else

        # Add the last column from metadata2 to metadata1
        metadata_joined_filtered[!, :Anachronicity] = metadata1_filtered[!, end]

        # Write the combined metadata to a TSV file
        CSV.write("double_candidate_metadata.tsv", metadata_joined_filtered, delim = '\t')

    end

    # Open FASTA files and perform the filtering operation
    FastaReader(seqs) do fr
        FastaWriter("double_candidates.fasta", "w") do fw
            for (name, seq) in fr
                if split(name, " ")[1] in common_accessions
                    writeentry(fw, name, seq)
                end
            end
        end

    end

end

"""
Here, "double candidates" are defined as sequences that are both highly evolved,
with a large number of nucleotide substitutions relative to co-occurring sequences,
and are anachronistic, occurring long after a lineage is circulating and prevalent.
This method looks for double candidates from three anachronistic sequence detectors.
"""
function find_double_candidates(
    metadata1::DataFrame,
    metadata2::DataFrame,
    metadata3::DataFrame,
)

    # use assertions to catch runtime errors more informatively and make
    # assumptions about the data explicit
    @assert ncol(metadata1) > 0
    @assert ncol(metadata2) > 0
    @assert ncol(metadata3) > 0
    @assert "Accession" in names(metadata1)
    @assert "Accession" in names(metadata2)
    @assert "Accession" in names(metadata3)

    # join the metadata based on anachronistic accessions
    if nrow(metadata2) >= nrow(metadata3)
        metadata_joined = leftjoin(metadata2, metadata3, on = :Accession)
    else
        metadata_joined = leftjoin(metadata3, metadata2, on = :Accession)
    end

    # Get the intersecting accessions
    common_accessions = intersect(metadata1[!, :Accession], metadata_joined[!, :Accession])

    # Filter metadata to only common accessions
    metadata1_filtered = filter(row -> row[1] in common_accessions, metadata1)
    metadata_joined_filtered = filter(row -> row[1] in common_accessions, metadata_joined)

    # Sort both dataframes by Accession
    sort!(metadata1_filtered, :Accession)
    sort!(metadata_joined_filtered, :Accession)

    if occursin("Anachronicity", names(metadata_joined_filtered))

        # Add the last column from metadata2 to metadata1
        metadata1_filtered[!, :Anachronicity] = metadata_joined_filtered[!, end]

        # Write the combined metadata to a TSV file
        CSV.write("double_candidate_metadata.tsv", metadata1_filtered, delim = '\t')

    else

        # Add the last column from metadata2 to metadata1
        metadata_joined_filtered[!, :Anachronicity] = metadata1_filtered[!, end]

        # Write the combined metadata to a TSV file
        CSV.write("double_candidate_metadata.tsv", metadata_joined_filtered, delim = '\t')

    end

end

"""
Here, "double candidates" are defined as sequences that are both highly evolved,
with a large number of nucleotide substitutions relative to co-occurring sequences,
and are anachronistic, occurring long after a lineage is circulating and prevalent.
This method looks for double candidates from two anachronistic sequence detectors.
"""
function find_double_candidates(metadata1::DataFrame, metadata2::DataFrame)

    # use assertions to catch runtime errors more informatively
    @assert ncol(metadata1) > 0
    @assert ncol(metadata2) > 0
    @assert "Accession" in names(metadata1)
    @assert "Accession" in names(metadata2)
    @assert "Accession" in names(metadata1)
    @assert "Accession" in names(metadata2)

    # Get the intersecting accessions
    common_accessions = intersect(metadata1[!, :Accession], metadata2[!, :Accession])

    # Filter metadata to only common accessions
    metadata1_filtered = filter(row -> row[1] in common_accessions, metadata1)
    metadata2_filtered = filter(row -> row[1] in common_accessions, metadata2)

    # Sort both dataframes by Accession
    sort!(metadata1_filtered, :Accession)
    sort!(metadata2_filtered, :Accession)

    if occursin("Anachronicity", names(metadata2_filtered))

        # Add the last column from metadata2 to metadata1
        metadata1_filtered[!, :Anachronicity] = metadata2_filtered[!, end]

        # Write the combined metadata to a TSV file
        CSV.write("double_candidate_metadata.tsv", metadata1_filtered, delim = '\t')

    else

        # Add the last column from metadata2 to metadata1
        metadata2_filtered[!, :Anachronicity] = metadata1_filtered[!, end]

        # Write the combined metadata to a TSV file
        CSV.write("double_candidate_metadata.tsv", metadata2_filtered, delim = '\t')

    end

end

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
main(METADATA_FILES, FASTA_FILES)
