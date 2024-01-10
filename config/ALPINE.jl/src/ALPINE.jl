module ALPINE

# load dependencies
using DelimitedFiles,
    DataFrames,
    CSV,
    Arrow,
    FASTX,
    FastaIO,
    FileIO,
    Dates,
    BioSequences,
    Distances,
    Statistics,
    Pipe,
    CodecZstd,
    CodecZlib,
    FLoops,
    Missings,
    RCall,
    Scratch,
    Statistics,
    StringDistances
import Base.Threads

export distance_matrix,
    set_to_uppercase,
    weight_by_cluster_size,
    date_pango_calls,
    create_rarity_lookup,
    assign_anachron,
    find_anachron_seqs,
    find_double_candidates,
    estimate_prevalence

### FUNCTION(S) TO COMPUTE WEIGHTED DISTANCE MATRICES FOR GENBANK FASTA FILES ###
### ------------------------------------------------------------------------- ###
# replace lowercase n symbols with uppercase Ns
"""
The `vsearch --cluster_fast` algorithm produces some unusual patterns in its
outputs, including lowercase "n" symbols for masked bases. To avoid confusing
algorithms downstream, this function opens the input FASTA provided by the
argument `fasta_filename::String`, sets all its characters to uppercase, and
writes them to a temporary FASTA file with a name supplied by the argument
`temp_filename::String`.
"""
function set_to_uppercase(fasta_filename::String, temp_filename::String)
    open(fasta_filename) do infile
        open(temp_filename, "w") do outfile
            # Read, replace, and write one line at a time
            for line in eachline(infile)
                if startswith(line, ">")
                    println(outfile, line)
                else
                    line = uppercase(line)
                    println(outfile, line)
                end
            end
        end
    end
end

# determine how many sequences are in each cluster to compute weighted distances
"""
For the sequence identified by argument `seq_name::String` in the dataframe `dist_df::DataFrame`
and the `vsearch --cluster_fast` output table provided by argument `cluster_table::DataFrame`,
use the size of that sequence's cluster to "weight" its distances from other clusters. This 
weight, adjusted by `stringency::String`, will penalize large clusters and advantage small
clusters
"""
function weight_by_cluster_size(
    seq_name::String,
    stringency::String,
    dist_df::DataFrame,
    cluster_table::DataFrame,
)

    # filter down to centroid rows only
    centroids = filter(:1 => x -> x == "C", cluster_table)

    # make sure centroids are in the same order as the distance matrix
    centroids = centroids[indexin(names(dist_df), centroids.Column9), :]
    @assert centroids[:, 9] == names(dist_df)

    # extract a vector of cluster sizes that are in the same order as the 
    # sequences in the distance matrix
    all_sizes = centroids[:, 3]

    # filter down to hits only
    hits = filter(:1 => x -> x != "C", cluster_table)
    month_total = nrow(hits) == 0 ? 1 : nrow(hits)

    # find the number of sequences for the accession in question
    cluster_freq = filter(:9 => x -> x == seq_name, centroids)[:, 3][1] / month_total

    # return the weight for this centroid's distance
    weights =
        stringency == "strict" ? ((all_sizes .* -log(cluster_freq)) ./ month_total) :
        ((all_sizes .* (1 - cluster_freq)) ./ month_total)
    @assert length(weights) == nrow(dist_df)
    return weights

end

# compute cluster-size-weighted weighted distance matrices

"""
Compute a pairwise Jaccard distance matrix between all the sequences in the FASTA provided
by argument `temp_filename::String`. Use the information in argument `cluster_table::DataFrame`
and `stringency::String` to weight the distances by sequence ID, and use arguments `count::Int`,
`majority_centroid::String`, and `yearmonth::String` to structure and name output CSV files.
"""
function distance_matrix(
    temp_filename::String,
    cluster_table::DataFrame,
    yearmonth::String,
    stringency::String,
)

    # Collect both names and sequences
    seqs = [seq for (_, seq) in FastaReader(temp_filename)]
    seq_names = [name for (name, _) in FastaReader(temp_filename)]

    # Convert the sequences to BioSequence objects, replacing 'N' characters
    filtered_seqs = [replace(seq, 'N' => '-') for seq in seqs]
    seq_vectors = [LongSequence{DNAAlphabet{4}}(seq) for seq in filtered_seqs]

    # Compute the Hamming distance matrix
    dist_matrix = pairwise(Jaccard(51), seq_vectors, seq_vectors)

    # Convert the distance matrix to a DataFrame
    dist_df = DataFrame(dist_matrix, :auto)
    rename!(dist_df, seq_names)

    # weight distance estimates by relative cluster size
    for seq in names(dist_df)
        col_weights = weight_by_cluster_size(seq, stringency, dist_df, cluster_table)
        dist_df[!, seq] = dist_df[!, seq] .* col_weights
    end

    # Add a column for the sequence names
    dist_df[!, :Sequence_Name] = seq_names

    # Move the Sequence_Name column to the front
    select!(dist_df, :Sequence_Name, :)

    # Write the distance matrix to a CSV file
    CSV.write("$yearmonth-dist-matrix.csv", dist_df)

end

### ------------------------------------------------------------------------- ###


### FUNCTION THAT FINDS THE OVERLAP BETWEEN HIGH-DISTANCE AND ANACHRONISTIC SEQUENCES ###
### --------------------------------------------------------------------------------- ###

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
    prevalence = join((candidate_count / sample_size) * 100, ";")

    return (prevalence, sample_size)

end

### ---------------------------------------------------------------------------------------- ###

end # LongInfectionFinder module
