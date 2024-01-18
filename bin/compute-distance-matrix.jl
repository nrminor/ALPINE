#!/usr/bin/env -S julia -t auto

# loading necessary packages
using FastaIO, BioSequences, Distances, StringDistances, Statistics, DataFrames, CSV

# parse supplied command line arguments to locate files and set parameters
const FASTA_PATH::String = islink(ARGS[1]) ? readlink(ARGS[1]) : ARGS[1]
const TABLE_PATH::String = islink(ARGS[2]) ? readlink(ARGS[2]) : ARGS[2] # this will be used in the future
const YEARMONTH::String = ARGS[3]
const STRINGENCY::String = ARGS[4]

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

"""
Compute a pairwise Damerau-Levenshtein distance matrix between all the sequences in the FASTA provided
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
    dist_matrix = pairwise(DamerauLevenshtein(), seq_vectors)

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

# run main
main(FASTA_PATH, TABLE_PATH, YEARMONTH, STRINGENCY)
