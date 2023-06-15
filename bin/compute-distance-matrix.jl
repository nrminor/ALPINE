#!/usr/bin/env julia

# loading necessary packages
using FastaIO, BioSequences, Distances, Statistics, DataFrames, CSV

# parse supplied command line arguments to locate files and set parameters
const fasta_file = islink(ARGS[1]) ? readlink(ARGS[1]) : ARGS[1]
const table_path = islink(ARGS[2]) ? readlink(ARGS[2]) : ARGS[2] # this will be used in the future
const yearmonth = ARGS[3]
const count = parse(Int, ARGS[4])
const majority_centroid = ARGS[5]

# replace lowercase n symbols with uppercase Ns
const tmp = "tmp.fasta"
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
touch(tmp)
set_to_uppercase(fasta_file,tmp)

# determine how many sequences are in each cluster to compute weighted Distances
function weight_by_cluster_size(seq_name::String, cluster_table::DataFrame)

    # filter down to centroid rows only
    centroids = filter(:1 => x -> x == "C", cluster_table)

    # filter down to hits only
    hits = filter(:1 => x -> x != "C", cluster_table)
    if nrow(hits) == 0
        month_total = 1
    else
        month_total = nrow(hits)
    end

    # find the number of rows for the accession in question
    cluster_size = nrow(filter(:9 => x -> x == seq_name, centroids))

    # return the weight for this centroid's distance
    weight = cluster_size / month_total
    return weight

end

# nest the most intensive steps in a function for faster compilation
function distance_matrix(temp_filename::String, cluster_table::DataFrame, yearmonth::String)

    # Collect both names and sequences
    seqs = [seq for (name, seq) in FastaReader(temp_filename)]
    seq_names = [name for (name, seq) in FastaReader(temp_filename)]

    # Convert the sequences to BioSequence objects, replacing 'N' characters
    filtered_seqs = [replace(seq, 'N' => '-') for seq in seqs]
    seq_vectors = [LongSequence{DNAAlphabet{4}}(seq) for seq in filtered_seqs]

    # Compute the Hamming distance matrix
    dist_matrix = pairwise(Hamming(), seq_vectors, seq_vectors)

    # Convert the distance matrix to a DataFrame
    dist_df = DataFrame(dist_matrix, :auto)
    rename!(dist_df, seq_names)

    # weight distance estimates by relative cluster size
    for seq in names(dist_df)
        col_weight = weight_by_cluster_size(seq, cluster_table)
        dist_df[!, seq] = dist_df[!, seq] * col_weight
    end

    # Add a column for the sequence names
    dist_df[!, :Sequence_Name] = seq_names

    # Move the Sequence_Name column to the front
    select!(dist_df, :Sequence_Name, :)

    if count > 2

        # Write the distance matrix to a CSV file
        CSV.write("$yearmonth-dist-matrix.csv", dist_df)
    
    else

        # constrain down to a simple one by one with the non-majority cluster 
        # as the only row
        filter!(:Sequence_Name => !=(majority_centroid), dist_df)
        select!(dist_df, Not(Symbol(dist_df[1,:1])))

        # Write the distance matrix to a CSV file
        CSV.write("$yearmonth-dist-matrix.csv", dist_df)

    end

end

# read in the cluster table
cluster_table = CSV.read(table_path, DataFrame, delim="\t", header=false)

# run the function
distance_matrix(tmp, cluster_table, yearmonth)
