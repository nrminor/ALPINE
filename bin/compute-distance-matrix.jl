#!/usr/bin/env julia

# loading necessary packages
using FastaIO, BioSequences, Distances, Statistics, DataFrames, CSV

# parse supplied command line arguments to locate files and set parameters
const fasta_file = islink(ARGS[1]) ? readlink(ARGS[1]) : ARGS[1]
const yearmonth = ARGS[2]

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

# nest the most intensive steps in a function for faster compilation
function distance_matrix(temp_filename::String, yearmonth::String)

    # Collect both names, sequences, and masked base counts
    seqs = [seq for (name, seq) in FastaReader(temp_filename)]
    counts = [count(x -> x == 'N', seq) for seq in seqs]
    names = [name for (name, seq) in FastaReader(temp_filename)]

    # Convert the sequences to BioSequence objects, filtering out 'N' characters
    filtered_seqs = [replace(seq, 'N' => '-') for seq in seqs]
    seq_vectors = [LongSequence{DNAAlphabet{4}}(seq) for seq in filtered_seqs]

    # Compute the Hamming distance matrix
    dist_matrix = pairwise(Hamming(), seq_vectors, seq_vectors)

    # Find the sequence that is the highest distance, on average, from all the other sequences
    avg_dists = mean(dist_matrix, dims=1)
    max_avg_dist_index = argmax(avg_dists)[1]

    # Convert the distance matrix to a DataFrame
    dist_df = DataFrame(dist_matrix, :auto)
    rename!(dist_df, names)

    # Add a column for the sequence names
    dist_df[!, :Sequence_Name] = names

    # Move the Sequence_Name column to the front
    select!(dist_df, :Sequence_Name, :)

    # Write the distance matrix to a CSV file
    CSV.write("$yearmonth-dist-matrix.csv", dist_df)

end

# run the function
distance_matrix(tmp, yearmonth)