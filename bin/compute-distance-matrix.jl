#!/usr/bin/env julia

# loading necessary packages
using FastaIO, BioSequences, Distances, Statistics, DataFrames, CSV

# parse supplied command line arguments to locate files and set parameters
fasta_file = ARGS[1]
if islink(fasta_file)
    fasta_file = readlink(fasta_file)
end
yearmonth = String(ARGS[2])

# Collect both names and sequences
seqs = [seq for (name, seq) in FastaReader(fasta_file)]
names = [name for (name, seq) in FastaReader(fasta_file)]

# Convert the sequences to BioSequence objects
seq_vectors = [LongSequence{DNAAlphabet{5}}(seq) for seq in seqs]

# Compute the Hamming distance matrix
dist_matrix = pairwise(Hamming(), seq_vectors, seq_vectors)

# Find the sequence that is the highest distance, on average, from all the other sequences
avg_dists = mean(dist_matrix, dims=1)
max_avg_dist_index = argmax(avg_dists)[1]
println("The sequence with the highest average distance from all others is sequence number $max_avg_dist_index")

# Convert the distance matrix to a DataFrame
dist_df = DataFrame(dist_matrix)
rename!(dist_df, names)

# Add a column for the sequence names
dist_df[!, :Sequence_Name] = names

# Move the Sequence_Name column to the front
select!(dist_df, :Sequence_Name, :)

# Write the distance matrix to a CSV file
CSV.write("$yearmonth-dist-matrix.csv", dist_df)