#!/usr/bin/env julia

# loading necessary packages
using FastaIO, FileIO

# saving command line arguments supplied by nextflow
fasta_path = ARGS[1]
output_filename = ARGS[2]

# Check if the input FASTA file is a symlink, and if it is, follow the symlink
if islink(fasta_path)
    fasta_path = readlink(fasta_path)
end

# create the output file
touch(output_filename)

# iterate through each FASTA record
FastaWriter(output_filename) do fa
    # creating a loop that goes through sequence records, replaces '-' with 'N',
    # and writes the modified records to the output file
    FastaReader(fasta_path) do fr
        for (name, seq) in fr
            # Replace '-' symbols with 'N' letters in the sequence
            modified_seq = replace(seq, '-' => 'N')
            
            # Write the entry to the output file
            writeentry(fa, name, modified_seq)
        end
    end
end