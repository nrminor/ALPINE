#!/usr/bin/env julia

# loading necessary packages
using FastaIO, FileIO
import Base.Threads

# Check if the input FASTA file is a symlink, and if it is, follow the symlink
const fasta_path = islink(ARGS[1]) ? readlink(ARGS[1]) : ARGS[1]

# create the output file
const out = "filtered-by-n.fasta"

# nest the most intensive steps in a function for faster compilation
function filter_by_n(input_fasta_path::String, output_filename::String)

    touch(output_filename)

    # create a file lock to prevent threads from corrupting the output file
    u = ReentrantLock()

    # creating a loop that goes through sequence records and writes out
    # any sequences that have less than the minimum N count
    FastaWriter(output_filename , "a") do fa
        FastaReader(input_fasta_path) do fr
            @sync for (name, seq) in fr
                Threads.@spawn begin
                    max_n_count = floor(length(seq) * 0.1)
                    n_count = count("N", convert(String, seq))
                    if n_count < max_n_count
                        Threads.lock(u) do
                            writeentry(fa, name, seq)
                        end
                    end
                end
            end
        end
    end

end

# run the function
filter_by_n(fasta_path, out)