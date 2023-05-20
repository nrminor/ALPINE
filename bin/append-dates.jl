#!/usr/bin/env julia

# loading necessary packages
using CSV, DataFrames, FastaIO, FileIO
import Base.Threads

# Check if the input metadata file is a symlink, and if it is, follow the symlink
const metadata_path = islink(ARGS[1]) ? readlink(ARGS[1]) : ARGS[1]

# Check if the input FASTA file is a symlink, and if it is, follow the symlink
const fasta_path = islink(ARGS[2]) ? readlink(ARGS[2]) : ARGS[2]

# define output file name
const output = "dated-seqs.fasta"

# nest the most intensive steps in a function for faster compilation
function append_dates(input_tsv_path::String, input_fasta_path::String, output_handle::String)

    # prepare empty output file to append records to
    touch(output_handle)

    # create a file lock to prevent threads from corrupting the output file
    u = ReentrantLock()

    # read in metadata
    metadata_df = CSV.read(input_tsv_path, DataFrame, delim="\t")

    # create lookup of dates and accessions
    accession_to_date = Dict(zip(metadata_df[!,"Accession"], metadata_df[!,"Isolate Collection date"]))

    # add date for each FASTA record
    FastaWriter(output_handle, "a") do fa
        FastaReader(input_fasta_path) do reader 
            @sync for (name, seq) in reader
                Threads.@spawn begin
                    # Parse out the accession number from the defline
                    accession = split(name, ' ')[1]
                    # Look up the collection date for the accession number
                    collection_date = get(accession_to_date, accession, "")
                    # Add the collection date to the defline after a pipe symbol
                    new_name = string(name, "|", collection_date)
                    # Write the modified record to a new FASTA file
                    Threads.lock(u) do
                        writeentry(fa, new_name, seq)
                    end
                end
            end
        end
    end

end

# run the function
append_dates(metadata_path,fasta_path,output)
