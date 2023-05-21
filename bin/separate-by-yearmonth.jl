#!/usr/bin/env julia

# loading necessary packages
using FastaIO, FileIO, Dates

# Check if the input FASTA file is a symlink, and if it is, follow the symlink
const fasta_path = islink(ARGS[1]) ? readlink(ARGS[1]) : ARGS[1]

# locate metadata
const input_tsv_path = islink(ARGS[2]) ? readlink(ARGS[1]) : ARGS[1]

# read in metadata
metadata_df = CSV.read(input_tsv_path, DataFrame, delim="\t")

# create lookup of dates and accessions
accession_to_date = Dict(zip(metadata_df[!,"Accession"], metadata_df[!,"Isolate Collection date"]))

# define a function that accesses the collection date for each record name
function lookup_date(record_name::String, lookup::Dict)

    # Look up the collection date for the accession number
    date = get(lookup, record_name, "")
    try
        return Date(date)
    catch
        return nothing
    end
end

# nest the most intensive steps in a function for faster compilation
function separate_by_month(input_fasta::String)

    # creating a dictionary of FastaWriters for each year_month
    open_writers = Dict{String, FastaWriter}()

    # create a finalizer to close all writers when script ends
    atexit(() -> foreach(close, values(open_writers)))

    # process each sequence in the input FASTA file
    FastaReader(input_fasta) do fr
        for (name, seq) in fr
            record_date = lookup_date(name,accession_to_date)

            if record_date === nothing
                @warn "Skipping record with missing or unparseable date:" name
                continue
            end

            year_month = Dates.format(record_date, "yyyy-mm")

            # get the writer for this year_month, or create a new one if it doesn't exist yet
            writer = get!(open_writers, year_month) do
                output_filename = string(year_month, ".fasta")
                touch(output_filename)
                FastaWriter(output_filename, "a")
            end

            writeentry(writer, name, seq)
        end
    end
end

# run the function
separate_by_month(fasta_path)