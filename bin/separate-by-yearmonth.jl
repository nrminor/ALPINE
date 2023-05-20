#!/usr/bin/env julia

# loading necessary packages
using FastaIO, FileIO, Dates

# Check if the input FASTA file is a symlink, and if it is, follow the symlink
const fasta_path = islink(ARGS[1]) ? readlink(ARGS[1]) : ARGS[1]

# defining a function that parses the dates that were added to the defline
function parse_date_from_defline(defline::String)
    components = split(defline, '|')
    date_string = components[end]
    try
        return Date(date_string)
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
            record_date = parse_date_from_defline(name)

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

            # split the defline by the pipe symbol and keep the first split
            name = split(name, '|')[1]

            writeentry(writer, name, seq)
        end
    end
end

# run the function
separate_by_month(fasta_path)