#!/usr/bin/env julia

# loading necessary packages
using FastaIO, FileIO, Dates

# locate fasta
fasta_path = ARGS[1]

# Check if the input FASTA file is a symlink, and if it is, follow the symlink
if islink(fasta_path)
    fasta_path = readlink(fasta_path)
end

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

# defining a function that makes a dictionary of record-yearmonth pairs
function group_records_by_year_month(fr::FastaReader)
    grouped_records = Dict{String, Vector{Tuple{String, String}}}()

    for (name, seq) in fr
        record_date = parse_date_from_defline(name)

        if record_date === nothing
            @warn "Skipping record with missing or unparseable date:" name
            continue
        end

        year_month = Dates.format(record_date, "yyyy-mm")

        if !haskey(grouped_records, year_month)
            grouped_records[year_month] = []
        end

        push!(grouped_records[year_month], (name, seq))
    end

    return grouped_records
end

# defining a function that will make a FASTA for each yearmonth group of records
function write_grouped_records(grouped_records::Dict{String, Vector{Tuple{String, String}}})
    for (year_month, records) in grouped_records
        output_filename = string(year_month, ".fasta")

        # make sure a file of that name is available for writing
        touch(output_filename)

        FastaWriter(output_filename, "a") do fa
            for (name, seq) in records
                writeentry(fa, name, seq)
            end
        end
    end
end

# grouping records
grouped_records = FastaReader(fasta_path) do fr
    group_records_by_year_month(fr)
end

# writing new FASTA for each yearmonth group of records
write_grouped_records(grouped_records)
