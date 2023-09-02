#!/usr/bin/env -S julia -t auto

# loading packages
using ALPINE, DelimitedFiles, DataFrames, CSV, Dates, Pipe

# saving command line arguments as static-typed constants
if length(ARGS) < 4
    throw(ArgumentError("At least four arguments are required"))
end
const METADATA_PATH::String = islink(ARGS[1]) ? readlink(ARGS[1]) : ARGS[1]
const MAX_DATE::String = ARGS[2]
const MIN_DATE::String = ARGS[3]
const GEOGRAPHY::String = join(ARGS[4:end], " ")

# run the function based on the provided dates with multiple dispatch
function main(metadata_path::String, geo::String, max_date::String, min_date::String)

    if min_date == "none" && max_date == string(Dates.today())
        filter_metadata_by_geo(metadata_path, geo)
    elseif min_date != "none" && max_date == string(Dates.today())
        filter_metadata_by_geo(metadata_path, geo, min_date)
    elseif min_date != "none" && max_date != string(Dates.today())
        filter_metadata_by_geo(metadata_path, geo, min_date, max_date)
    end

end

# precompile main to eek out a bit more performance
precompile(main, (String,String,String,String));

# run main
main(METADATA_PATH, GEOGRAPHY, MAX_DATE, MIN_DATE)
