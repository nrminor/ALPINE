#!/usr/bin/env julia

# loading packages
using LongInfectionFinder, DelimitedFiles, DataFrames, CSV, Dates

# saving command line arguments as static-typed constants
const metadata_path = islink(ARGS[1]) ? readlink(ARGS[1]) : ARGS[1]
const geo = ARGS[2]
const max_date = ARGS[3]
const min_date = ARGS[4]

# run the function based on the provided dates with multiple dispatch
function main(metadata_path::String, geo::String, max_date::String, min_date::String)

    if min_date == ""
        filter_metadata_by_geo(metadata_path, geo)
    elseif min_date != "" & Date.Dates(max_date) == today()
        filter_metadata_by_geo(metadata_path, geo, min_date)
    elseif min_date != "" & Date.Dates(max_date) != today()
        filter_metadata_by_geo(metadata_path, geo, min_date, max_date)
    end

end
main(metadata_path, geo, max_date, min_date)
