#!/usr/bin/env julia

# loading necessary packages
using LongInfectionFinder, FileIO, DataFrames, Dates, Arrow, Pipe

# locate metadata
const input_table_path = islink(ARGS[1]) ? readlink(ARGS[1]) : ARGS[1]

# validate metadata
valid_meta = validate_metadata(input_table_path)

# Writing validated metadata
Arrow.write("validated-metadata.arrow", valid_meta)
