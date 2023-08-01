#!/usr/bin/env julia

# loading necessary packages
using LongInfectionFinder, FileIO, DataFrames, Dates, Arrow, Pipe

# locate metadata
const input_table_path = islink(ARGS[1]) ? readlink(ARGS[1]) : ARGS[1]

# retrieve arrow database
arrow_path = represent_in_arrow(input_table_path)

# validate metadata
valid_meta = validate_metadata(arrow_path)

# Writing validated metadata
Arrow.write("validated-metadata.arrow", valid_meta)
