#!/usr/bin/env julia

# loading necessary packages
using LongInfectionFinder, FileIO, DataFrames, Arrow

# locate metadata
const input_table_path = islink(ARGS[1]) ? readlink(ARGS[1]) : ARGS[1]

# retrieve arrow database
metadata_df = represent_in_arrow(input_table_path)

# validate metadata
valid_meta = validate_metadata(metadata_df)

# Writing validated metadata
Arrow.write("validated-metadata.arrow", valid_meta)
