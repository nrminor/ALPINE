#!/usr/bin/env python3

"""
FINESS COLUMN NAMES AND CONVERT TO APACHE ARROW IPC
--------------------------------------------------

In the advanced finder pipeline, this script is run as soon
as a) the metadata is supplied by the user, or b) the 
decompressed metadata is extracted from the NCBI datasets 
download. When run, the script double-checks typing for all
columns but especially for collection date, which is crucial
to parse at multiple junctures downstream. Is uses the Polars
data frame library to quickly read and efficiently represent
the metadata in memory. The output metadata in ZStandard-
compressed Arrow IPC format, will have the same benefits
for future scripts.
"""

import argparse
import polars as pl

def main(metadata_path: str):
    """
    This script reads a very large TSV file using Apache Arrow
    in memory-representation. If the TSV is GISAID metadata, it
    renames some of the columns to conform with scripts that 
    expect Genbank metadata downstream. It also ensures that 
    the column of dates is typed and stored as date objects.
    Finally, it writes the validated metadata to Apache Arrow
    IPC format on disk for downstream usage.

    Args:
    metadata_path: The path to metadata in TSV format.

    Returns:
    None
    """

    metadata = pl.read_csv(metadata_path, separator="\t")

    if "GC-Content" in metadata.columns:

        # rename "Accession ID", "Collection date", "Location", and "Pango lineage"
        metadata = metadata.rename({
            "Virus name": "Accession",
            "Accession ID": "EPI_ISL",
            "Collection date": "Isolate Collection date",
            "Location": "Geographic Location",
            "Pango lineage": "Virus Pangolin Classification"
        })

    # Double check the column name for geographic locations
    if "Geographic location" in metadata.columns:
        metadata.rename({"Geographic location": "Geographic Location"})

    # make sure date column is typed as dates
    metadata = metadata.with_columns(
        metadata.select(pl.col("Isolate Collection date").str.strptime(pl.Date, format="%Y-%m-%d", strict=False))
    )
    metadata = metadata.drop_nulls("Isolate Collection date")

    # convert to arrow IPC file on disk
    metadata.write_ipc("validated-metadata.arrow", compression="zstd")

# end main def

if __name__ == "__main__":

    # parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("metadata_path", help="The path to the large metadata TSV.")
    args = parser.parse_args()

    # run main
    main(args.metadata_path)
