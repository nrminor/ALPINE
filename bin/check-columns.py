#!/usr/bin/env python3

"""
FINESSE COLUMN NAMES AND CONVERT TO APACHE ARROW IPC
---------------------------------------------------

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
import sys

import polars
from loguru import logger


def parse_command_line_args():
    """parse command line arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument("metadata_path", help="The path to the large metadata TSV.")
    args = parser.parse_args()
    return args.metadata_path


def main():
    """
    This script uses Polars to read a very large TSV file
    using Apache Arrow in memory-representation. It saves
    additional memory and CPU effort by reading it as a
    LazyFrame instead of a DataFrame. It renames some of the
    columns to conform with scripts that expect Genbank metadata
    downstream. It also ensures that the column of dates is typed
    and stored as date objects. Finally, it writes the validated
    metadata to Apache Arrow IPC format on disk for downstream usage.

    Args (parsed from the command line):
    metadata_path: The path to metadata in TSV format.

    Returns:
    None
    """

    logger.add(sys.stderr, backtrace=True, diagnose=True, colorize=True)

    # parse command line arguments
    metadata_path = parse_command_line_args()

    # Scan the metadata into a LazyFrame
    if ".arrow" in metadata_path:
        logger.opt(lazy=True).debug("Arrow-formatted metadata detected.")
        metadata = polars.scan_ipc(metadata_path)
    elif ".parquet" in metadata_path:
        logger.opt(lazy=True).debug("Parquet-formatted metadata detected.")
        metadata = polars.scan_parquet(metadata_path)
    elif ".csv" in metadata_path:
        logger.opt(lazy=True).debug("CSV-formatted metadata detected.")
        metadata = polars.scan_csv(metadata_path)
    elif ".tsv" in metadata_path:
        logger.opt(lazy=True).debug("TSV-formatted metadata detected.")
        metadata = polars.scan_csv(metadata_path, separator="\t")
    else:
        print("Could not parse the input metadata file type.")
        sys.exit("Please only input CSV, TSV, or Apache Arrow/IPC files.")

    # Run some pseudo-eager evaluations
    if "GC-Content" in metadata.columns:
        logger.opt(lazy=True).debug("GISAID metadata detected.")
        # double check that the expected column names are present
        assert "Virus name" in metadata.columns
        assert "Accession ID" in metadata.columns
        assert "Collection date" in metadata.columns
        assert "Location" in metadata.columns
        assert "Pango lineage" in metadata.columns

        # rename "Accession ID", "Collection date", "Location", and "Pango lineage"
        metadata = metadata.rename(
            {
                "Virus name": "Accession",
                "Accession ID": "EPI_ISL",
                "Collection date": "Isolate Collection date",
                "Location": "Geographic Location",
                "Pango lineage": "Virus Pangolin Classification",
            }
        )

    # Double check the column name for geographic locations
    if "Geographic Location" in metadata.columns:
        metadata = metadata.rename({"Geographic Location": "Geographic location"})

    # Correct date typing
    assert "Isolate Collection date" in metadata.columns
    metadata = metadata.with_columns(
        polars.col("Isolate Collection date")
        .str.strptime(polars.Date, format="%Y-%m-%d", strict=False)
        .alias("Isolate Collection date")
    ).filter(polars.col("Isolate Collection date").is_not_null())

    logger.opt(lazy=True).debug("Dates successfully formatted.")
    logger.opt(lazy=True).debug(
        "Now exporting as ZStandard-compressed Apache Arrow IPC file."
    )

    # evaluate and sink into compressed Arrow file in batches
    metadata.sink_ipc("validated.arrow", compression="zstd")


# end main def

if __name__ == "__main__":
    main()
