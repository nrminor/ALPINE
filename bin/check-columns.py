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
import asyncio
import sys
from enum import Enum, auto
from pathlib import Path

import polars as pl
from loguru import logger


class FileType(Enum):
    ARROW = auto()
    PARQUET = auto()
    CSV = auto()
    TSV = auto()
    UNKNOWN = auto()


def parse_command_line_args():
    """parse command line arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "metadata_path", type=Path, help="The path to the large metadata TSV."
    )
    args = parser.parse_args()
    return args.metadata_path


async def _determine_file_type(metadata_path: Path) -> FileType:
    path_str = str(metadata_path)

    if ".arrow" in path_str:
        return FileType.ARROW
    elif ".parquet" in path_str:
        return FileType.PARQUET
    elif ".csv" in path_str:
        return FileType.CSV
    elif ".tsv" in path_str:
        return FileType.TSV

    return FileType.UNKNOWN


@logger.catch
async def _read_metadata(metadata_path: Path) -> pl.LazyFrame:
    filetype = await _determine_file_type(metadata_path)

    # Scan the metadata into a LazyFrame
    if filetype == FileType.ARROW:
        logger.opt(lazy=True).debug("Arrow-formatted metadata detected.")
        metadata = pl.scan_ipc(metadata_path)
    elif filetype == FileType.PARQUET:
        logger.opt(lazy=True).debug("Parquet-formatted metadata detected.")
        metadata = pl.scan_parquet(metadata_path)
    elif filetype == FileType.CSV:
        logger.opt(lazy=True).debug("CSV-formatted metadata detected.")
        metadata = pl.scan_csv(metadata_path)
    elif filetype == FileType.TSV:
        logger.opt(lazy=True).debug("TSV-formatted metadata detected.")
        metadata = pl.scan_csv(metadata_path, separator="\t")
    else:
        print("Could not parse the input metadata file type.")
        sys.exit("Please only input CSV, TSV, or Apache Arrow/IPC files.")

    return metadata


@logger.opt(lazy=True).catch
async def reconcile_columns(metadata: pl.LazyFrame) -> pl.LazyFrame:
    """
    Reconcile column names between provided GISAID data and NCBI data
    so that the steps downstream can assume the required information is
    under the column headers it expects.
    """

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
        pl.col("Isolate Collection date")
        .str.strptime(pl.Date, format="%Y-%m-%d", strict=False)
        .alias("Isolate Collection date")
    ).filter(pl.col("Isolate Collection date").is_not_null())

    logger.opt(lazy=True).debug("Dates successfully formatted.")
    logger.opt(lazy=True).debug(
        "Now exporting as ZStandard-compressed Apache Arrow IPC file."
    )

    return metadata


async def main():
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
    metadata = await _read_metadata(metadata_path)

    new_cols_lf = await reconcile_columns(metadata)

    # evaluate and sink into compressed Arrow file in batches
    new_cols_lf.sink_ipc("validated.arrow", compression="zstd")


if __name__ == "__main__":
    asyncio.run(main())
