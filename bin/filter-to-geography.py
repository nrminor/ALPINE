#!/usr/bin/env python3

"""
usage: filter-to-geography.py [-h] --metadata METADATA [--geography GEOGRAPHY] [--max_date MAX_DATE] [--min_date MIN_DATE]

options:
  -h, --help            show this help message and exit
  --metadata METADATA, -m METADATA
                        The metadata file in Apache Arrow format.
  --geography GEOGRAPHY, -g GEOGRAPHY
                        A substring to use for filtering to a geography of interest, e.g., 'USA'
  --max_date MAX_DATE, -b MAX_DATE
                        The maximum date of interest. Defaults to today.
  --min_date MIN_DATE, -s MIN_DATE
                        The minimum date of interest. Defaults to None.
"""

import argparse
import asyncio
import os
import sys
from datetime import date
from enum import Enum
from inspect import cleandoc
from pathlib import Path
from typing import Optional

import polars as pl
from loguru import logger
from pydantic import validate_call
from pydantic.dataclasses import dataclass
from result import Err, Ok, Result


@dataclass(frozen=True, kw_only=True)
class FilterParams:
    """
    The FilterParams dataclass contains the three relevant
    filters for subsetting the input dataset: geography of
    interest, mininum date, and maximum date within a date
    range of interest. The minimum date is optional. If no
    minimum date is specified, it will be of type `None`.
    """

    geography: str
    min_date: Optional[date]
    max_date: date


class FileType(Enum):
    """
    `FileType` is a sum type/enum containing the file formats
    supported by the script. If the provided file name contains
    a supported extension, the data will be lazily loaded.
    """

    ARROW = (".arrow", pl.scan_ipc)
    PARQUET = (".parquet", pl.scan_parquet)
    CSV = (".csv", pl.scan_csv)
    TSV = (".tsv", lambda path: pl.scan_csv(path, separator="\t"))

    def __init__(self, file_extension, load_function):
        self.file_extension = file_extension
        self.load_function = load_function

    @staticmethod
    def determine_file_type(metadata_path: str):
        """
        This static method returns the file-type variant that
        matches the provided metadata file.
        """
        for file_type in FileType:
            if file_type.file_extension in metadata_path:
                return file_type
        return None


@validate_call
def parse_command_line_args() -> Result[argparse.Namespace, str]:
    """
        Parse command line arguments while passing errors onto main.

    Args:
        `None`

    Returns:
        `Result[argparse.Namespace, str]`: A Result type containing an `argparse` namespace
        or an error message string.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--metadata",
        "-m",
        type=Path,
        required=True,
        help="The metadata file in Apache Arrow format.",
    )
    parser.add_argument(
        "--geography",
        "-g",
        type=str,
        required=False,
        default="USA",
        help="A substring to use for filtering to a geography of interest, e.g., 'USA'",
    )
    parser.add_argument(
        "--max_date",
        "-b",
        type=str,
        required=False,
        default=date.today(),
        help="The maximum date of interest. Defaults to today.",
    )
    parser.add_argument(
        "--min_date",
        "-s",
        type=str,
        required=False,
        default=None,
        help="The minimum date of interest. Defaults to None.",
    )
    args = parser.parse_args()
    if len(vars(args)) < 1:
        return Err("Zero arguments could be parsed.")

    return Ok(args)


@logger.catch(reraise=True)
async def read_metadata(metadata_path: Path) -> Result[pl.LazyFrame, str]:
    """
    Read metadata handles reading the metadata based on its file
    extension, alongside any errors.
    """
    path_str = str(metadata_path)

    # Determine the file type
    file_type = FileType.determine_file_type(path_str)

    if file_type:
        logger.info("{}-formatted metadata detected.", file_type.name)
        return Ok(file_type.load_function(metadata_path))

    return Err(
        cleandoc(
            """
        Could not parse the input metadata file type. Please double check
        that input is CSV, TSV, or Apache Arrow/IPC files.
        """
        )
    )


@logger.opt(lazy=True).catch(reraise=True)
async def filter_metadata(
    meta_lf: pl.LazyFrame, filters: FilterParams
) -> Result[pl.LazyFrame, str]:
    """
        Function `filter_metadata` performs filtering with Polars
        expressions on a Polars LazyFrame query. This enables the
        script to process larger-than-memory datasets. This function
        carefully handles errors with the `result` library and
        double-checks if the necessary columns are provided.

    Args:
        `metadata_path: Path`
        `filters: FilterParams`

    Returns:
        `Result[pl.LazyFrame, str]`
    """

    if "Geographic location" not in meta_lf.columns:
        return Err(
            "Geographic location column not found in provided metadata column header."
        )
    if "Isolate Collection date" not in meta_lf.columns:
        return Err(
            "Isolate Collection date column not found in provided metadata column header."
        )

    if filters.min_date is None:
        logger.info("No minimum date provided.")
        logger.info("Filtering to {} before {}.", filters.geography, filters.max_date)
        return Ok(
            meta_lf.filter(
                pl.col("Geographic location").str.contains(filters.geography)
            ).filter(pl.col("Isolate Collection date") <= pl.lit(filters.max_date))
        )

    logger.info(
        "Filtering to {} between {} and {}.",
        filters.geography,
        filters.min_date,
        filters.max_date,
    )
    return Ok(
        meta_lf.filter(pl.col("Geographic location").str.contains(filters.geography))
        .filter(pl.col("Isolate Collection date") >= pl.lit(filters.min_date))
        .filter(pl.col("Isolate Collection date") <= pl.lit(filters.max_date))
    )


@logger.opt(lazy=True).catch(reraise=True)
async def write_out_accessions() -> Result[None, str]:
    """
        If the above functions complete successfully, the
        function `write_out_accessions()` quickly queries the
        filtered dataset, separates out the accession column,
        and writes it to a simple text file without the header.
    Args:
        None

    Returns:
        `Result[None, str]`
    """

    if not os.path.isfile("filtered-to-geography.arrow"):
        return Err(
            cleandoc(
                """
                The filtered dataset in Apache Arrow format 'filtered-to-geography.arrow' is
                missing in the current working directory.
                """
            )
        )

    filtered_lf = pl.scan_ipc("filtered-to-geography.arrow", memory_map=False)

    if "Accession" not in filtered_lf:
        return Err(
            cleandoc(
                """
                Column names may have been corrupted or otherwise interfered with after filtering;
                column 'Accession' is missing
                """
            )
        )

    (
        filtered_lf.select("Accession")
        .collect()
        .write_csv("accessions.txt", separator="\t", include_header=False)
    )

    return Ok(None)


async def main() -> None:
    """
    Function `main()` controls the flow of data through the above functions.
    """

    # initialize logger
    logger.remove()
    logger.add(sys.stderr, backtrace=True, diagnose=True, colorize=True)

    # parse desired filters out of command line arguments
    args_attempt = parse_command_line_args()
    if isinstance(args_attempt, Err):
        sys.exit("Command line argument parsing failed.")
    args = args_attempt.unwrap()  # pylint: disable=E1111
    metadata_path = args.metadata
    run_min_date = None if args.min_date in ("null", "") else args.min_date
    filters = FilterParams(
        geography=args.geography,
        min_date=run_min_date,
        max_date=args.max_date,
    )

    # Scan the metadata into a LazyFrame and return any errors
    metadata_attempt = await read_metadata(metadata_path)
    if isinstance(metadata_attempt, Err):
        sys.exit(metadata_attempt.unwrap_err())

    # filter the metadata
    lz_attempt = await filter_metadata(metadata_attempt.unwrap(), filters)
    if isinstance(lz_attempt, Err):
        sys.exit(
            f"Metadata filtering failed with the following error:\n{lz_attempt.unwrap_err()}"
        )

    # export to filtered arrow file with compression
    lz_attempt.unwrap().collect().write_ipc(
        "filtered-to-geography.arrow", compression="zstd"
    )

    # separate out accessions
    acc_attempt = await write_out_accessions()
    if isinstance(acc_attempt, Err):
        sys.exit(f"Failed to separate out accessions with:\n{acc_attempt.unwrap_err()}")


if __name__ == "__main__":
    asyncio.run(main())
