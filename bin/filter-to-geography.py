#!/usr/bin/env python3

"""
usage: filter-to-geography.py [-h] --metadata METADATA [--geography GEOGRAPHY] [--max_date MAX_DATE]
                              [--min_date MIN_DATE] [--and_include AND_INCLUDE] [--or_include OR_INCLUDE]

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
  --and_include AND_INCLUDE, -a AND_INCLUDE
                        Accessions to include among those filtered by date and geography
  --or_include OR_INCLUDE, -o OR_INCLUDE
                        Accessions to include in addition to those filtered by date and geography
"""

import argparse
import asyncio
import os
import sys
from dataclasses import dataclass
from datetime import date
from enum import Enum
from inspect import cleandoc
from pathlib import Path
from typing import Optional

import polars as pl
from loguru import logger
from pydantic import validate_call
from pydantic.dataclasses import dataclass as pyd_dataclass
from result import Err, Ok, Result


@pyd_dataclass(frozen=True, kw_only=True)
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
    max_date: Optional[date]


@pyd_dataclass(frozen=True, kw_only=True)
class ExplicitInclusions:
    """
    The ExplicitInclusions dataclass contains two optional
    paths to so-called "and inclusions" and "or inclusions,"
    as specified in the command line.
    """

    and_inclusions: Optional[Path]
    or_inclusions: Optional[Path]


@dataclass
class InclusionDFs:
    """
    The InclusionDFs dataclass contains two optional data
    frames associated with any provided exclusion paths.
    """

    and_lf: Optional[pl.LazyFrame]
    or_lf: Optional[pl.LazyFrame]


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
    parser.add_argument(
        "--and_include",
        "-a",
        type=Path,
        required=False,
        default=None,
        help="Accessions to include among those filtered by date and geography",
    )
    parser.add_argument(
        "--or_include",
        "-o",
        type=Path,
        required=False,
        default=None,
        help="Accessions to include in addition to those filtered by date and geography",
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


@logger.catch(reraise=True)
async def read_inclusions(inclusions: ExplicitInclusions) -> Result[InclusionDFs, str]:
    """
    Read inclusion list handles, alongside any errors.
    """

    # initialize lazyframe option types
    and_lf: Optional[pl.LazyFrame]
    or_lf: Optional[pl.LazyFrame]

    # check if inclusions lists are provided and scan them with Polars
    if inclusions.and_inclusions is not None:
        and_lf = pl.scan_csv(
            inclusions.and_inclusions,
            has_header=False,
            separator="\t",
            new_columns=["Accession"],
        )
    else:
        and_lf = None
    if inclusions.or_inclusions is not None:
        or_lf = pl.scan_csv(
            inclusions.or_inclusions,
            has_header=False,
            separator="\t",
            new_columns=["Accession"],
        )
    else:
        or_lf = None

    # do some error handling
    if and_lf is not None and and_lf.collect().shape[0] == 0:
        return Err("And inclusion file that exists was specified, but it is empty.")
    if or_lf is not None and or_lf.collect().shape[0] == 0:
        return Err("OrAnd inclusion file that exists was specified, but it is empty.")

    return Ok(
        InclusionDFs(
            and_lf=and_lf,
            or_lf=or_lf,
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


async def handle_explicit_inclusions(
    input_lf: pl.LazyFrame, meta_lf: pl.LazyFrame, includes: InclusionDFs
) -> Result[pl.LazyFrame, str]:
    """
    `handle_explicit_inclusions()` does just that. It handles any provided
    "and" and "or" inclusions lists by transforming the filtered metadata
    to include those accessions. Note that this step does currently involve
    mutating and executing the filtering query.
    """

    # start with any provided "and inclusions", where the user chooses to filter
    # with the provided filter parameters and with a provided pre-curated accession list.
    # Under these conditions, and accession must be on the include list *and* match
    # the provided filters
    if includes.and_lf is not None:
        logger.info("Filtering metadata to provided 'and inclusions'.")
        meta_lf = includes.and_lf.join(meta_lf, on="Accession", how="left")

    # next, handle any "or inclusions", where the filtered metadata may contain any
    # accessions that match the provided filters *or* be in the provided list
    if includes.or_lf is not None:
        logger.info("Filtering metadata to provided 'or inclusions'.")
        or_df = includes.or_lf.join(input_lf, on="Accession", how="left").collect()
        meta_lf = meta_lf.collect().vstack(or_df).lazy()

    return Ok(meta_lf)


@logger.opt(lazy=True).catch(reraise=True)
async def write_out_accessions(filtered_lf: pl.LazyFrame) -> Result[None, str]:
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

    # unpack the metadata
    metadata_path = args.metadata
    assert os.path.isfile(metadata_path), "Provided metadata file does not exist."

    # unpack any inclusion lists
    if args.and_include is not None:
        assert os.path.isfile(
            args.and_inclusions
        ), "Provided and inclusion list file not found."
    if args.or_include is not None:
        assert os.path.isfile(
            args.or_inclusions
        ), "Provided or inclusion list file not found."
    inclusions = ExplicitInclusions(
        and_inclusions=args.and_inclusions,
        or_inclusions=args.or_inclusions,
    )

    # re-pack filter parameters provided in the command line
    filters = FilterParams(
        geography=args.geography,
        min_date=None if args.min_date in ("null", "", None) else args.min_date,
        max_date=date.today() if args.max_date in (None, "null", "") else args.max_date,
    )

    # Scan the metadata into a LazyFrame and return any errors
    logger.info("Reading metadata.")
    metadata_attempt = await read_metadata(metadata_path)
    if isinstance(metadata_attempt, Err):
        sys.exit(metadata_attempt.unwrap_err())

    # read any inclusion lists
    logger.info("Reading inclusion lists, if provided.")
    include_attempt = await read_inclusions(inclusions)
    if isinstance(include_attempt, Err):
        sys.exit(
            f"Parsing of include files failed with the following error:\n{include_attempt.unwrap_err()}"
        )

    # filter the metadata
    logger.info("Filtering metadata with provided filters.")
    lz_attempt = await filter_metadata(metadata_attempt.unwrap(), filters)
    if isinstance(lz_attempt, Err):
        sys.exit(
            f"Metadata filtering failed with the following error:\n{lz_attempt.unwrap_err()}"
        )

    # modify the polars query to handle inclusion lists
    final_attempt = await handle_explicit_inclusions(
        metadata_attempt.unwrap(), lz_attempt.unwrap(), include_attempt.unwrap()
    )
    if isinstance(final_attempt, Err):
        sys.exit(
            f"Transforming with include lists failed with the following error:\n{final_attempt.unwrap_err()}"
        )

    # write out accessions
    logger.info("Writing out simple list of accessions as text file.")
    acc_attempt = await write_out_accessions(final_attempt.unwrap())
    if isinstance(acc_attempt, Err):
        sys.exit(f"Failed to separate out accessions with:\n{acc_attempt.unwrap_err()}")

    # export to filtered arrow file with compression
    logger.info("Writing filtered metadata into a compressed Arrow IPC file.")
    lz_attempt.unwrap().collect().write_ipc(
        "filtered-to-geography.arrow", compression="zstd"
    )


if __name__ == "__main__":
    asyncio.run(main())
