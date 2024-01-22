#!/usr/bin/env python3

"""
USE PRE-ASSIGNED PANGO LINEAGES FROM NCBI METADATA TO IDENTIFY ANACHRONISTICS
-----------------------------------------------------------------------------

At each new release of the Pangolin tool, NCBI updates the pango lineage calls
associated with all SARS-CoV-2 sequences. This script uses these sequences,
along with a file reporting when each lineage was designated in Pangolin, to
identify lineages that are collected long after a lineage was prevalent. This
is the fastest branch of the advanced-virus-finder workflow, as it does not
involve running Pangolin or clustering.
"""

# make necessary modules available
import argparse
import multiprocessing
from pathlib import Path
import sys
from functools import lru_cache
from typing import Optional, Tuple

import numpy
import pandas
import pyarrow
from loguru import logger
from pydantic import validate_call


@validate_call
def parse_command_line_args() -> Tuple[Path, Path, int, int]:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "metadata_path", type=Path, help="The path to the metadata file."
    )
    parser.add_argument(
        "dates_path", type=Path, help="The path to the pango lineage dates file."
    )
    parser.add_argument(
        "infection_cutoff", type=int, help="The minimum infection duration in days."
    )
    parser.add_argument(
        "cores", type=int, help="The number of cores to use for parallel processing."
    )
    args = parser.parse_args()
    return args.metadata_path, args.dates_path, args.infection_cutoff, args.cores


@logger.catch
@lru_cache
def add_designation_date(
    lineage: str, dates: pandas.DataFrame
) -> Optional[numpy.datetime64]:
    """
    This function collects a pangolin lineage and returns the date
    it was designated in pangolin.

    Args:
    lineage: The current lineage of interest
    dates: a pandas dataframe of lineage designation dates

    Returns:
    Union[numpy.datetime64, None] of the designation date
    """

    date_values = numpy.array(
        dates.loc[dates["lineage"] == lineage, "designation_date"].values
    )
    return date_values[0] if date_values.size > 0 else None


@logger.catch
@lru_cache
def add_infection_duration(
    i: int, ncbi_dates: list, desig_dates: list
) -> Optional[int]:
    """
    This function takes the difference of the collection date
    and associated lineage designation date to identify whether
    a given sample is anachronistic.

    Args:
    i: Integer referring to the current row index.

    Returns:
    Union[int, None] representing the anachronicity
    """

    if pandas.isna(ncbi_dates[i]) or pandas.isna(desig_dates[i]):
        return None

    anachronicity = (ncbi_dates[i] - desig_dates[i]).days
    return anachronicity


@logger.opt(lazy=True).catch
def main():
    """
    This script reads in filtered metadata and pango lineage dates,
    and then identifies long infection candidates.

    Args (parsed from the command line):
    metadata_path: The path to the metadata file.
    dates_path: The path to the pango lineage dates file.
    infection_cutoff: The minimum infection duration in days.
    cores: The number of cores to use for parallel processing.

    Returns:
    None
    """

    logger.add(sys.stderr, backtrace=True, diagnose=True, colorize=True)

    # parse command line arguments
    metadata_path, dates_path, infection_cutoff, cores = parse_command_line_args()

    # Memory-map arrow IPC-formatted metadata
    with pyarrow.memory_map(str(metadata_path), "r") as source:
        metadata_arrays = pyarrow.ipc.open_file(source).read_all()
        metadata = metadata_arrays.to_pandas()

        # read in dates
        dates = pandas.read_csv(dates_path, na_values=["", "NA"])

        # make sure the expected columns are present
        assert "Isolate Collection date" in metadata.columns

        # Ensure dates are properly formatted
        metadata["Isolate Collection date"] = pandas.to_datetime(
            metadata["Isolate Collection date"], errors="coerce"
        )
        dates.designation_date = pandas.to_datetime(
            dates.designation_date, errors="ignore", format="%Y-%m-%d"
        ).fillna("2021-02-18")

        # Parallelize looping through filtered metadata to merge designation dates
        metadata["lineage_designation"] = numpy.NAN
        metadata["infection_duration"] = 0

        # prepare the arguments as a list of tuples
        args1 = [
            (lineage, dates) for lineage in metadata["Virus Pangolin Classification"]
        ]

        # iterate through rows in parallel
        with multiprocessing.Pool(cores) as p:
            metadata["lineage_designation"] = p.starmap(add_designation_date, args1)
            ncbi_dates = metadata["Isolate Collection date"]
            desig_dates = metadata["lineage_designation"]
            args2 = [(i, ncbi_dates, desig_dates) for i in range(len(ncbi_dates))]
            metadata["infection_duration"] = p.starmap(add_infection_duration, args2)

        # Filter down to long infection candidates
        long_infections = metadata.loc[
            (metadata["infection_duration"] >= infection_cutoff)
            & ~metadata["Accession"].isna()
        ].sort_values("infection_duration", ascending=False)

        # Write anachronistic candidates to TSV
        long_infections.to_csv(
            "anachronistic_metadata_only_candidates.tsv", index=False, sep="\t"
        )


if __name__ == "__main__":
    main()
