#!/usr/bin/env python3

# make necessary modules available
import argparse
from typing import Optional
import multiprocessing as mp
import pandas as pd
import numpy as np
# import pyarrow as arrow

# Defining functions for multiprocessing
def add_designation_date(lineage: str, dates: pd.DataFrame) -> Optional[np.datetime64]:
    """
    This function collects a pangolin lineage and returns the date
    it was designated in pangolin.

    Args:
    lineage: The current lineage of interest
    dates: a pandas dataframe of lineage designation dates

    Returns:
    Union[np.datetime64, None] of the designation date
    """

    date_values = dates.loc[dates['lineage'] == lineage, 'designation_date'].values
    return date_values[0] if date_values.size > 0 else None
# end add_designation_date()

def add_infection_duration(i: int, ncbi_dates: list, desig_dates: list) -> Optional[int]:
    """
    This function takes the difference of the collection date
    and associated lineage designation date to identify whether
    a given sample is anachronistic.

    Args:
    i: Integer referring to the current row index.

    Returns:
    Union[int, None] representing the anachronicity
    """

    if pd.isna(ncbi_dates[i]) or pd.isna(desig_dates[i]):
        return None

    anachronicity = (ncbi_dates[i] - desig_dates[i]).days
    return anachronicity
# end add_infection_duration()

# define main function
def main(metadata_path: str, dates_path: str, infection_cutoff: int, cores: int):
    """
    This script reads in filtered metadata and pango lineage dates,
    and then identifies long infection candidates.

    Args:
    metadata_path: The path to the metadata file.
    dates_path: The path to the pango lineage dates file.
    infection_cutoff: The minimum infection duration in days.
    cores: The number of cores to use for parallel processing.

    Returns:
    None
    """

    # read in data
    metadata = pd.read_csv(metadata_path,
                           delimiter='\t',
                           na_values=["", "NA"],
                           dtype={"Virus Pangolin Classification": str})

    # read in dates
    dates = pd.read_csv(dates_path,
                        na_values=["", "NA"])
    
    # Ensure dates are properly formatted
    metadata["Isolate Collection date"] = pd.to_datetime(metadata["Isolate Collection date"], errors='coerce')
    dates.designation_date = pd.to_datetime(dates.designation_date, errors="ignore", format="%Y-%m-%d").fillna("2021-02-18")

    # Parallelize looping through filtered metadata to merge designation dates
    metadata["lineage_designation"] = np.NAN
    metadata["infection_duration"] = 0

    # prepare the arguments as a list of tuples
    args1 = [(lineage, dates) for lineage in metadata['Virus Pangolin Classification']]

    # iterate through rows in parallel
    with mp.Pool(cores) as p:
        metadata['lineage_designation'] = p.starmap(add_designation_date, args1)
        ncbi_dates = metadata['Isolate Collection date']
        desig_dates = metadata['lineage_designation']
        args2 = [(i, ncbi_dates, desig_dates) for i in range(len(ncbi_dates))]
        metadata['infection_duration'] = p.starmap(add_infection_duration, args2)
    
    # Filter down to long infection candidates
    long_infections = metadata.loc[(metadata["infection_duration"] >= infection_cutoff) & ~metadata["Accession"].isna()].sort_values("infection_duration", ascending=False)

    # Write anachronistic candidates to TSV
    long_infections.to_csv("anachronistic_metadata_only_candidates.tsv",
                           index=False,
                           sep="\t")
# end main function def

# run the script if not imported as a module
if __name__ == "__main__":

    # parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("metadata_path",
                        help="The path to the metadata file.")
    parser.add_argument("dates_path",
                        help="The path to the pango lineage dates file.")
    parser.add_argument("infection_cutoff",
                        type=int,
                        help="The minimum infection duration in days.")
    parser.add_argument("cores",
                        type=int,
                        help="The number of cores to use for parallel processing.")
    args = parser.parse_args()

    # run main
    main(args.metadata_path, args.dates_path, args.infection_cutoff, args.cores)
