#!/usr/bin/env python3

"""
REPORT EVOLUTIONARILY ADVANCED SEQUENCES FOR ADDITIONAL SCRUTINY
----------------------------------------------------------------

This script collates a large number of clustering, genetic distance, 
and sequence files to identify sequences that are suitable for 
additional review. To do so efficiently, it uses the Polars
library to represent data frames and strictly separates file
reading and writing from CPU-intensive computations. See function
docstrings below for more granular explanations.
----------------------------------------------------------------
"""


# bring in the modules used in this namespace
import os
import argparse
import numpy as np
import polars as pl
import pyarrow
from polars.testing import assert_frame_equal


# define functions:
def quantify_stringency(stringency: str) -> int:
    """
    Define strictness level as a quantile.

    Args:
    stringency: a string supplied by the user or the pipeline specifying
    the level of stringency to apply when choosing which candidates to
    retain.

    Returns:
    Integer representing a quantile to apply out of 1,000. This quantile
    will be used to select a retention threshold base on the true
    distribution of distance scores in the data in question.
    """

    if stringency == "strict":
        strict_quant = 995
    elif stringency == "intermediate":
        strict_quant = 990
    elif stringency == "lenient":
        strict_quant = 980
    else:
        strict_quant = 995

    return strict_quant

def read_metadata_files(metadata_path: str,
                        yearmonths: list,
                        workingdir: str) -> tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame]:
    """
    Metadata file reading is relegated to this function so that other functions
    can be CPU-bound alone instead of interleaving IO and additional computations.
    The dataframes returned by this function have been organized so that it will
    be trivial to use joins with them downstream.

    Args:
    yearmonths: a list of strings that look like "2021-01"
    workingdir: the working directory. By default, it is the current `pwd` working directory.

    Returns:
    A tuple of three polars dataframes, one of which is the full database metadata, one
    of which records the distance scores for each accession, and one of which records 
    clustering metadata for each accession. These metadata include 1) whether an accession
    is a "hit" or a "centroid," 2) the cluster index for the month each accession was 
    collected it, 3) the size of each cluster, and 4) the sequence accessions themselves.
    """

    # Make some empty dataframes to vstack on top of
    dist_scores = pl.DataFrame({
        "Accession": None,
        "Distance Score": None
    }, schema={
        "Accession": pl.Utf8, "Distance Score": pl.Float64
    })
    dist_scores = dist_scores.clear()
    cluster_meta = pl.read_csv(f"{workingdir}/{yearmonths[0]}-clusters.uc",
                                separator="\t", has_header=False, columns=[0,1,2,8],
                                new_columns=["Type", "Index", "Size", "Accession"],
                                n_rows=1)
    cluster_meta = cluster_meta.with_columns(pl.lit("2023-08").alias("Month"))
    cluster_meta = cluster_meta.clear()

    # go through each month of data and bring in its various files
    for yearmonth in yearmonths:

        # sum up distance score for each accession and prep to vstack on
        # the distance score data frame, which we will use in some joins
        # in a later, non-IO-bound function.
        distmat = pl.read_csv(f"{workingdir}/{yearmonth}-dist-matrix.csv")
        distmat = distmat.with_columns(
            Distance_Score=pl.Series(
            [distmat.select(col).sum().item() for col in distmat.columns[1:]]
            )
        ).rename({
            "Sequence_Name": "Accession",
            "Distance_Score": "Distance Score"
        }).select(pl.col(
            ["Accession", "Distance Score"]
        ))

        # Now, bring in clustering metadata and do the same kind of thing
        cluster_table = pl.read_csv(f"{workingdir}/{yearmonth}-clusters.uc",
                                    separator = "\t", has_header=False, columns=[0,1,2,8],
                                    new_columns=["Type", "Index", "Size", "Accession"])
        cluster_table = cluster_table.with_columns(
            Month=pl.repeat(yearmonth, n=cluster_table.shape[0])
        ).filter(
            pl.col("Type") != "S"
        )

        # correct any corrupted rows so the indices can be properly typed
        if cluster_table.dtypes[1] == pl.Utf8:
            cluster_table = cluster_table.filter(
                (pl.col("Type") == "C") |
                (pl.col("Type") == "H") &
                (pl.col("Index") != "*")
            ).with_columns(
                pl.col("Index").cast(pl.Int64)
            )

        # Use an assertion to check for possible silent errors before vstacking and
        # returning the dataframes
        assert_frame_equal(
            left=distmat.select(pl.col("Accession")).unique().sort(by="Accession"),
            right=cluster_table.filter(
            pl.col("Type")=="C"
            ).select(
            pl.col("Accession")).unique().sort(by="Accession"
            )
        )

        # vstack (i.e., append) them onto the growing metadata dataframes
        dist_scores.vstack(distmat, in_place=True)
        cluster_meta.vstack(cluster_table, in_place=True)

    # parse the full database metadata as a data frame to return
    metadata = pl.read_ipc(f"{workingdir}/{metadata_path}", use_pyarrow = True)

    return (metadata, dist_scores, cluster_meta)

# end read_metadata_files def

def collate_metadata(metadata: pl.DataFrame,
                     dist_scores: pl.DataFrame,
                     cluster_meta: pl.DataFrame,
                     stringency: int) -> pl.DataFrame:
    """
    This function takes the dataframes read by read_metadata_files and 
    runs a number of computations on them. The thrust of all these
    computations is to join them together with the whole-database 
    metadata while maximizing helpful information for each candidate
    cluster. Notably, this function handles the CPU-bound tasks that
    are required with the provided metadata. 
    
    Future updates may implement multiprocessing for some of the slower
    steps here. They may also implement LazyFrames instead of dataframes
    to handle very large metadata files.

    Args:
    metadata - a polars dataframe that contains the whole-dataset metadata
    dist_scores - a polars dataframe that contains distance scores for each
                  cluster centroid
    cluster_meta - a polars dataframe that contains metadata from the
                   clustering algorithm.
    """

    # join distance scores onto the cluster metadata, which will add
    # scores to the centroids from each cluster
    cluster_meta = cluster_meta.join(dist_scores, on = "Accession", how="left")

    # merge the month and cluster index columns to make a unique identifier
    cluster_meta = cluster_meta.with_columns([
        pl.format("{}_{}", "Month", "Index").alias("Month Index")
    ])

    # Use a series of Polars expressions to:
    # 1 - replace sequence length entries in "H" entries with the cluster
    # size entries from the associated "C" entries.
    # 2 - Spread the distance scores from each centroid to all the hit
    # sequences in their clusters, and
    # 3 - select only the cluster size, distance, accession, and month
    # columns for joining with the big metadata
    centroid_data = cluster_meta.filter(
        pl.col("Type") == "C"
    ).select([
        "Month Index", "Size", "Distance Score"
    ])
    cluster_meta = cluster_meta.join(
        centroid_data, on="Month Index", how="left", suffix="_right"
    ).select([
        "Accession", "Month Index", "Size_right", "Distance Score_right"
    ]).rename({
        "Size_right": "Cluster Size",
        "Distance Score_right": "Distance Score"
    })

    # filter out any clusters with only one entry (this shouldn't
    # be necessary when run in the context of the pipeline, but
    # we include it here just to be safe.)
    cluster_meta = cluster_meta.groupby("Month Index").agg(
        pl.col("Month Index").count().alias("count")
    ).join(
        cluster_meta, on="Month Index", how="left"
    ).filter(
        pl.col("count") > 1
    ).select(
        cluster_meta.columns
    )

    # join with the big metadata
    high_dist_meta = metadata.join(cluster_meta, on="Accession", how="left")

    # sort the new high distance metadata by distance score
    # in descending order and remove any null rows
    high_dist_meta = high_dist_meta.sort(
        by = "Distance Score", descending=True, nulls_last=True
    ).filter(
        pl.col("Distance Score").is_not_null()
    )

    # And finally, filter down to distances above the retention threshold
    retention_threshold = np.quantile(high_dist_meta['Distance Score'], (stringency / 1000))
    high_dist_meta = high_dist_meta.filter(
        pl.col("Distance Score") >= retention_threshold
    )

    return high_dist_meta

# end collate_metadata def

def main(metadata_path: str, stringency: str, workingdir: str):
    """
    Main function. Main ties together all the above functions if 
    they are being invoked as a script.
    """

    # Retrive a quantile to use to define a retention threshold downstream
    strict_quant = quantify_stringency(stringency)

    # List all files in current directory that end with '-dist-matrix.csv'
    files = [f for f in os.listdir(workingdir) if f.endswith('-dist-matrix.csv')]

    # Remove '-dist-matrix.csv' from each file name
    yearmonths = [f.replace('-dist-matrix.csv', '') for f in files]

    # retrieve all candidate metadata
    metadata, dist_scores, cluster_meta = read_metadata_files(metadata_path, yearmonths, workingdir)

    # pile up all metadata from clustering, distance-calling, and the original database
    high_dist_metadata = collate_metadata(metadata, dist_scores, cluster_meta, strict_quant)

    # write the metadata
    high_dist_metadata.write_csv(f"{workingdir}/high_distance_candidates.tsv", separator="\t")

# end main def

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("metadata",
                        default="filtered-to-geography.arrow",
                        help="The name of the database metadata file to join other files to.")
    parser.add_argument("stringency",
                        help="The chosen level of strictness for retaining candidates.")
    parser.add_argument("workingdir",
                        default=".",
                        help="Directory where the script should search for files.")
    args = parser.parse_args()
    main(args.metadata, args.stringency, args.workingdir)
