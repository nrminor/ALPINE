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
import argparse
import asyncio
import os
import subprocess
import sys
from dataclasses import dataclass
from glob import glob
from typing import Tuple

import numpy
import polars as pl
import seaborn  # type: ignore
from loguru import logger
from matplotlib import pyplot
from pydantic import validate_call


@validate_call
def parse_command_line_args() -> Tuple[str, str, str, str]:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--sequences",
        "-s",
        default="filtered-to-geography.fasta",
        type=str,
        help="The name of the database FASTA file to filter for candidates.",
    )
    parser.add_argument(
        "--metadata",
        "-m",
        default="filtered-to-geography.arrow",
        type=str,
        help="The name of the database metadata file to join other files to.",
    )
    parser.add_argument(
        "--stringency",
        "-str",
        type=str,
        required=True,
        help="The chosen level of strictness for retaining candidates.",
    )
    parser.add_argument(
        "--workingdir",
        "-wd",
        default=".",
        type=str,
        help="Directory where the script should search for files.",
    )
    args = parser.parse_args()
    return args.metadata, args.sequences, args.stringency, args.workingdir


@logger.catch(reraise=True)
@validate_call
async def quantify_stringency(stringency: str) -> int:
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

    logger.info(
        "Stringency of {} quantified to the {}th quantile.", stringency, strict_quant
    )

    return strict_quant


@logger.opt(lazy=True).catch(reraise=True)
async def visualize_distance_scores(metadata: pl.DataFrame, threshold: float) -> None:
    """
    This function uses Matplotlib and Seaborn to visualize the distribution
    of distance scores in the provided metadata. It also plots the retention
    threshold that was computed based on the data. The graphic is then written
    to a PDF file called "distance_score_distribution.pdf".

    Args:
    metadata - a Polars DataFrame with a column called "Distance Score"
    threshold - a floating point value specifying the distance score retention threshold

    Returns:
    None
    """

    # retrieve the distance scores from the metadata dateframe
    metadata_pd = metadata.to_pandas()
    assert "Distance Score" in metadata_pd.columns
    assert metadata_pd["Distance Score"].nunique() > 1

    # Set the style and size of the plot and generate the histogram
    pyplot.figure(figsize=(7, 5.5))
    seaborn.set_style("whitegrid")
    seaborn.histplot(
        metadata_pd["Distance Score"], kde=True, color="lightblue", element="step"
    )
    max_count = max(pyplot.hist(metadata_pd["Distance Score"], bins=10, alpha=0)[0])

    # Add a vertical line at the retention threshold
    pyplot.axvline(x=threshold, color="red", lw=3)

    # Add a text label for the retention threshold
    pyplot.text(threshold + 20, max_count / 2, f"Retention Threshold:\n{threshold}")

    # Add labels
    pyplot.xlabel("Distance Score")
    pyplot.ylabel("Frequency")
    pyplot.title("Frequency Distribution of Nucleotide Distances")

    # Save the plot to a file
    pyplot.savefig("distance_score_distribution.pdf")


@dataclass
class ClusterTables:
    """
    In-memory tables of information for all clusters/accessions, depending
    on whether preclustering was performed
    """

    full_meta: pl.DataFrame
    dist_scores: pl.DataFrame
    cluster_meta: pl.DataFrame
    whether_preclust: bool


@logger.opt(lazy=True).catch(reraise=True)
@validate_call
async def read_metadata_files(
    metadata_filename: str, yearmonths: list, workingdir: str
) -> ClusterTables:
    """
    Metadata file reading is relegated to this function so that other functions
    can be CPU-bound alone instead of interleaving IO and additional computations.
    The dataframes returned by this function have been organized so that it will
    be trivial to use joins with them downstream.

    Args:
    yearmonths: a list of strings that look like "2021-01"
    workingdir: the working directory. By default, it is the current `pwd` working directory.

    Returns:
    A tuple of three pl.dataframes, one of which is the full database metadata, one
    of which records the distance scores for each accession, and one of which records
    clustering metadata for each accession. These metadata include 1) whether an accession
    is a "hit" or a "centroid," 2) the cluster index for the month each accession was
    collected it, 3) the size of each cluster, and 4) the sequence accessions themselves.
    """

    # start by checking if there are any VSEARCH clustering metadata files
    precluster = len(glob("*.uc")) > 0
    logger.info("Whether preclustering was performed: {}", precluster)

    # prepare to compile distances with an empty dataframe
    dist_scores = pl.DataFrame(
        {"Accession": None, "Distance Score": None},
        schema={"Accession": pl.Utf8, "Distance Score": pl.Float64},
    ).clear()

    # read full input metadata
    metadata = pl.read_ipc(f"{workingdir}/{metadata_filename}", use_pyarrow=True)

    # conditionally prepare to compile VSEARCH cluster metadata
    if precluster:
        cluster_meta = (
            pl.read_csv(
                f"{workingdir}/{yearmonths[0]}-clusters.uc",
                separator="\t",
                has_header=False,
                columns=[0, 1, 2, 8],
                new_columns=["Type", "Index", "Size", "Accession"],
                n_rows=1,
            )
            .with_columns(pl.lit("2023-08").alias("Month"))
            .clear()
        )
    else:
        cluster_meta = metadata.select("Accession")

    # go through each month of data and bring in its various files
    for yearmonth in yearmonths:
        distmat = pl.read_csv(f"{workingdir}/{yearmonth}-dist-matrix.csv")

        # sum up distance score for each accession and prep to vstack on
        # the distance score data frame, which we will use in some joins
        # in a later, non-IO-bound function.
        assert "Sequence_Name" in distmat.columns
        samples = [sample for sample in distmat.columns if sample != "Sequence_Name"]
        sample_size = len(samples)
        sums = [distmat[sample].sum() for sample in samples]
        logger.info(
            "Found {} sequences for {}, with a median distance score of {}",
            sample_size,
            yearmonth,
            (numpy.median(sums) / sample_size),
        )
        distmat = distmat.with_columns(pl.Series(sums).alias("Distance Score")).rename(
            {"Sequence_Name": "Accession"}
        )

        # continue if preclustering was not performed
        if not precluster:
            distmat = (
                distmat.with_columns(pl.lit(sample_size).alias("Sample Size"))
                .with_columns(
                    (pl.col("Distance Score") / pl.col("Sample Size")).alias(
                        "Distance Score"
                    ),
                )
                .select(pl.col(["Accession", "Distance Score"]))
            )

            dist_scores.vstack(distmat, in_place=True)
            continue

        dist_scores.vstack(
            distmat.select(pl.col(["Accession", "Distance Score"])), in_place=True
        )

        # Now, bring in clustering metadata and do the same kind of thing
        cluster_table = (
            pl.read_csv(
                f"{workingdir}/{yearmonth}-clusters.uc",
                separator="\t",
                has_header=False,
                columns=[0, 1, 2, 8],
                new_columns=["Type", "Index", "Size", "Accession"],
            )
            .with_columns(pl.lit(yearmonth).alias("Month"))
            .filter(pl.col("Type") != "S")
        )

        # correct any corrupted rows so the indices can be properly typed
        if cluster_table.dtypes[1] == pl.Utf8:
            cluster_table = cluster_table.filter(
                (pl.col("Type") == "C")
                | (pl.col("Type") == "H") & (pl.col("Index") != "*")
            ).with_columns(pl.col("Index").cast(pl.Int64))

        # vstack (i.e., append) them onto the growing metadata dataframes
        cluster_meta.vstack(cluster_table, in_place=True)

    return ClusterTables(
        full_meta=metadata,
        dist_scores=dist_scores,
        cluster_meta=cluster_meta,
        whether_preclust=precluster,
    )


@logger.opt(lazy=True).catch(reraise=True)
async def collate_metadata(
    cluster_tables: ClusterTables,
    stringency: int,
) -> tuple[pl.DataFrame, pl.DataFrame]:
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

    Args (parsed from the command line):
    metadata - a pl.dataframe that contains the whole-dataset metadata
    dist_scores - a pl.dataframe that contains distance scores for each
                cluster centroid
    cluster_meta - a pl.dataframe that contains metadata from the
                clustering algorithm.
    """

    # unpack required inputs
    metadata = cluster_tables.full_meta
    dist_scores = cluster_tables.dist_scores
    cluster_meta = cluster_tables.cluster_meta

    # join distance scores onto the cluster metadata, which will add
    # scores to the centroids from each cluster
    assert "Accession" in cluster_meta.columns
    cluster_meta = cluster_meta.join(dist_scores, on="Accession", how="left")

    # Expand cluster centroid distance to all cluster members if preclustering was performed
    if cluster_tables.whether_preclust:
        # merge the month and cluster index columns to make a unique identifier
        cluster_meta = cluster_meta.with_columns(
            [pl.format("{}_{}", "Month", "Index").alias("Month Index")]
        )

        # Use a series of Polars expressions to:
        # 1 - replace sequence length entries in "H" entries with the cluster
        # size entries from the associated "C" entries.
        # 2 - Spread the distance scores from each centroid to all the hit
        # sequences in their clusters, and
        # 3 - select only the cluster size, distance, accession, and month
        # columns for joining with the big metadata
        assert "Type" in cluster_meta.columns
        assert "Month Index" in cluster_meta.columns
        assert "Size" in cluster_meta.columns
        assert "Distance Score" in cluster_meta.columns
        centroid_data = cluster_meta.filter(pl.col("Type") == "C").select(
            ["Month Index", "Size", "Distance Score"]
        )
        cluster_meta = (
            cluster_meta.join(
                centroid_data, on="Month Index", how="left", suffix="_right"
            )
            .select(["Accession", "Month Index", "Size_right", "Distance Score_right"])
            .rename(
                {"Size_right": "Cluster Size", "Distance Score_right": "Distance Score"}
            )
        )

        # filter out any clusters with only one entry (this shouldn't
        # be necessary when run in the context of the pipeline, but
        # we include it here just to be safe.)
        cluster_meta = (
            cluster_meta.group_by("Month Index")
            .agg(pl.col("Month Index").count().alias("count"))
            .join(cluster_meta, on="Month Index", how="left")
            .filter(pl.col("count") > 1)
            .select(cluster_meta.columns)
        )

    # join with the big metadata
    high_dist_meta = metadata.join(cluster_meta, on="Accession", how="left")

    # sort the new high distance metadata by distance score
    # in descending order and remove any null rows
    high_dist_meta = (
        high_dist_meta.sort(by="Distance Score", descending=True, nulls_last=True)
        .filter(pl.col("Distance Score").is_not_null())
        .unique()
    )

    # And finally, filter down to distances above the retention threshold while
    # also visualizing the distance score distribution
    retention_threshold = float(
        numpy.quantile(high_dist_meta["Distance Score"], (stringency / 1000))
    )
    try:
        await visualize_distance_scores(high_dist_meta, retention_threshold)
    except Exception as traceback:  # pylint: disable=W0718
        logger.debug("Plotting the distance score distribution failed: {}", traceback)
    logger.info("Distance score retention threshold: {}", retention_threshold)
    logger.info(
        "Max-observed distance score: {}",
        max(high_dist_meta.select("Distance Score").to_series().to_list()),
    )
    high_dist_meta = high_dist_meta.filter(
        pl.col("Distance Score") >= retention_threshold
    )

    # separate out the accessions
    accessions = high_dist_meta.select("Accession")

    return high_dist_meta, accessions


async def main():
    """
    Main function. Main ties together all the above functions if
    they are being invoked as a script. It also works to "glue"
    in a subprocess command using Seqkit to handle FASTA filtering
    (Seqkit, written in Go, is much faster at filtering large
    numbers of FASTA records than Python).
    """

    # initialize logger
    logger.remove()
    logger.add(sys.stderr, backtrace=True, diagnose=True, colorize=True, enqueue=True)

    # retrieve file paths and settings from keyword command line arguments
    metadata_name, fasta_path, stringency, workingdir = parse_command_line_args()

    # check that the provided input files exist
    assert os.path.isfile(
        metadata_name
    ), "Provided metadata path does not point to a file that exists."
    assert os.path.isfile(
        fasta_path
    ), "Provided FASTA path does not point to a file that exists."
    assert os.path.isdir(workingdir), "Provided working directory path does not exist."

    # Retrieve a quantile to use to define a retention threshold downstream
    strict_quant = await quantify_stringency(stringency)
    logger.info("Selected quantile stringency: {}", strict_quant)

    # List all files in current directory that end with '-dist-matrix.csv'
    logger.info("Collecting available distance matrix files.")
    files = [
        file for file in os.listdir(workingdir) if file.endswith("-dist-matrix.csv")
    ]

    # Remove '-dist-matrix.csv' from each file name
    yearmonths = [file.replace("-dist-matrix.csv", "") for file in files]

    # retrieve all candidate metadata
    logger.info("Reading candidate metadata files.")
    cluster_tables = await read_metadata_files(metadata_name, yearmonths, workingdir)

    # pile up all metadata from clustering, distance-calling, and the original database
    logger.info(
        "Piling up all metadata from clustering, distance-calling, and the original dataset"
    )
    high_dist_meta, accessions = await collate_metadata(
        cluster_tables,
        strict_quant,
    )

    # Use an assertion to spare the computer some effort if no candidates were found
    assert accessions.shape[0] > 0, "No candidate sequences identified."

    # write the metadata and the accessions
    logger.info("Writing high-distance metadata and accession.")
    high_dist_meta.write_csv(
        f"{workingdir}/high_distance_candidates.tsv", separator="\t"
    )
    accessions.write_csv(
        f"{workingdir}/high-dist-accessions.txt", include_header=False, separator="\t"
    )

    # Use Seqkit to separate out high distance sequences based on the metadata
    logger.info("Filtering FASTA records to high distance accessions only.")
    cmd = [
        "seqkit",
        "grep",
        "-f",
        "high-dist-accessions.txt",
        fasta_path,
        "-o",
        f"{workingdir}/high_distance_candidates.fasta",
    ]
    subprocess.run(cmd, check=True)

    # remove the accessions file that is now unnecessary
    os.remove("high-dist-accessions.txt")


if __name__ == "__main__":
    asyncio.run(main())
