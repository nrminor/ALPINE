#!/usr/bin/env python3

"""
REPORT REPEAT LINEAGE SEQUENCES AND METADATA
--------------------------------------------

Accessions flagged in this script clustered together in the
so-called "meta-clustering", where high-distance sequences
were run through vsearch --cluster_fast with an even tighter
identity threshold, 99.99%. Sequences that cluster together
are likely to represent the same, high-distance viral lineage
that is persistently infecting one host, though other scenarios
are possible.
"""

import argparse
import asyncio
import os
import sys
from functools import lru_cache
from pathlib import Path
from typing import Tuple

import Bio.SeqIO
import pandas as pd  # type: ignore
from loguru import logger
from pydantic import validate_call


@validate_call
def parse_command_line_args() -> Tuple[Path, Path]:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "cluster_table_path", help="The path to the cluster table file."
    )
    parser.add_argument(
        "metadata_path", type=Path, help="The path to the metadata file."
    )
    args = parser.parse_args()
    return args.cluster_table_path, args.metadata_path


@logger.catch(reraise=True)
async def read_cluster_table(cluster_table_path: Path) -> pd.DataFrame:
    """Read the VSEARCH cluster table from a 'meta-cluster' on ALPINE candidates."""
    return pd.read_csv(cluster_table_path, sep="\t", header=None, na_values=["", "NA"])


@logger.catch(reraise=True)
async def get_unique_clusters(cluster_df: pd.DataFrame) -> pd.DataFrame:
    """Filter the cluster table to the members of all unique clusters"""
    unique_df = (
        cluster_df.loc[cluster_df[0] != "C"]
        .groupby(1)
        .filter(lambda x: len(x) > 1)
        .sort_values(1)
    )

    # Check that each accession is only present in one cluster
    assert unique_df[8].nunique() == len(unique_df)
    logger.info("All clusters confirmed to be unique.")

    return unique_df


@logger.catch(reraise=True)
@lru_cache
async def fasta_input(cluster: str) -> list:
    """Read the cluster FASTA and convert records to a list"""

    logger.info("Assessing repeat lineages in cluster {cluster}")
    return list(Bio.SeqIO.parse(f"meta-cluster-seqs{cluster}", "fasta"))


@logger.catch(reraise=True)
async def fasta_output(fasta_seq: list, cluster: str) -> None:
    """Write the cluster FASTA"""

    with open(f"repeat-lineage-{cluster}.fasta", "w", encoding="utf-8") as outfile:
        Bio.SeqIO.write(fasta_seq, outfile, "fasta")


@logger.catch(reraise=True)
async def read_metadata(metadata_path: Path) -> pd.DataFrame:
    """Read the full-dataset metadata for candidates"""

    return pd.read_csv(metadata_path, delimiter="\t", na_values=["", "NA"])


@logger.catch(reraise=True)
async def merge_metadata(
    unique_df: pd.DataFrame, metadata: pd.DataFrame
) -> pd.DataFrame:
    """Join repeat-lineage metadata onto the full-dataset metadata"""

    repeat_clusts = pd.DataFrame(
        {"Accession": unique_df[8], "Repeat Cluster Number": unique_df[1]}
    )
    merged_metadata = pd.merge(metadata, repeat_clusts, on="Accession", how="inner")

    return (
        merged_metadata.loc[~merged_metadata["Repeat Cluster Number"].isna()]
        .sort_values("Repeat Cluster Number")
        .drop(["Distance Score", "Cluster Size"], axis=1)
    )


async def main():
    """
    This function reads in a cluster table and metadata file,
    and then compiles FASTA sequences for each repeat cluster and
    collates metadata about these new clusters.

    Args (parsed from the command line):
    cluster_table_path: The path to the cluster table file.
    metadata_path: The path to the metadata file.

    Returns:
    None
    """

    # initialize logger settings
    logger.remove()
    logger.add(sys.stderr, colorize=True, backtrace=True, diagnose=True)

    # parse command line arguments
    cluster_table_path, metadata_path = parse_command_line_args()

    # check that the provided input files exist
    assert os.path.isfile(
        cluster_table_path
    ), "Provided cluster table path does not point to a file that exists."
    assert os.path.isfile(
        metadata_path
    ), "Provided metadata path does not point to a file that exists."

    # read the metadata
    logger.info("Processing VSEARCH clustering metadata.")
    cluster_df = await read_cluster_table(cluster_table_path)

    # make sure only unique clusters are represented
    unique_df = await get_unique_clusters(cluster_df)

    # Compile FASTA sequences for each cluster
    for cluster in set(unique_df[1].tolist()):
        logger.info("Assessing repeat lineages in cluster {}", cluster)
        fasta_seq = await fasta_input(cluster)
        await fasta_output(fasta_seq, cluster)

    # Collate metadata about these new clusters
    metadata = await read_metadata(metadata_path)

    # Create new table to store accessions and repeat cluster number
    logger.info("Merging repeat-lineage cluster metadata with dataset metadata.")
    merged_metadata = await merge_metadata(unique_df, metadata)

    # Write repeat lineage metadata
    logger.info("Writing repeat lineage metadata.")
    merged_metadata.to_csv(
        "repeat-lineage-metadata.tsv", na_rep="", index=False, sep="\t"
    )


if __name__ == "__main__":
    asyncio.run(main())
