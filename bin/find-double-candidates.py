#!/usr/bin/env python3

"""
`find-double-candidates.py` does just that: find the Venn overlap between any
provided anachronistic data and any provided high-distance data, i.e., double
candidates. It filters both metadata and sequence data so that users may focus
on the highest-interest accessions in the provided dataset.
"""

import asyncio
import glob
import inspect
import os
import sys
from typing import List, Set, Tuple

import polars as pl
from Bio import SeqIO
from loguru import logger
from pydantic import validate_call


@logger.catch(reraise=True)
@validate_call
async def list_metadata_files() -> Tuple[str, ...]:
    """
    List the available TSV files in the current working directory. This
    assumes that the only available TSV files are ALPINE metadata outputs.
    As such, it is potentially error-prone to use this script outside the
    context of the ALPINE Nextflow wrapper.
    """
    return tuple(glob.glob("*.tsv"))


@logger.catch(reraise=True)
@validate_call
async def list_fasta_files() -> Tuple[str, ...]:
    """
    List the available FASTA files in the current working directory. This
    assumes that the only available TSV files are ALPINE metadata outputs.
    As such, it is potentially error-prone to use this script outside the
    context of the ALPINE Nextflow wrapper.
    """
    return tuple(glob.glob("*.fasta"))


@logger.catch(reraise=True)
@validate_call
async def sort_meta_files(meta_files: Tuple[str, ...]) -> List[str]:
    """
    Sort the metadata files by size and keep the two largest. This provides
    maximal grist for Venn-overlapping potential double-candidates.
    """

    assert len(meta_files) > 1, "Too few metadata files provided for comparison."

    return sorted(meta_files, key=os.path.getsize, reverse=True)[0:2]


@logger.catch(reraise=True)
@validate_call
async def sort_fasta_files(fasta_files: Tuple[str, ...]) -> List[str]:
    """
    Sort the sequence files by size and keep the two largest. This provides
    maximal grist for Venn-overlapping potential double-candidates.
    """

    assert len(fasta_files) > 0, "Too few FASTA files provided for comparison."

    if len(fasta_files) > 1:
        return sorted(fasta_files, key=os.path.getsize, reverse=True)[0:2]

    return sorted(fasta_files, key=os.path.getsize, reverse=True)


@logger.catch(reraise=True)
async def read_metadata(sorted_meta: List[str]) -> Tuple[pl.LazyFrame, pl.LazyFrame]:
    """
    Asynchronously read metadata files.
    """

    left_meta: pl.LazyFrame = pl.scan_csv(sorted_meta[0], separator="\t")
    right_meta: pl.LazyFrame = pl.scan_csv(sorted_meta[1], separator="\t")

    return left_meta, right_meta


@logger.catch(reraise=True)
async def overlap_metadata(
    left_meta: pl.LazyFrame, right_meta: pl.LazyFrame
) -> pl.LazyFrame:
    """
    Inner-join the provided metadata files by the "Accession" column while
    remaining sensitive to possible anachronicity or infection duration
    columns.
    """

    # run some assertions to make sure we don't get into trouble here
    assert (
        left_meta.width > 0
    ), "No columns found in the larger metadata. They may either be empty or not tab-delimited."
    assert (
        right_meta.width > 0
    ), "No columns found in the smaller metadata. They may either be empty or not tab-delimited."
    left_cols = left_meta.collect().columns
    right_cols = right_meta.collect().columns
    assert "Accession" in left_cols, "Accession column not found in larger metadata."
    assert "Accession" in right_cols, "Accession column not found in smaller metadata."

    # initialize double_candidates so scoping is clear
    double_candidates: pl.LazyFrame
    if "Anachronicity (days)" in right_meta:
        double_candidates = left_meta.join(
            right_meta.select("Accession", "Anachronicity (days)"),
            on="Accession",
            how="inner",
        )

    if "Anachronicity (days)" in left_meta:
        double_candidates = right_meta.join(
            left_meta.select("Accession", "Anachronicity (days)"),
            on="Accession",
            how="inner",
        )

    if "infection_duration" in right_meta:
        double_candidates = left_meta.join(
            right_meta.select("Accession", "infection_duration"),
            on="Accession",
            how="inner",
        )

    if "infection_duration" in left_meta:
        double_candidates = right_meta.join(
            left_meta.select("Accession", "infection_duration"),
            on="Accession",
            how="inner",
        )

    if (
        "Anachronicity (days)" not in left_meta
        and "infection_duration" not in left_meta
    ) and (
        "Anachronicity (days)" not in right_meta
        and "infection_duration" not in right_meta
    ):
        logger.info("No columns to allow reconciling two datasets found.")
        return double_candidates

    return double_candidates


@logger.catch(reraise=True)
async def get_common_accessions(double_candidates: pl.LazyFrame) -> Set[str]:
    """
    Execute a Polars query to extract a unique set of double candidate
    accessions. These will be used to filter the provided FASTA records
    down to those that are double candidates.
    """

    return set(double_candidates.select("Accession").collect().to_series().to_list())


@logger.opt(lazy=True).catch(reraise=True)
@validate_call
async def filter_fasta(candidate_set: Set[str], fasta_path: str) -> None:
    """
    Filter FASTA records so that the output FASTA only contains accessions
    that are double candidates.
    """

    # Make sure the FASTA file still exists
    assert os.path.isfile(fasta_path), "FASTA name provided not found."

    # fasta record ticker
    ticker: int = 0

    with open(fasta_path, "r", encoding="utf-8") as input_handle, open(
        "double_candidates.fasta", "w", encoding="utf-8"
    ) as output_handle:
        for record in SeqIO.parse(input_handle, "fasta"):
            # Get the first part of the sequence identifier
            accession = record.id.split(" ")[0]

            # ignore non-double-candidates
            if accession not in candidate_set:
                continue

            # write accessions in the set
            ticker += 1
            SeqIO.write(record, output_handle, "fasta")

    assert (
        len(candidate_set) == ticker
    ), "Could not find all accessions in the double candidate set in the largest provided FASTA."


async def main() -> None:
    """
    Main manages the asynchronous runtime that flows data through
    the above defined functions.
    """

    logger.remove()
    logger.add(sys.stderr, backtrace=True, diagnose=True, colorize=True)

    # retrieve available metadata filenames
    logger.info("Listing all available metadata and FASTA files.")
    meta_files = await list_metadata_files()
    logger.info("Found the following metadata files: {}", meta_files)

    # retrieve available FASTA filenames
    fasta_files = await list_fasta_files()
    logger.info("Found the following FASTA file(s): {}", fasta_files)

    assert len(meta_files) >= 2, "Not enough metadata files for comparison."

    # find large files that will be joined to
    logger.info("Determining two largest FASTA and metadata files to inner-join.")
    sorted_meta = await sort_meta_files(meta_files)
    sorted_fastas = await sort_fasta_files(fasta_files)

    assert (
        len(sorted_meta) == 2
    ), "Failed to find at least two metadata files to compare."
    assert (
        len(sorted_fastas) >= 1
    ), "Failed to find at least one FASTA files to compare."

    # read metadata
    logger.info("Reading metadata files.")
    left_meta, right_meta = await read_metadata(sorted_meta)

    # venn-overlap metadata
    logger.info("Venn overlapping the accessions in the provided metadata files.")
    double_candidates = await overlap_metadata(left_meta, right_meta)

    # exit if no data was found
    if not isinstance(double_candidates, pl.LazyFrame):
        sys.exit(
            inspect.cleandoc(
                """
                The script was unable to find columns 'Anachronicity (days)' or 'infection_duration'
                in the two largest provided metadata
                """
            )
        )

    # get accessions to filter out of the FASTA
    common_accessions = await get_common_accessions(double_candidates)

    # filter the FASTA as the last step
    logger.info("Filtering input FASTA records to double candidates.")
    await filter_fasta(common_accessions, sorted_fastas[0])


if __name__ == "__main__":
    asyncio.run(main())
