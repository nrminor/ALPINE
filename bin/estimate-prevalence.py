#!/usr/bin/env python3

"""
`estimate-prevalence.py` simply computes and prints the prevalence of ALPINE
'double candidates'.

```
usage: estimate-prevalence.py [-h] [--early_stats EARLY_STATS] [--late_stats LATE_STATS]

options:
  -h, --help            show this help message and exit
  --early_stats EARLY_STATS, --e EARLY_STATS
                        A seqkit stats table from after filtering but before clustering.
  --late_stats LATE_STATS, --l LATE_STATS
                        A seqkit stats table from after clustering.
```
"""

import argparse
import os
import sys
from functools import lru_cache
from pprint import pprint
from typing import Tuple

import numpy as np
import polars as pl
from loguru import logger
from pydantic import validate_call


@validate_call
def parse_command_line_args() -> Tuple[str, str]:
    """parse command line arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--early_stats",
        "--e",
        required=False,
        default="early_stats.tsv",
        help="A seqkit stats table from after filtering but before clustering.",
    )
    parser.add_argument(
        "--late_stats",
        "--l",
        required=False,
        default="late_stats.tsv",
        help="A seqkit stats table from after clustering.",
    )
    args = parser.parse_args()
    return args.early_stats, args.late_stats


def _pipe(data, *funcs):
    for func in funcs:
        data = func(data)
    return data


@lru_cache
def _calculate_prevalence(data):
    candidate_count, sample_size = data
    return (candidate_count / sample_size) * 100


@lru_cache
def _round_prevalence(prevalence):
    return np.round(prevalence, 5)


@logger.opt(lazy=True).catch(reraise=True)
def estimate_prevalence(
    early_stats: pl.LazyFrame, late_stats: pl.LazyFrame
) -> Tuple[float, int]:
    """
    Function `estimate_prevalence` uses the outputs from `seqkit stats` to
    estimate the proportion of sequences flagged at the end of ALPINE in the
    sequences used as the workflow's input.
    """

    assert (
        "num_seqs" in early_stats.columns
    ), "Required column 'num_seqs' not found in provided seqkit early stats table."
    assert (
        "num_seqs" in late_stats.columns
    ), "Required column 'num_seqs' not found in provided seqkit late stats table."

    # make sure that candidates were found at all before pulling out
    # the final sample size
    candidate_count: int = (
        0
        if late_stats.fetch(1).shape[0]
        else late_stats.select(pl.col("num_seqs")).collect().to_series().to_list()[0]
    )

    # pull out the input sample size
    sample_size: int = (
        early_stats.select(pl.col("num_seqs")).collect().to_series().to_list()[0]
    )

    # compute the prevalence estimate
    prevalence = _pipe(
        (candidate_count, sample_size), _calculate_prevalence, _round_prevalence
    )

    logger.info(
        "{candidate_count} candidates found among {sample_size} input sequences for a prevalence of {prevalence}%."
    )

    return (prevalence, sample_size)


def main() -> None:
    """
    Main controls dataflow.
    """

    logger.remove()
    logger.add(sys.stderr, backtrace=True, diagnose=True, colorize=True)

    early_file, late_file = parse_command_line_args()

    if not os.path.isfile(early_file):
        sys.exit(f"Early stats file, {early_file}, does not exist.")
    if not os.path.isfile(early_file):
        sys.exit(f"Late stats file, {late_file}, does not exist.")

    early_stats = pl.scan_csv(early_file, separator="\t")
    late_stats = pl.scan_csv(late_file, separator="\t")

    prevalence, sample_size = estimate_prevalence(early_stats, late_stats)

    pprint(
        f"{prevalence}% of ${sample_size} sequences were flagged as highly evolved and anachronistic."
    )


if __name__ == "__main__":
    main()
