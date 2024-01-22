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
from pathlib import Path
from pprint import pprint
from typing import Tuple

import numpy as np
import polars as pl
from loguru import logger
from pydantic import validate_call


@validate_call
def parse_command_line_args() -> Tuple[Path, Path]:
    """parse command line arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--early_stats",
        "--e",
        type=Path,
        required=True,
        default="early_stats.tsv",
        help="A seqkit stats table from after filtering but before clustering.",
    )
    parser.add_argument(
        "--late_stats",
        "--l",
        type=Path,
        required=True,
        default="late_stats.tsv",
        help="A seqkit stats table from after clustering.",
    )
    args = parser.parse_args()
    return args.early_stats, args.late_stats


def _pipe(data, *funcs):
    for func in funcs:
        data = func(data)
    return data


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
        "{} candidates found among {} input sequences ({}%).",
        candidate_count,
        sample_size,
        prevalence,
    )

    return (prevalence, sample_size)


def main() -> None:
    """
    Main controls dataflow.
    """

    # initialize new colorized logger after removing the default
    logger.remove()
    logger.add(sys.stderr, backtrace=True, diagnose=True, colorize=True)

    # parse out early and late files
    early_file, late_file = parse_command_line_args()

    # check that the provided files actually exist
    assert os.path.isfile(
        early_file
    ), f"Early stats file, {early_file}, does not exist."
    assert os.path.isfile(early_file), f"Late stats file, {late_file}, does not exist."

    # scan the stats files for lazy evaluation
    early_stats = pl.scan_csv(early_file, separator="\t")
    late_stats = pl.scan_csv(late_file, separator="\t")

    # compute prevalence percentage and sample size
    prevalence, sample_size = estimate_prevalence(early_stats, late_stats)

    # pretty-print for forwarding up to the workflow level
    pprint(
        f"{prevalence}% of ${sample_size} sequences were flagged as highly evolved & anachronistic."
    )


if __name__ == "__main__":
    main()
