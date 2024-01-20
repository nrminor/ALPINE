# ALPINE: Anachronistic Lineage and Persistent INfection Explorer
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.1-23aa62.svg)](https://www.nextflow.io/) [![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/) [![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/) [![Docker CI/CD](https://github.com/nrminor/ALPINE/actions/workflows/docker-image.yaml/badge.svg)](https://github.com/nrminor/ALPINE/actions/workflows/docker-image.yaml) [![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)

## Overview
ALPINE is a high-throughput NextFlow pipeline for mutation-agnostic discovery of evolutionarily advanced pathogens and prolonged infection candidates in public sequence databases. It was designed to probe large numbers of pathogen sequences, numbering in the millions to tens of millions, for _highly evolved_ and/or _anachronistic_ sequences. These two concepts are core to understanding the pipeline's outputs, as they encompass two ways of identifying pathogens that may be greater public health risks than other co-occuring pathogens.

### What the pipeline looks for
**Highly evolved or evolutionarily advanced pathogens** are the bread and butter of this pipeline. It considers a sequence highly evolved if it carries an exceptional number of substitutions compared to co-circulating sequences. When placed on a phylogenetic tree for the pathogen of interest, these sequences end up on unusually long branches. In the case of SARS-CoV-2, these long branches may be reminiscent of the original Omicron lineage, which was highly advanced compared the co-circulating lineages at the time.

This pipeline uses clustering and distance matrices (see Workflow Steps; manuscript in prep) to flag advanced candidate sequences along with their associated metadata. It also generates a variety of visualizations that can help conceptualize these candidates in mutational space. In addition, the pipeline examines the similarity of all highly evolved candidates to flag sequences that are particularly likely to have come from the same prolonged infection.

Here, we define an **anachronistic lineage** as a lineage that is detected at a time when that it is otherwise extinct, or at least no longer prevalently circulating. These lineages are likely to come from long-term, single-host infections such as SARS-CoV-2 prolonged infections, though they may also represent reverse spillbacks from animal reservoirs. In addition to flagging highly evolved candidates, the pipeline also flags these anachronistic sequences and their metadata.

When both highly evolved and anachronistic candidates are identified, the pipeline identifies any overlapping or **"double candidates"**. These candidates, which have undergone a great deal of genomic evolution _and_ are no longer a member of a circulating lineage, are the most likely to stem from prolonged infections. These prolonged infections have either lasted an exceptionally long time, or, for reasons that are under active research, have been characterized by an exceptionally high substitution rate.

Overall, the candidates that are flagged in each of the above approaches are meant to be useful in a variety of downstream analyses, including dN/dS inference, tree-building, amino acid mutation analysis, and more. They can also be used in real time by public health agencies and research labs to track the emergence of potentially high-risk pathogens. Our hope is that sequences flagged here can then spur more rapid epidemiological investigation in cases where an outbreak has been identified, and more active research in cases where a candidate is identified retroactively.

### How it works
The pipeline uses a variety of tools to run quickly, use a minimum of memory and disk space, and work the same way across many kinds of computers.
- Clustering with [_vsearch_](https://github.com/torognes/vsearch), an open-source alternative to USEARCH written in C++ by Torbjørn Rognes and colleagues. This tool makes it possible to reduce the dimensionality of immense sequence datasets.
- [_Nextflow_](https://www.nextflow.io/) as the workflow language and underlying task manager. Nextflow makes it possible to parallelize and scale tasks across any number of computer cores. It also allocates resources and runs containers to make sure each step of the pipeline has the software it needs to run.
- [_Docker_](https://www.docker.com/) (or Apptainer) as the container engine. These tools ensure that the pipeline runs the same way on any machine, regardless of whether its software requirements are installed.
- [_Apache Arrow_](https://arrow.apache.org/) for data tables. Using the Arrow spec instead of traditional CSVs makes it faster to read and write large tables of data with millions of rows. It also makes working with the data in those tables more memory-efficient. This is part of why the pipeline does not require an HPC cluster to run; in principle, there's no reason it couldn't run on a laptop.
- [_ZStandard_](https://facebook.github.io/zstd/) for sequence data compression. ZStandard is a fast and effective compression algorithm that is used on sequence data throughout the pipeline. This makes it tractable to probe millions of sequences without requiring terabytes of available disk space.
- [_Polars DataFrames_](https://pola.rs/), a blazingly fast dataframe frontend and query engine written in Rust. We use Polars to filter and transform data tables throughout each executation of ALPINE.
- Core command line interface written in [_Rust_](https://www.rust-lang.org/). Rust is a statically-compiled, lower-level language that has been gaining traction in the life sciences and in data science writ large. It comes with modern tooling, guaranteed memory safety, and the ability to build an API that feels high-level without sacrificing performance or granular memory management. ALPINE's utilities for calling distance matrices on large sequences, filtering large FASTAs based on the number of masked bases (N's), and sorting each FASTA entry into a separate file for each month are all written in Rust. And most Python scripts in the project `bin/` directory are similarly fast because of their use of Rust-based Polars DataFrames (or more specifically, [LazyFrame queries](https://docs.pola.rs/user-guide/lazy/using/)).

For more details on exactly what the pipeline does, see the Workflow Overview section below.

## Quick Start

If Docker and NextFlow are already installed on your system, simply run the following command in the directory of your choice to execute the workflow with default settings (more below):

```
nextflow run nrminor/ALPINE
```

## Detailed Setup Instructions

## Configuration

## Workflow Overview

## Acknowledgements
