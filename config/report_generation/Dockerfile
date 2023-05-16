# Use an official Miniconda runtime as a parent image
FROM continuumio/miniconda3:23.3.1-0

# Set the working directory to /app
WORKDIR /app

# Update the system and install 'ps' utility
RUN apt-get update && apt-get install -y procps

# Install mamba
RUN conda install -y -c conda-forge mamba

# Install R channel packages
RUN mamba install -y -c conda-forge -c r -c bioconda r-base r-tidyverse r-ape r-adephylo 

# Activate the environment
ENV PATH /opt/conda/bin:$PATH
