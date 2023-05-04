FROM ubuntu:20.04

# Set environment variables
ENV DEBIAN_FRONTEND noninteractive
ENV TZ America/New_York

RUN apt-get update && \
    apt-get install -y \
    wget \
    python3 \
    python3-pip \
    unzip \
    libz-dev \
    libbz2-dev \
    liblzma-dev \
    build-essential \
    curl \
    git \
    libssl-dev \
    libcurl4-openssl-dev \
    default-jre \
    default-jdk \
    libnss-sss \
    procps

# Install miniconda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh
ENV PATH /opt/conda/bin:$PATH

# Install psutil
RUN pip3 install psutil

# Install vsearch
RUN conda install -c bioconda vsearch