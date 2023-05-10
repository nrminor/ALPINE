FROM ubuntu:20.04

# Install required packages
RUN apt-get update && \
    apt-get install -y \
    build-essential \
    wget \
    curl \
    gcc \
    make \
    zlib1g-dev \
    libz-dev \
    libbz2-dev \
    liblzma-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    muscle \
    python3-pip \
    libssl-dev \
    libcurl4-openssl-dev \
    default-jre \
    default-jdk \
    libnss-sss \
    procps

# Install psutil
RUN pip3 install psutil

# Install Biopython
RUN pip3 install biopython

# Install Biopython
RUN pip3 install pysam

# Install pandas
RUN pip3 install pandas

# Install polars
RUN pip3 install polars

# Install miniconda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh
ENV PATH /opt/conda/bin:$PATH

# Install vsearch
RUN conda install -c bioconda vsearch

# Install seqtk
RUN wget https://github.com/lh3/seqtk/archive/refs/tags/v1.3.tar.gz && \
    tar -zxvf v1.3.tar.gz && \
    cd seqtk-1.3 && \
    make && \
    mv seqtk /usr/local/bin/ && \
    cd .. && \
    rm -rf seqtk-1.3 v1.3.tar.gz && \
    apt-get remove -y wget gcc make zlib1g-dev && \
    apt-get autoremove -y && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Set working directory
WORKDIR /app