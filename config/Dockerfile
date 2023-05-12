FROM ubuntu:20.04

# Set environment variables
ENV DEBIAN_FRONTEND noninteractive
ENV TZ America/New_York

# Install dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    build-essential \
    ca-certificates \
    build-essential \
    curl \
    wget \
    unzip \
    gzip \
    make \
    gnupg \
    libz-dev \
    libbz2-dev \
    liblzma-dev \
    libssl-dev \
    zlib1g-dev \
    libncurses5-dev \
    libgdbm-dev \
    libnss3-dev \
    libreadline-dev \
    libffi-dev \
    libsqlite3-dev \
    libbz2-dev \
    liblzma-dev \
    libxml2-dev \
    libxslt-dev \
    libcairo2-dev \
    libpango1.0-dev \
    libpng-dev \
    libxt-dev \
    libxaw7-dev \
    libglu1-mesa-dev \
    libglfw3-dev \
    xvfb \
    muscle \
    python3-pip \
    git \
    libssl-dev \
    libcurl4-openssl-dev \
    default-jre \
    default-jdk \
    libnss-sss \
    procps \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install miniconda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh
ENV PATH /opt/conda/bin:$PATH

# Install psutil
RUN conda install -c conda-forge psutil

# Install Biopython
RUN conda install -c conda-forge biopython

# Install pandas
RUN conda install -c anaconda pandas

# Install polars
RUN conda install -c conda-forge polars

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
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Install Julia
RUN wget https://julialang-s3.julialang.org/bin/linux/x64/1.9/julia-1.9.0-linux-x86_64.tar.gz && \
    tar -xzf julia-1.9.0-linux-x86_64.tar.gz -C /opt && \
    ln -s /opt/julia-1.9.0/bin/julia /usr/local/bin/julia && \
    rm julia-1.9.0-linux-x86_64.tar.gz

# Install Julia packages
RUN julia -e 'using Pkg; Pkg.add("BioSequences"); Pkg.add("CSV"); Pkg.add("DataFrames"); Pkg.add("BioAlignments"); Pkg.add("PopGen"); Pkg.add("FASTX"); Pkg.add("XAM")'
# ; Pkg.add("SubstitutionModels"); Pkg.add("KmerAnalysis"); Pkg.add("GFF3"); Pkg.add("GenomicFeatures"); Pkg.add("BioServices"); Pkg.add("BioFetch"); Pkg.add("BED"); Pkg.add("BioTools"); Pkg.add("NaturalSelection"); Pkg.add("Phylogenies")'

# Set working directory
WORKDIR /app
