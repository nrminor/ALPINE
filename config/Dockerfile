FROM ubuntu:20.04

# Set environment variables
ENV DEBIAN_FRONTEND noninteractive
ENV TZ America/New_York

# Install dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    build-essential \
    ca-certificates \
    curl \
    wget \
    unzip \
    gzip \
    make \
    zstd \
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
ENV PATH=/opt/conda/bin:$PATH

# Install mamba
RUN conda install -y -n base -c conda-forge mamba libarchive

# Install R channel packages
RUN mamba install -y -c conda-forge -c r -c bioconda r-base \
    r-tidyverse \
    r-ape \
    r-plotly \
    r-fst \
    r-devtools \
    r-units \
    r-sf \
    r-arrow \
    r-duckdb \
    bioconductor-biostrings

# Install outbreak.info and cdc R packages
RUN Rscript -e "install.packages('cdcfluview', repos = 'https://cloud.r-project.org/', lib='/opt/conda/lib/R/library')"
RUN Rscript -e "devtools::install_github('outbreak-info/R-outbreak-info', lib='/opt/conda/lib/R/library')"

# Install Rust and Apache Arrow  CLI tools
RUN curl https://sh.rustup.rs -sSf | sh -s -- -y && \
    cargo install csv2arrow && \
    cargo install csv2parquet

# Install Julia
RUN wget https://julialang-s3.julialang.org/bin/linux/x64/1.9/julia-1.9.0-linux-x86_64.tar.gz && \
    tar -xzf julia-1.9.0-linux-x86_64.tar.gz -C /opt && \
    ln -s /opt/julia-1.9.0/bin/julia /usr/local/bin/julia && \
    rm julia-1.9.0-linux-x86_64.tar.gz

# Copy Julia functions to precompile as a module
ENV JULIA_DEPOT_PATH=/root/.julia
COPY LongInfectionFinder.jl/ /root/.julia/environments/v1.9
RUN julia -e 'using Pkg; \
            Pkg.activate(joinpath(DEPOT_PATH[1], "environments", "v1.9")); \
            Pkg.instantiate(); \
            Pkg.precompile()'

# Make sure nothing is amiss and that julia is still in the path
ENV PATH=$PATH:/opt/julia-1.9.0/bin:/root/.julia/compiled/v1.9

# Install python libraries and other packages
RUN mamba install -y -c bioconda -c anaconda -c conda-forge \
    psutil \
    libzip \
    biopython pandas polars numpy numba \
    vsearch \
    ncbi-datasets-cli=14.2.2 \
    seqkit \
    pangolin \
    usher \
    nextstrain \
    python-duckdb

# Install nextflow
RUN wget -qO- https://get.nextflow.io | bash && \
    mv nextflow /usr/local/bin/

# make sure bin files are executable
RUN chmod +x /usr/local/bin/* && \
    chmod +x /opt/julia-1.9.0/bin/* && \
    chmod +x /opt/conda/bin/* && \
    chmod +rw /root/ && \
    chmod -R +rwx /root/ && \
    chmod -R +rwx /root/.julia/ && \
    chmod -R +rwx /root/.julia/compiled/ && \
    chmod -R +rwx /root/.julia/compiled/v1.9/*

# Set working directory
WORKDIR /scratch
