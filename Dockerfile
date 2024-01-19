FROM ubuntu:20.04

# Set working directory
WORKDIR /scratch

# Set environment variables
ENV DEBIAN_FRONTEND noninteractive
ENV TZ America/New_York

# Install dependencies
# Use `dpkg -l > apt-get.lock` inside a container to print package versions after a successful build
COPY ./config/apt-get.lock /tmp/apt-get.lock
RUN apt-get update && \
    awk '/^ii/ { printf("apt-get install -y %s=%s\n", $2, $3) }' /tmp/apt-get.lock > /tmp/install_packages.sh && \
    chmod +x /tmp/install_packages.sh && \
    /tmp/install_packages.sh && \
    rm /tmp/install_packages.sh && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install miniconda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py311_23.5.2-0-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh
ENV PATH=/opt/conda/bin:$PATH

# Install mamba
RUN conda config --add channels conda-forge && \
    conda update -n base --all && \
    conda install -y -n base \
    mamba=1.2.0=py311h51e3938_0 \
    libmamba=1.2.0=habaa8ee_0 \
    libmambapy=1.2.0=py311h2b46443_0 \
    conda-libmamba-solver=23.1.0=pyhd8ed1ab_0

# Install Python dependencies
# conda env export > environment.yml
# conda list -e > requirements.txt
COPY ./config/environment.yml /tmp/environment.yml
RUN mamba env update --file /tmp/environment.yml && \
    mamba clean --all && \
    rm /tmp/environment.yml
# COPY ./config/requirements.txt /tmp/requirements.txt
# RUN python3 -m pip install -r /tmp/requirements.txt && \
#     rm /tmp/requirements.txt
ENV NXF_HOME=/scratch/.nextflow

# install Rust
RUN mkdir -m777 /opt/rust /opt/.cargo
ENV RUSTUP_HOME=/opt/rust CARGO_HOME=/opt/.cargo PATH=/opt/.cargo/bin:$PATH
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | bash -s -- -y && \
    bash "/opt/.cargo/env"

# Install Rust tools
RUN cd /opt && \
    git clone https://github.com/nrminor/ALPINE-core.git && \
    cd ALPINE-core && \
    cargo install --path . \
    --root /opt/.cargo && \
    cargo clean
RUN cargo install qsv \
    --locked --root /opt/.cargo \
    --features=apply,foreach,polars,to,to_parquet,self_update,feature_capable
# RUN cargo install nu --features=dataframe --root /opt/.cargo

# Install outbreak.info and cdc R packages
RUN Rscript -e "install.packages('cdcfluview', repos = 'https://cloud.r-project.org/', lib='/opt/conda/lib/R/library', clean = TRUE)"
RUN Rscript -e "devtools::install_github('outbreak-info/R-outbreak-info', lib='/opt/conda/lib/R/library', clean = TRUE)"

# Install the still spreadsheet unit-tester
RUN cd /opt && \
    git clone https://github.com/danielecook/still.git && \
    cd still && \
    go build
ENV PATH=$PATH:/opt/still

# Make sure nothing is amiss and that everything is still in the path
ENV NXF_HOME=/scratch/.nextflow

# make sure shells are bash
CMD ["/bin/bash"]
