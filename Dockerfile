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

# Install all other dependencies, save for the two R packages below and the julia packages
COPY ./config/conda_env.yaml /tmp/conda_env.yaml
RUN mamba env update --file /tmp/conda_env.yaml && \
    mamba install -y -c conda-forge pydantic && \
    mamba clean --all && \
    rm /tmp/conda_env.yaml
ENV NXF_HOME=/scratch/.nextflow

# install Rust
RUN apt-get update && \
    apt-get install -y curl && \
    apt-get clean
RUN mkdir -m777 /opt/rust /opt/.cargo
ENV RUSTUP_HOME=/opt/rust CARGO_HOME=/opt/.cargo PATH=/opt/.cargo/bin:$PATH
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | bash -s -- -y && \
    bash "/opt/.cargo/env"

# Install Rust tools
RUN cd /opt && \
    git clone https://github.com/nrminor/ALPINE-core.git && \
    cd ALPINE-core && \
    cargo install --path . \
    --root /opt/.cargo
RUN cargo install prqlc \
    --root /opt/.cargo
RUN apt-get update && \
    apt-get install -y clang && \
    apt-get clean && \
    cargo install qsv \
    --locked --root /opt/.cargo \
    --features=apply,foreach,polars,to,to_parquet,self_update,feature_capable
# RUN cargo install nu --features=dataframe --root /opt/.cargo

# Install outbreak.info and cdc R packages
RUN Rscript -e "install.packages('cdcfluview', repos = 'https://cloud.r-project.org/', lib='/opt/conda/lib/R/library', clean = TRUE)"
RUN Rscript -e "devtools::install_github('outbreak-info/R-outbreak-info', lib='/opt/conda/lib/R/library', clean = TRUE)"

# Install Julia
RUN wget https://julialang-s3.julialang.org/bin/linux/x64/1.9/julia-1.9.0-linux-x86_64.tar.gz && \
    tar -xzf julia-1.9.0-linux-x86_64.tar.gz -C /opt && \
    ln -s /opt/julia-1.9.0/bin/julia /usr/local/bin/julia && \
    rm julia-1.9.0-linux-x86_64.tar.gz

# Copy Julia functions to precompile as a module
ENV JULIA_DEPOT_PATH=/root/.julia
ENV JULIA_SCRATCH_TRACK_ACCESS=0
COPY ./config/ALPINE.jl/ /root/.julia/environments/v1.9
RUN julia -e 'using Pkg; \
            Pkg.activate(joinpath(DEPOT_PATH[1], "environments", "v1.9")); \
            Pkg.instantiate(); \
            using PackageCompiler; \
            create_sysimage(:ALPINE, sysimage_path="/opt/.julia/alpine.so")'

# Install go and still
RUN apt-get update && \
    apt-get install -y golang-go && \
    apt-get clean
RUN cd /opt && \
    git clone https://github.com/danielecook/still.git && \
    cd still && \
    go build
ENV PATH=$PATH:/opt/still

# Make sure nothing is amiss and that everything is still in the path
ENV NXF_HOME=/scratch/.nextflow
ENV PATH=$PATH:/opt/julia-1.9.0/bin:/scratch/.julia/compiled/v1.9

# make sure bin files are executable
RUN chmod +x /usr/local/bin/* && \
    chmod +x /opt/julia-1.9.0/bin/* && \
    rm -f /opt/conda/bin/cpp && \
    chmod +x /opt/conda/bin/* && \
    chmod +rw /root/ && \
    chmod -R +rwx /root/ && \
    chmod -R +rwx /root/.julia/ && \
    chmod -R +rw /root/.julia/logs/ && \
    chmod -R +rwx /root/.julia/logs/* && \
    chmod -R +rwx /root/.julia/compiled/ && \
    chmod -R +rwx /root/.julia/compiled/v1.9/* && \
    chmod -R +rwx /root/.julia/logs/* && \
    chmod -R +rwx /root/.julia/packages/ && \
    chmod -R 755 /root/.julia/scratchspaces/

# make sure shells are bash
CMD ["/bin/bash"]
