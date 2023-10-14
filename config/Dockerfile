FROM ubuntu:20.04

# Set working directory
WORKDIR /scratch

# Set environment variables
ENV DEBIAN_FRONTEND noninteractive
ENV TZ America/New_York

# Install dependencies
COPY ./apt-get-packages.txt /tmp/package_list.txt
RUN apt-get update && \
    awk '/^ii/ { printf("apt-get install -y %s=%s\n", $2, $3) }' /tmp/package_list.txt > /tmp/install_packages.sh && \
    chmod +x /tmp/install_packages.sh && \
    /tmp/install_packages.sh && \
    rm -rf /tmp && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install miniconda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh
ENV PATH=/opt/conda/bin:$PATH

# Install mamba
RUN conda install -y -n base -c conda-forge mamba=1.5.1

# Install all other dependencies, save for the two R packages below and the julia packages
COPY ./conda_env.yaml /tmp/conda_env.yaml
RUN mamba env update --file /tmp/conda_env.yaml && \
    mamba clean --all
ENV NXF_HOME=/scratch/.nextflow

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
COPY ALPINE.jl/ /root/.julia/environments/v1.9
RUN julia -e 'using Pkg; \
            Pkg.activate(joinpath(DEPOT_PATH[1], "environments", "v1.9")); \
            Pkg.instantiate(); \
            Pkg.build(); \
            Pkg.precompile(); \
            using TimeZones; \
            using ALPINE'

# Make sure nothing is amiss and that julia is still in the path
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
