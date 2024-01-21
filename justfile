just default:
    just --list

install-conda:
    wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py311_23.5.2-0-Linux-x86_64.sh -O ./miniconda.sh
    -/bin/bash ~/miniconda.sh -b
    rm miniconda.sh
alias conda := install-conda

install-mamba:
    conda config --add channels conda-forge
    conda update -n base --all
    conda install -y -n base \
    mamba=1.2.0=py311h51e3938_0 \
    libmamba=1.2.0=habaa8ee_0 \
    libmambapy=1.2.0=py311h2b46443_0 \
    conda-libmamba-solver=23.1.0=pyhd8ed1ab_0
alias mamba := install-mamba

python-tools:
    -just conda
    -just mamba
    mamba env update --file config/environment.yml
    mamba clean --all
alias py := python-tools

rust-tools:
    cargo install alpine
    cargo install nu --features=dataframe
    cargo install --locked qsv --features=apply,foreach,polars,to,to_parquet,self_update,feature_capable
alias rs := rust-tools

r-environment:
    #!/usr/bin/env Rscript
    install.packages('cdcfluview', repos = 'https://cloud.r-project.org/', lib='/opt/conda/lib/R/library', clean = TRUE)
    devtools::install_github('outbreak-info/R-outbreak-info', lib='/opt/conda/lib/R/library', clean = TRUE)
alias r := r-environment

go-tools:
    (
        -mkdir ~/bioinformatics && cd ~/bioinformatics && git clone https://github.com/danielecook/still.git
        cd ~/bioinformatics/still && go build
    )
    -echo 'export PATH=$PATH:~/bioinformatics/still' >> ~/.zprofile
    -echo 'export PATH=$PATH:~/bioinformatics/still' >> ~/.bash_profile	
alias go := go-tools

all-macos:
    -/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
    brew install gzip
    brew install unzip
    brew install gcc
    brew install pre-commit
    brew install go
    brew install r
    brew install curl
    brew install wget
    brew install llvm
    brew install make
    brew install cmake
    brew install zstd
    brew install seqkit
    brew install docker --cask
    just py
    just rs
    just r
    just go
alias mac := all-macos

all-ubuntu:
    apt-get update && \
    awk '/^ii/ { printf("apt-get install -y %s=%s\n", $2, $3) }' config/apt-get.lock > config/install_packages.sh && \
    chmod +x config/install_packages.sh && \
    config/install_packages.sh && \
    rm config/install_packages.sh && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*
    just py
    just rs
    just r
    just go
alias ub := all-ubuntu
