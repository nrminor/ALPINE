FROM ubuntu:20.04

# Install required packages
RUN apt-get update && apt-get install -y wget gcc make zlib1g-dev muscle python3-pip

# Install Biopython
RUN pip3 install biopython

# Install pandas
RUN pip3 install pandas

# Install polars
RUN pip3 install polars

# Install seeqtk
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