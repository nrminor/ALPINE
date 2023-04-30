FROM ubuntu:20.04

# Install required packages
RUN apt-get update && apt-get install -y muscle python3-pip

# Install Biopython
RUN pip3 install biopython

# Set working directory
WORKDIR /app