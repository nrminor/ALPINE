# Use an official R runtime as a parent image
FROM r-base:latest

# Set the working directory to /app
WORKDIR /app

# Install required R packages
RUN R -e "install.packages('tidyverse', repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('ape', repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('adephylo', repos='http://cran.us.r-project.org')"
