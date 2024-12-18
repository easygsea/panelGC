# Base Image
FROM rocker/r-ver:4.3

# Metadata
LABEL base_image="rocker/r-ver:4.3"
LABEL version="1"
LABEL software="panelGC"
LABEL software.version="1.0.1"
LABEL about.summary="An open source tool for quantifying and monitoring GC bias in sequencing panels"
LABEL about.home="https://github.com/easygsea/panelGC"
LABEL about.license="SPDX:GPL-3.0"
LABEL about.tags="GC bias"
LABEL maintainer="Murathan Goktas <mgoktas@bcgsc.ca>"

# System dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    apt-utils \
    libcurl4-openssl-dev \
    gcc gdb libxml2-dev libz-dev liblzma-dev libbz2-dev libgit2-dev \
    pkg-config fortran77-compiler byacc automake curl bedtools \
    openjdk-11-jre-headless git unzip && apt-get clean

# R packages
RUN install2.r --error BiocManager argparser tidyverse
RUN Rscript -e 'requireNamespace("BiocManager"); BiocManager::install("GenomicRanges"); BiocManager::install("rtracklayer")'
# Install Nextflow
ENV NEXTFLOW_VERSION=23.04.0
RUN curl -s https://get.nextflow.io | bash && \
    mv nextflow /usr/local/bin/ && \
    chmod 755 /usr/local/bin/nextflow
