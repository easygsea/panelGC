# Base Image
FROM rocker/r-ver:4.3

# Metadata
LABEL base_image="rocker/r-ver:4.3"
LABEL version="1"
LABEL software="panelGC"
LABEL software.version="1.0.0"
LABEL about.summary="an open source tool for quantifying and monitoring GC bias in sequencing panels"
LABEL about.home="https://github.com/easygsea/panelGC"
LABEL about.license="SPDX:GPL-3.0"
LABEL about.tags="GC bias"
LABEL maintainer="Murathan Goktas <mgoktas@bcgsc.ca>"

RUN apt-get update \
        && apt-get install -y --no-install-recommends apt-utils \
        && apt-get install -y --no-install-recommends \
        libcurl4-openssl-dev \
        ## Basic deps
        gcc \
        gdb \
        libxml2-dev \
        libz-dev \
        liblzma-dev \
        libbz2-dev \
        libgit2-dev \
        ## sys deps from bioc_full
        pkg-config \
        fortran77-compiler \
        byacc \
        automake \
        curl \
        ## bedtools
        bedtools \
        && apt-get clean

RUN install2.r \
        --error \
        BiocManager \
        argparser \
        data.table \
        readxl \
        tidyverse

RUN Rscript -e 'requireNamespace("BiocManager"); BiocManager::install("GenomicRanges"); BiocManager::install("rtracklayer");'
