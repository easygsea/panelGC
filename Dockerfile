# Base Image
FROM rocker/r-ver:4.3

# Metadata
LABEL base_image="rocker/r-ver:4.3"
LABEL version="1"
LABEL software="panelGC"
LABEL software.version="1.2.0"
LABEL about.summary="An open source tool for quantifying and monitoring GC bias in sequencing panels"
LABEL about.home="https://github.com/easygsea/panelGC"
LABEL about.license="SPDX:GPL-3.0"
LABEL about.tags="GC bias"

# System dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    apt-utils \
    libcurl4-openssl-dev \
    gcc gdb libxml2-dev libz-dev liblzma-dev libbz2-dev libgit2-dev \
    pkg-config fortran77-compiler byacc automake curl bedtools \
    openjdk-11-jre-headless git unzip tar bzip2 \
    # for samtools tview
    libncurses5-dev libncursesw5-dev \
    && rm -rf /var/lib/apt/lists/* \
    && apt-get clean

# Install HTSlib and samtools
ENV HTSLIB_VERSION=1.19
RUN curl -sL "https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2" -o "htslib-${HTSLIB_VERSION}.tar.bz2" && \
    tar -xjf "htslib-${HTSLIB_VERSION}.tar.bz2" && \
    cd "htslib-${HTSLIB_VERSION}" && \
    ./configure --prefix=/usr/local && \
    make && make install && \
    cd .. && \
    curl -sL "https://github.com/samtools/samtools/releases/download/${HTSLIB_VERSION}/samtools-${HTSLIB_VERSION}.tar.bz2" -o "samtools-${HTSLIB_VERSION}.tar.bz2" && \
    tar -xjf "samtools-${HTSLIB_VERSION}.tar.bz2" && \
    cd "samtools-${HTSLIB_VERSION}" && \
    ./configure --prefix=/usr/local && \
    make && make install && \
    cd .. && \
    rm -rf htslib-* samtools-*

# R packages
RUN install2.r --error argparser tidyverse data.table
