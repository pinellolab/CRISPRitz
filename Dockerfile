# CRISPRitz standalone image — Python 3.11.
#
# This image no longer `conda install`s the published crispritz package. That
# package was pinned to python=3.8, which is unsatisfiable against the current
# continuumio/miniconda3 base and broke the build on BOTH linux/amd64 and
# linux/arm64 at `conda install python=3.8` (see the history in
# .github/workflows/docker.yml).
#
# Instead this builds the C++ binaries (buildTST / searchTST / searchBruteForce)
# from THIS source tree and lays them out in the layout the root crispritz.py
# expects for a non-source install: with crispritz.py in <prefix>/bin and the
# binaries + Python scoring scripts in <prefix>/opt/crispritz (crispritz.py
# computes that as os.path.dirname(__file__)[:-3] + "opt/crispritz/").
#
# The Python runtime is 3.11 and the scoring stack is pinned to
# scikit-learn 1.1.* + numpy<2: the vendored (CRISPR-HAWK) azimuth Doench-2016
# model only unpickles on that combo. See
# sourceCode/Python_Scripts/Scores/azimuth/ATTRIBUTION.md. bcftools + bedtools
# (needed by the scoring / VCF paths) come from Debian.
#
# Builds and runs natively on linux/amd64 and linux/arm64 (Apple Silicon).
FROM python:3.11-slim

LABEL maintainer="Pinello Lab"

ENV DEBIAN_FRONTEND=noninteractive
SHELL ["/bin/bash", "-c"]

# System build + runtime deps: g++/OpenMP + GSL to compile the C++ binaries,
# and bcftools + bedtools for the VCF and scoring (getfasta) paths.
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        g++ \
        make \
        libgomp1 \
        libgsl-dev \
        bcftools \
        bedtools \
        procps \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Python scoring/analysis stack. scikit-learn is pinned to 1.1.* and numpy to
# <2 so the vendored azimuth model unpickles; do NOT bump sklearn past 1.1.x
# without re-serializing the model (see azimuth/ATTRIBUTION.md).
RUN pip install --no-cache-dir \
        'scikit-learn==1.1.3' \
        'numpy==1.24.4' \
        'scipy==1.10.1' \
        'pandas==2.0.3' \
        biopython \
        more-itertools \
        statsmodels \
        matplotlib \
        intervaltree

# Build the C++ binaries from source. buildTST is compiled from mainParallel.cpp
# (the multi-threaded index builder); searchTST and searchBruteForce from their
# respective translation units.
WORKDIR /opt/crispritz-src
COPY . /opt/crispritz-src

RUN g++ -std=c++11 -O3 -fopenmp \
        sourceCode/CRISPR-Cas-Tree/mainParallel.cpp -o buildTST \
    && g++ -std=c++11 -O3 -fopenmp \
        sourceCode/CRISPR-Cas-Tree/searchOnTST.cpp \
        sourceCode/CRISPR-Cas-Tree/detailedOutput.cpp \
        sourceCode/CRISPR-Cas-Tree/convert.cpp \
        -I sourceCode/CRISPR-Cas-Tree/include -o searchTST \
    && g++ -std=c++11 -O3 -fopenmp \
        sourceCode/CRISPRofiler/main.cpp \
        sourceCode/CRISPRofiler/profiling.cpp \
        sourceCode/CRISPRofiler/guide_searching.cpp \
        sourceCode/CRISPRofiler/pam_searching.cpp \
        sourceCode/CRISPRofiler/pre_computation.cpp \
        sourceCode/CRISPRofiler/reading.cpp \
        sourceCode/CRISPRofiler/analysis.cpp -o searchBruteForce

# Install into the layout crispritz.py resolves in a non-source install:
#   /usr/local/bin/crispritz.py           -> executable on PATH
#   /usr/local/opt/crispritz/{buildTST,searchTST,searchBruteForce}
#   /usr/local/opt/crispritz/Python_Scripts/...
# crispritz.py derives the opt path as dirname(__file__)[:-3] + "opt/crispritz/"
# i.e. "/usr/local/bin"[:-3] -> "/usr/local/" + "opt/crispritz/".
RUN mkdir -p /usr/local/opt/crispritz \
    && cp crispritz.py /usr/local/bin/crispritz.py \
    && chmod +x /usr/local/bin/crispritz.py \
    && cp buildTST searchTST searchBruteForce /usr/local/opt/crispritz/ \
    && cp -R sourceCode/Python_Scripts /usr/local/opt/crispritz/

WORKDIR /root

CMD ["crispritz.py"]
