# Set the base image to anaconda python 3.6
FROM continuumio/miniconda3

# File Author / Maintainer
MAINTAINER Samuele Cancelleri

ENV SHELL bash

RUN conda config --add channels defaults
RUN conda config --add channels conda-forge
RUN conda config --add channels bioconda

RUN apt-get update && apt-get install gsl-bin libgsl0-dev -y && apt-get install libgomp1 -y && apt-get install g++-5 -y && apt-get clean

#Install crispritz package
RUN conda install crispritz -y && conda clean --all -y

