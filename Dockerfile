# Set the base image to anaconda python 3.6
FROM continuumio/miniconda3

# File Author / Maintainer
MAINTAINER Samuele Cancelleri

ENV SHELL bash

#update conda channel with bioconda and conda-forge
RUN conda config --add channels defaults 
RUN conda config --add channels conda-forge
RUN conda config --add channels bioconda

#update packages of the docker system
RUN apt-get update && apt-get install gsl-bin libgsl0-dev -y && apt-get install libgomp1 -y && apt-get clean

#Install crispritz package (change the dafault version of python to 3.6)
RUN conda update -n base -c defaults conda
RUN conda install python=3.8 -y
RUN conda install xlrd -y
RUN conda install crispritz -y && conda clean --all -y
