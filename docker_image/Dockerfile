#Download base image ubuntu 16.04
FROM ubuntu:16.04

# Update Ubuntu Software repository
RUN apt-get update

# Upgrade Ubuntu Software repository
RUN apt-get upgrade -y

# Install all required libraries and packages
RUN apt-get install g++-5 -y
RUN apt-get install python3 -y
RUN apt-get install wget -y
RUN apt-get install unzip -y
RUN apt-get install libboost-all-dev -y
RUN apt-get install bcftools -y
RUN apt-get install curl -y
RUN apt-get install python3-tk -y
RUN apt-get install make -y

# Set case-insesitive bash
RUN echo set completion-ignore-case on | tee -a /etc/inputrc

# Set g++-5 default g++ and python3 default python
RUN update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-5 50
RUN update-alternatives --install /usr/bin/python python /usr/bin/python3.5 50

# Pip libraries for python
RUN curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
RUN python get-pip.py
RUN pip install intervaltree
RUN pip install matplotlib
RUN pip install pandas
RUN pip install scipy
RUN pip install numpy
RUN pip install more-itertools
RUN pip install ecdf
RUN pip install statsmodels

# Prepare CRISPRitz to be tested and used
RUN wget https://github.com/InfOmics/CRISPRitz/archive/master.zip
RUN unzip master.zip
RUN rm master.zip
RUN mv CRISPRitz-master/ CRISPRitz
#RUN bash CRISPRitz/download_test_files.sh
RUN mv CRISPRitz/ root/
#RUN mv chr22_* root/CRISPRitz/
RUN rm root/CRISPRitz/download_test_files.sh root/CRISPRitz/README.md
RUN chmod -R 700 root/CRISPRitz/*
RUN make -f root/CRISPRitz/docker_image/Makefile_docker
RUN mv root/CRISPRitz/docker_image/crispritz.py buildTST searchBruteForce searchTST /usr/bin/
RUN cp -R root/CRISPRitz/sourceCode/Python_Scripts/ /usr/bin/
#RUN rm root/CRISPRitz/Makefile root/CRISPRitz/meta.yaml root/CRISPRitz/build.sh root/CRISPRitz/LICENSE root/CRISPRitz/Dockerfile root/CRISPRitz/crispritz.py
#RUN rm -rf root/CRISPRitz/docker_image
#RUN rm -rf root/CRISPRitz/sourceCode
RUN rm -rf root/CRISPRitz

# Set workdir
WORKDIR root/CRISPRitz

# Set expose to port 80 and 443
EXPOSE 80 443


