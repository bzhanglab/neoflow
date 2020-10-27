#specifying the base image
FROM ubuntu:19.10

MAINTAINER Kai Li <kail@bcm.edu>

#USER root

ARG PYTHON=python3
ARG PIP=pip3

# See http://bugs.python.org/issue19846
ENV LANG C.UTF-8

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    wget \
    tcsh \
    git \
    unzip \
    pkg-config \
    curl \
    ${PYTHON} \
    ${PYTHON}-pip && \
    apt-get -y  install --fix-missing openjdk-8-jre && \
    apt-get clean

RUN ${PIP} --no-cache-dir install --upgrade \
    pip \
    setuptools

# Some TF tools expect a "python" binary
RUN ln -s $(which ${PYTHON}) /usr/local/bin/python

RUN ${PIP} install numpy pandas biopython

#install MSGFPlus
RUN wget https://github.com/MSGFPlus/msgfplus/releases/download/v2019.07.03/MSGFPlus_v20190703.zip \
	&& unzip MSGFPlus_v20190703.zip \
	&& rm MSGFPlus_v20190703.zip \
	&& mv MSGFPlus.jar /opt/
#install X!Tandem
RUN wget ftp://ftp.thegpm.org/projects/tandem/source/tandem-linux-17-02-01-4.zip \
	&& unzip tandem-linux-17-02-01-4.zip \
	&& rm tandem-linux-17-02-01-4.zip \
	&& chmod 755 tandem-linux-17-02-01-4/bin/tandem.exe \
	&& mv tandem-linux-17-02-01-4 /opt/tandem-linux-17-02-01-4
#install Comet
RUN wget https://sourceforge.net/projects/comet-ms/files/comet_2018014.zip/download \
	&& mv download comet_2018014.zip \
	&& unzip comet_2018014.zip \
	&& rm comet_2018014.zip \
	&& mv comet.2018014.linux.exe /opt/ \
	&& chmod 755 /opt/comet.2018014.linux.exe
#install mzidlib-1.7
RUN wget http://www.proteoannotator.org/datasets/releases/ProteoAnnotator-1.7.86.zip \
	&& unzip ProteoAnnotator-1.7.86.zip \
	&& rm ProteoAnnotator-1.7.86.zip \
	&& mv ProteoAnnotator-1.7.86/mzidlib-1.7 /opt/mzidlib-1.7/
#install pepquery-1.3
RUN wget http://www.pepquery.org/data/PepQuery_v1.3.0.tar.gz \
	&& tar -xzvf PepQuery_v1.3.0.tar.gz \
	&& rm PepQuery_v1.3.0.tar.gz \
	&& mv PepQuery_v1.3.0 /opt/PepQuery_v1.3.0/

RUN wget -O /opt/customprodbj.jar https://github.com/wenbostar/Customprodbj/releases/download/v1.1.0/customprodbj.jar

RUN wget -O /opt/pepmap.jar https://github.com/wenbostar/pepmap/releases/download/v1.0.0/pepmap.jar

RUN git clone https://github.com/bzhanglab/neoflow/ /opt/neoflow/



#chmod of /home/user and change working directory
RUN chmod -R 755 /opt/

#change working directory
WORKDIR /opt/

#specify the command executed when the container is started
CMD ["/bin/bash"]
