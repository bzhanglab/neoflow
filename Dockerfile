#specifying the base image
FROM java:8
#FROM ubuntu:latest
FROM biocontainers/biocontainers:v1.0.0_cv4

MAINTAINER Kai Li <kail@bcm.edu>

USER root


RUN apt-get update \
  && apt-get install -y python3-pip python3-dev 

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
#install pepquery-1.4
RUN wget http://www.pepquery.org/data/PepQuery_v1.3.0.tar.gz \
	&& tar -xzvf PepQuery_v1.3.0.tar.gz \
	&& rm PepQuery_v1.3.0.tar.gz \
	&& mv PepQuery_v1.3.0 /opt/PepQuery_v1.3.0/

#chmod of /home/user and change working directory
RUN chmod -R 755 /opt/

#change working directory
WORKDIR /opt/

#specify the command executed when the container is started
CMD ["/bin/bash"]
