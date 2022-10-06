# syntax=docker/dockerfile:1

FROM rocker/tidyverse

MAINTAINER plutzer@wustl.edu

COPY r-installation.R r-installation.R
RUN Rscript r-installation.R

WORKDIR /work

ADD ./build/precompiled_linux/SAINTexpress-int /bin

RUN apt-get update && apt-get install -y python3-pip

COPY requirements.txt requirements.txt
RUN pip3 install -r requirements.txt

COPY *.R /
COPY *.py /

CMD /bin/bash