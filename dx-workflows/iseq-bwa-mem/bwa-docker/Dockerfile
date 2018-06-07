FROM ubuntu:18.04
MAINTAINER Katarzyna Kolanek, https://github.com/lltw

ADD https://github.com/GregoryFaust/samblaster/releases/download/v.0.1.24/samblaster-v.0.1.24.tar.gz .
ADD https://github.com/samtools/samtools/releases/download/1.8/samtools-1.8.tar.bz2 .
ADD https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2 .


RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y \
    build-essential \
    libz-dev \
    libncurses5-dev \
    libbz2-dev \
    liblzma-dev && \
    tar -xvf bwa-0.7.17.tar.bz2&& \
    rm bwa-0.7.17.tar.bz2 && \
    make -C /bwa-0.7.17 && \
    ln -s /bwa-0.7.17/bwa /usr/bin/bwa && \
    tar -xvf samtools-1.8.tar.bz2 && \
    rm samtools-1.8.tar.bz2 && \
    make -C /samtools-1.8  && \
    ln -s /samtools-1.8/samtools /usr/bin/samtools && \
    tar -xvf samblaster-v.0.1.24.tar.gz && \
    rm samblaster-v.0.1.24.tar.gz && \
    make -C /samblaster-v.0.1.24  && \
    ln -s /samblaster-v.0.1.24/samblaster /usr/bin/samblaster
