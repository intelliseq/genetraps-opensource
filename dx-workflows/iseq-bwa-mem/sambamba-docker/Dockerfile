FROM ubuntu:18.04
MAINTAINER Katarzyna Kolanek, https://github.com/lltw

ADD https://github.com/biod/sambamba/releases/download/v0.6.7/sambamba_v0.6.7_linux.tar.bz2 .

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y && \
    tar -xvf sambamba_v0.6.7_linux.tar.bz2 && \
    rm sambamba_v0.6.7_linux.tar.bz2 && \
    ln -s /sambamba /usr/bin/sambamba
