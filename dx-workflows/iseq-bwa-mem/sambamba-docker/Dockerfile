FROM ubuntu:18.04
MAINTAINER Katarzyna Kolanek, https://github.com/lltw
MAINTAINER Marcin piechota, https://github.com/marpiech

ADD https://github.com/biod/sambamba/releases/download/v0.6.7/sambamba_v0.6.7_linux.tar.bz2 .

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y \ 
    bzip2 && \
    tar -xvf sambamba_v0.6.7_linux.tar.bz2 && \
    rm sambamba_v0.6.7_linux.tar.bz2 && \
    ln -s /sambamba /usr/bin/sambamba
