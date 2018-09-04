FROM ubuntu:18.04
MAINTAINER Marcin Piechota, https://github.com/marpiech

ENV package_name bwa-0.7.17

ADD http://downloads.sourceforge.net/project/bio-bwa/${package_name}.tar.bz2 /tmp/

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y build-essential && \
    apt-get install -y software-properties-common && \
    apt-get install -y libz-dev && \
    apt-get install -y libncurses5-dev libncursesw5-dev && \
    cd /tmp/ && tar xjvf ${package_name}.tar.bz2 && \
    cd /tmp/${package_name} && \ 
    make && \
    mv /tmp/${package_name}/bwa /usr/bin

ENTRYPOINT ["/usr/bin/bwa"]
