FROM ubuntu:18.04
MAINTAINER Marcin Piechota, https://github.com/marpiech

ENV package_name fastqc_v0.11.7.zip

ADD http://www.bioinformatics.babraham.ac.uk/projects/fastqc/${package_name} /tmp/

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y  software-properties-common && \
    add-apt-repository ppa:webupd8team/java -y && \
    apt-get update && \
    echo oracle-java8-installer shared/accepted-oracle-license-v1-1 select true | /usr/bin/debconf-set-selections && \
    apt-get install -y oracle-java8-installer && \
    apt-get install -y unzip && \
    apt-get install -y perl && \
    apt-get clean && \
    cd /opt/ && unzip /tmp/${package_name} \
    && rm /tmp/${package_name} \
    && chmod 750 /opt/FastQC/fastqc

CMD ["-help"]
ENTRYPOINT ["/opt/FastQC/fastqc"]
