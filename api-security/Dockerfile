FROM ubuntu:18.04

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y  software-properties-common && \
    apt-get install -y  gradle && \
    add-apt-repository ppa:linuxuprising/java -y && \
    apt-get update && \
    echo oracle-java11-installer shared/accepted-oracle-license-v1-2 select true | /usr/bin/debconf-set-selections && \
    apt-get install -y oracle-java11-installer && \
    apt-get clean


RUN echo "echo starting app" > start.sh
RUN echo "java -jar app.jar" >> start.sh

ADD build/libs/app.jar app.jar
ADD src/main/resources/ src/main/resources/

ENTRYPOINT ["sh", "./start.sh"]
