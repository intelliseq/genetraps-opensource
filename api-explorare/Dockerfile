FROM ubuntu:16.04

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y  software-properties-common && \
    add-apt-repository ppa:linuxuprising/java -y && \
    apt-get update && \
    echo oracle-java10-installer shared/accepted-oracle-license-v1-1 select true | /usr/bin/debconf-set-selections && \
    apt-get install -y oracle-java10-installer && \
    apt-get clean

#FROM openjdk:8-jdk-alpine

RUN echo "echo starting app" > start.sh
RUN echo "java -jar app.jar" >> start.sh
ARG JAR_FILE
ADD ${JAR_FILE} app.jar
ENTRYPOINT ["sh", "./start.sh"]
#ENTRYPOINT ["java","-Djava.security.egd=file:/dev/./urandom","-jar","/app.jar"]
