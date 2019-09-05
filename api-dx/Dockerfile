FROM ubuntu:18.04

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y  software-properties-common && \
    apt-get install -y  gradle && \
    add-apt-repository ppa:linuxuprising/java -y && \
    apt-get update && \
    apt-get install -y openjdk-11-jdk

ENV MYSQL_USER=mysql \
    MYSQL_VERSION=5.7.26 \
    MYSQL_DATA_DIR=/var/lib/mysql \
    MYSQL_RUN_DIR=/run/mysqld \
    MYSQL_LOG_DIR=/var/log/mysql

RUN apt-get update \
 && DEBIAN_FRONTEND=noninteractive apt-get install -y mysql-server \
 && rm -rf ${MYSQL_DATA_DIR} \
 && rm -rf /var/lib/apt/lists/*

COPY src/main/resources/bash/runmysql.sh /sbin/runmysql.sh

RUN chmod 755 /sbin/runmysql.sh

ARG BUCKET_TAG
ARG QUEUE_TAG

ADD "http://resouces.intelliseq.pl/public/intelliseqngs/workflows/cromwell/cromwell.jar" "cromwell.jar"
ADD "http://resouces.intelliseq.pl/public/intelliseqngs/workflows/cromwell/awsmysql.cfg" "aws.cfg"

RUN sed -i "s+BUCKET_TAG+${BUCKET_TAG}+" aws.cfg
RUN sed -i "s+QUEUE_TAG+${QUEUE_TAG}+" aws.cfg

RUN echo "echo starting mysql" > start.sh
RUN echo "/sbin/runmysql.sh" >> start.sh
RUN echo "echo starting cromwell server" >> start.sh
RUN echo "java -Dconfig.file=aws.cfg -jar cromwell.jar server &" >> start.sh
RUN echo "echo starting app" >> start.sh
RUN echo "java -jar app.jar" >> start.sh

ADD build/libs/api-dx-latest.jar app.jar
ADD src/main/resources/ src/main/resources/

ENTRYPOINT ["sh", "./start.sh"]
