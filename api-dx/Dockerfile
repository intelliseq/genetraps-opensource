FROM ubuntu:18.04

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y  software-properties-common && \
    apt-get install -y  gradle && \
    add-apt-repository ppa:linuxuprising/java -y && \
    apt-get update && \
    apt-get install -y openjdk-11-jdk

ARG BUCKET_TAG
ARG QUEUE_TAG

ADD "http://resouces.intelliseq.pl/public/intelliseqngs/workflows/cromwell/cromwell.jar" "cromwell.jar"
ADD "http://resouces.intelliseq.pl/public/intelliseqngs/workflows/cromwell/aws.cfg" "aws.cfg"

RUN sed -i "s+BUCKET_TAG+${BUCKET_TAG}+" aws.cfg
RUN sed -i "s+QUEUE_TAG+${QUEUE_TAG}+" aws.cfg

RUN echo "free -m" > start.sh
RUN echo "echo starting cromwell server" >> start.sh
RUN echo "java -Dconfig.file=aws.cfg -jar cromwell.jar server &" >> start.sh
RUN echo "echo starting app" >> start.sh
RUN echo "java -jar app.jar" >> start.sh

ADD build/libs/api-dx-latest.jar app.jar
ADD src/main/resources/ src/main/resources/

ENTRYPOINT ["sh", "./start.sh"]
