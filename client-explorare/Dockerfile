FROM ubuntu:18.04

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y npm && \
    apt-get update && \
    apt-get clean


RUN echo "echo starting app" > start.sh
RUN echo "npm start" >> start.sh

ADD . / ./

ENTRYPOINT ["sh", "./start.sh"]
