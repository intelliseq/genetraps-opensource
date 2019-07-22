FROM ubuntu:18.04

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y  software-properties-common && \
    apt-get install -y  gradle && \
    add-apt-repository ppa:linuxuprising/java -y && \
    apt-get update && \
    apt-get install -y openjdk-11-jdk && \
    apt-get clean && \
    apt-get install -y python-pip python-dev build-essential && \
    apt-get install -y python-dev libffi-dev libssl-dev && \
    pip install --upgrade --force-reinstall pip==9.0.3 && \
    pip install dxpy --upgrade

RUN echo "echo logging into dxpy" > start.sh
RUN echo "dx login --token "'$'"DNANEXUS_TOKEN_TEST --noprojects" >> start.sh
RUN echo "dx select genetraps-test" >> start.sh
RUN echo "echo starting app" >> start.sh
RUN echo "java -jar app.jar" >> start.sh

ADD build/libs/api-dx-latest.jar app.jar
ADD src/main/resources/ src/main/resources/

ENTRYPOINT ["sh", "./start.sh"]
