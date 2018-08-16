#!/bin/bash

### VERSIONS ###
ISEQ_FASTQC_VERSION=0.12
ISEQ_BWA_VERSION=0.1

### FUNCTIONS ###
function docker_tag_exists() {
    curl --silent -f -lSL https://index.docker.io/v1/repositories/$1/tags/$2 > /dev/null
}

### TOOLS ###
# DOCKERHUB #
echo $DOCKERHUB_PASSWORD | docker login -u $DOCKERHUB_USERNAME --password-stdin
# DXWDL #
if [ ! -f /tmp/dxWDL/dxWDL-0.65.jar ]; then
    wget https://github.com/dnanexus/dxWDL/releases/download/0.65/dxWDL-0.65.jar --directory-prefix /tmp/dxWDL
fi

### FASTQC ###

if docker_tag_exists intelliseq/fastqc $ISEQ_FASTQC_VERSION; then
    echo "Fastqc without changes..."
else
    cat fastqc/Dockerfile | docker build - -t intelliseq/fastqc:$ISEQ_FASTQC_VERSION
    docker tag intelliseq/fastqc:$ISEQ_FASTQC_VERSION intelliseq/fastqc:latest
    # test
    docker run -v "$PWD"/../../test-data/fastq/:/tmp/data -v /tmp/:/tmp/out intelliseq/fastqc -o /tmp/out /tmp/data/capn3.1.fq.gz
    java -jar /tmp/dxWDL/dxWDL-0.65.jar compile wdl-fastqc/iseq-fastqc.wdl -f
fi

### BWA ###

#if docker_tag_exists intelliseq/bwa $ISEQ_BWA_VERSION; then
#    echo "BWA without changes..."
#else
#    cat bwa/Dockerfile | docker build - -t intelliseq/bwa:$ISEQ_BWA_VERSION
#    docker tag intelliseq/bwa:$ISEQ_BWA_VERSION intelliseq/bwa:latest
#    # test
#    docker run -v "$PWD"/../../test-data/fastq/:/tmp/data -v /tmp/:/tmp/out intelliseq/fastqc -o /tmp/out /tmp/data/capn3.1.fq.gz
#    java -jar /tmp/dxWDL/dxWDL-0.65.jar compile wdl-fastqc/iseq-fastqc.wdl -f
#fi
