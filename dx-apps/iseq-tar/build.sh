#!/bin/bash

### VERSIONS ###
ISEQ_TAR_VERSION=0.1

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

### TAR ###

if docker_tag_exists intelliseq/bwa $ISEQ_TAR_VERSION; then
    echo "BWA without changes..."
else
    cat Dockerfile | docker build - -t intelliseq/tar:$ISEQ_TAR_VERSION
    docker tag intelliseq/"tar":$ISEQ_TAR_VERSION intelliseq/"tar":latest
    # test
    docker run -v "$PWD"/../../test-data/fastq/:/tmp/data -v /tmp/:/tmp/out intelliseq/tar -cf /tmp/out/archive.tar /tmp/data/capn3.1.fq.gz
    java -jar /tmp/dxWDL/dxWDL-0.65.jar compile iseq-tar.wdl -f
fi
