# genetraps-opensource

## Build status
[![Build Status](https://travis-ci.org/marpiech/genetraps-opensource.svg?branch=master)](https://travis-ci.org/marpiech/genetraps-opensource)

## genetraps system
# docker ports
8081 - client explorare
8082 - api explorae
8083 - client dnatoken
8084 - api dnatoken
8085 - client dx
8086 - api dx

## URL Tutorial

curl localhost:8080/hello

curl -X POST localhost:8080/upload?url=http://resources.intelliseq.pl/kamilant/test-data/fastq/capn3.1.fq.gz&tag=left \
**To check job id, and get file id** \
curl localhost:8080/describe/{**job1id**}

curl -X POST localhost:8080/upload?url=http://resources.intelliseq.pl/kamilant/test-data/fastq/capn3.2.fq.gz&tag=right \
curl localhost:8080/describe/{**job2id**}

curl -X POST localhost:8080/fastqc?fileId=**file1id** \
curl localhost:8080/describe/{**job3id**}

curl -X POST localhost:8080/fastqc?fileId=**file2id** \
curl localhost:8080/describe/{**job4id**}

**To get output folder, get folder json from describe, and delete "/rawdata" at the end** \
curl localhost:8080/describe/{**file1id**} \
curl -X POST localhost:8080/bwa?left=**file1id**&right=**file2id**&outputFolder=**outputFolder**
