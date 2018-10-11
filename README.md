# genetraps-opensource
System  for automated WGS and WES analysis
copyright by intelliseq

## Build status
Master [![Build Status](https://travis-ci.org/intelliseq/genetraps-opensource.svg?branch=master)](https://travis-ci.org/intelliseq/genetraps-opensource)\
Develop [![Build Status](https://travis-ci.org/intelliseq/genetraps-opensource.svg?branch=develop)](https://travis-ci.org/intelliseq/genetraps-opensource)

## Services status
[![CLIENT-INDEX NOT ONLINE](http://genetraps.intelliseq.pl/client--index-online-brightgreen.svg)](http://genetraps.intelliseq.pl)\
[![API-DX NOT ONLINE](http://genetraps.intelliseq.pl:8086/status?v=1)](http://genetraps.intelliseq.pl:8086/hello)\
[![API-SECURITY NOT ONLINE](http://genetraps.intelliseq.pl:8088/status)](http://genetraps.intelliseq.pl:8088/hello)

## Usage
Build project
```
gradle build docker -p api-dx -x test
docker run -d -p 8080:8080 -e "DNANEXUS_TOKEN="$DNANEXUS_TOKEN -t pl.intelliseq.genetraps.api.dx/api-dx:latest
./scripts/wait-for-service.sh localhost:8080/hello 60
```

To check if server is up
```
curl localhost:8086/hello
```

**!!!**To run any request below you need an authorisation token to authorize
```
TOKEN=$(curl -XPOST "web_app:secret@localhost:8088/oauth/token" -d grant_type=password -d client_id=web_app -d username=$STAGING_USERNAME -d password=$STAGING_PASSWORD | jq -r ".access_token")
```

To create lowest sample folder and get number
```
SAMPLEID=$(curl -H "Authorization: Bearer $TOKEN" localhost:8086/mkdir | jq -r ".response")
```

To upload file to sample
```
U1ID=$(curl -X POST -H "Authorization: Bearer $TOKEN" "localhost:8086/upload?url=http://rawgit.com/intelliseq/genetraps-opensource/master/test-data/fastq/capn3.1.fq.gz&tag=left&sampleid=$SAMPLEID" | jq -r ".id")
```

To check job id, and get file id
```
curl -H "Authorization: Bearer $TOKEN" "localhost:8086/describe/$U1ID"
curl -H "Authorization: Bearer $TOKEN" "localhost:8086/describe/$U1ID" | jq -r ".state"
```

To get file id if job is done
```
F1ID=$(curl -H "Authorization: Bearer $TOKEN" "localhost:8086/describe/$U1ID" | jq -r '.output.file | .["$dnanexus_link"]')
```

To perform fastqc analysis
```
FQC1=$(curl -X POST -H "Authorization: Bearer $TOKEN" localhost:8086/fastqc?fileId=$F1ID | jq -r ".id")
curl -H "Authorization: Bearer $TOKEN" "localhost:8086/describe/$FQC1"
curl -H "Authorization: Bearer $TOKEN" "localhost:8086/describe/$FQC1" | jq -r ".state"
```

To run bwa
*reference (in properties): genetraps-resources:reference/grch38-no-alt/grch38-no-alt.tar*  
With tags: "left" i "right"  
```
BWAID=$(curl -X POST -H "Authorization: Bearer $TOKEN" "localhost:8086/bwa?sampleid=$SAMPLEID" | jq -r ".id")
```
With files' ids:  
```
BWAID=$(curl -X POST -H "Authorization: Bearer $TOKEN" "localhost:8086/bwa?fastq_file_1=$F1ID&fastq_file_2=$F2ID" | jq -r ".id")
```

To run gatk hc
*reference and vcfs (their ids) in properties; interval not required*  
```
GATKID=$(curl -X POST -H "Authorization: Bearer $TOKEN" "localhost:8086/gatkhc?sampleid=$SAMPLEID&interval=chr15:42377802-42397802" | jq -r ".id")
```

To list out contents of a given sample's rawdata folder
```
curl -H "Authorization: Bearer $TOKEN" localhost:8086/sample/$SAMPLEID/ls
```
To list out contents of a given sample's rawdata folder (in reverse order)
```
curl -H "Authorization: Bearer $TOKEN" localhost:8086/sample/$SAMPLEID/revls
```

## Ports
```
8081 - client explorare
8082 - api explorare
8083 - client dnatoken
8084 - api dnatoken
8085 - client dx
8086 - api dx
8087 - client index
8088 - api security
```

## Project Setup
See `.travis.yml` and `travis-script.sh` for instructions

## Endpoints
[Endpoints](https://github.com/intelliseq/genetraps-opensource/blob/master/docs/endpoints.md)
