# genetraps-opensource

## Build status
[![Build Status](https://travis-ci.org/marpiech/genetraps-opensource.svg?branch=master)](https://travis-ci.org/marpiech/genetraps-opensource)

## Services status
### client-index service
![CLIENT-INDEX NOT ONLINE](http://genetraps.intelliseq.pl/client--index-online-brightgreen.svg)
### api-dx service
![API-DX NOT ONLINE](http://genetraps.intelliseq.pl/-client--index-online-brightgreen.svg)
### api-security
![API-SECURITY NOT ONLINE](http://genetraps.intelliseq.pl/-client--index-online-brightgreen.svg)

## development philosophy
# API dx
Project id should be determined by token.
As docker image port should be set to 8086.

## genetraps system
# docker ports
#### 8081 - client explorare
#### 8082 - api explorare
#### 8083 - client dnatoken
#### 8084 - api dnatoken
#### 8085 - client dx
#### 8086 - api dx
#### 8087 - client index
#### 8088 - api security

## URL Tutorial

`gradle build docker -p api-dx -x test`
`docker run -d -p 8080:8080 -e "DNANEXUS_TOKEN="$DNANEXUS_TOKEN -t pl.intelliseq.genetraps.api.dx/api-dx:latest`  
`./scripts/wait-for-service.sh localhost:8080/hello 60`  


**To check if server is up**  
`curl localhost:8080/hello`  

**To create lowest sample folder and get number**  
`SAMPLE_NUMBER=$(curl localhost:8080/mkdir | jq -r ".response")`  

**To upload file to sample**  
`U1ID=$(curl -X POST "localhost:8080/upload?url=http://resources.intelliseq.pl/kamilant/test-data/fastq/capn3.1.fq.gz&tag=left&sampleNumber=$SAMPLE_NUMBER" | jq -r ".id")`  

**To check job id, and get file id**  
`curl "localhost:8080/describe/$U1ID"`
`curl "localhost:8080/describe/$U1ID" | jq -r ".state"`

**To get file id if job is done**
`F1ID=$(curl "localhost:8080/describe/$U1ID" | jq -r '.output.file | .["$dnanexus_link"]')`

**To make fastqc analysis**
`FQC1=$(curl -X POST localhost:8080/fastqc?fileId=$F1ID | jq -r ".id")`
`curl "localhost:8080/describe/$FQC1"`
`curl "localhost:8080/describe/$FQC1" | jq -r ".state"`

# Project Setup
## Set environment variables
```
AWS_ACCESS_KEY_ID  
AWS_PROFILE  
AWS_REGION  
AWS_SECRET_ACCESS_KEY
```

## Set up tools
```
curl -s "https://get.sdkman.io" | bash
source "/home/marpiech/.sdkman/bin/sdkman-init.sh"
sdk install gradle 4.6
```

# Endpoints

**hello**
----
  Simple check if server is up

* **URL**

  /hello

* **Method:**

  `GET`

* **Success Response:**

  * **Code:** 200
    **Content:** `{"status": "up"}`

**mkdir**
----
  Creates directory for a new sample

* **URL**

  /mkdir

* **Method:**

  `GET`

* **Success Response:**

  * **Code:** 200
    **Content:** `{"response": 5}`

**upload**
----
  Upload a file to dnanexus server

* **URL**

  /upload

* **Method:**

  `POST`

*  **URL Params**

   **Required:**

   `url=[string]`
   `sampleNumber=[string]`

   **Optional**

   `tag=[string]` - can be used multiple times

* **Success Response:**

  * **Code:** 200
    **Content:** `{"id": "job-XXXXXXXXXXXXXXXXXXXXXXXX"}`

**describe**
----
  Get a full description of file/job/... in JSON

* **URL**

  /describe/:id

* **Method:**

  `GET`

*  **URL Params**

   **Required:**

   `id=[string]`

* **Success Response:**

  * **Code:** 200
    **Content:** huge json documented [here](https://wiki.dnanexus.com/API-Specification-v1.0.0/Applets-and-Entry-Points#API-method%3A-%2Fjob-xxxx%2Fdescribe)


**fastqc**
----
  Make fastqc job analysis

* **URL**

  /fastqc

* **Method:**

  `POST`

*  **URL Params**

   **Required:**

   `fileId=[string]`

* **Success Response:**

  * **Code:** 200
    **Content:** `{"id": "job-XXXXXXXXXXXXXXXXXXXXXXXX"}`
