### Endpoints
#### api-dx endpoints
[api-dx/hello](#api-dx-hello)  
[api-dx/mkdir](#api-dx-mkdir)  
[api-dx/upload](#api-dx-upload)  
[api-dx/describe](#api-dx-describe)  
[api-dx/fastqc](#api-dx-fastqc)  
[api-dx/bwa](#api-dx-bwa)  
[api-dx/gatkhc](#api-dx-gatkhc)  
[api-dx/ls](#api-dx-sample-ls)  
[api-dx/uploadfile](#api-dx-upload-file)  

----
#### api-dx hello
```
genetraps.intelliseq.pl:8086/hello
```
> Checks if server is up
* **URL:** `/hello`
* **Method:** `GET`
* **Success Response:**
  * **Code:** `200\`
  * **Content:** `{"status": "up"}`

----
#### api-dx mkdir
```
genetraps.intelliseq.pl:8086/mkdir
```
> Creates directory for a new sample
* **URL** `/mkdir`
* **Method:** `GET`
* **Success Response:**
  * **Code:** `200\`
  * **Content:** `{"response": 5}`

----
#### api-dx upload
```
genetraps.intelliseq.pl:8086/mkdir
```
> Uploads a file to dnanexus server
* **URL** `/upload`
* **Method:** `POST`
*  **URL Params**
   **Required:**  
   * `url=[string]`  
   * `sampleid=[string]`  
   **Optional**  
   * `tag=[string]` - can be used multiple times
* **Success Response:**
  * **Code:** `200`  
  * **Content:** `{"id": "job-XXXXXXXXXXXXXXXXXXXXXXXX"}`

----
#### api-dx describe
```
genetraps.intelliseq.pl:8086/describe
```
> Get a full description of file/job/... in JSON
* **URL** `/describe/:id`
* **Method:** `GET`
*  **URL Params**
   **Required:**
   * `id=[string]`
* **Success Response:**
  * **Code:** `200`
  * **Content:** huge json documented [here](https://wiki.dnanexus.com/API-Specification-v1.0.0/Applets-and-Entry-Points#API-method%3A-%2Fjob-xxxx%2Fdescribe)

----
#### api-dx fastqc
```
genetraps.intelliseq.pl:8086/fastqc
```
> Make fastqc job analysis
* **URL** `/fastqc`
* **Method:** `POST`
*  **URL Params**
   **Required:**
   * `fileId=[string]`
* **Success Response:**
  * **Code:** `200`
  * **Content:** `{"id": "job-XXXXXXXXXXXXXXXXXXXXXXXX"}`

----
#### api-dx bwa
```
genetraps.intelliseq.pl:8086/fastqc
```
> Make bwa job analysis
* **URL** `/bwa`
* **Method:** `POST`
*  **URL Params**
   **Required:**
   * `sampleid=[int]`  
   ***Or:***
   * `fastq_file_1(id)=[string]`
   * `fastq_file_2(id)=[string]`
* **Success Response:**
  * **Code:** `200`
  * **Content:** `{"id": "job-XXXXXXXXXXXXXXXXXXXXXXXX"}`

----
#### api-dx gatkhc
```
genetraps.intelliseq.pl:8086/gatkhc
```
>  Make gatk hc job analysis
* **URL** '/gatkhc'
* **Method:** `POST`
*  **URL Params**
   **Required:**
   * `sampleid=[string]`  
   **Not required:**
   * `interval=[string]`
* **Success Response:**
  * **Code:** `200`
  * **Content:** `{"id": "job-XXXXXXXXXXXXXXXXXXXXXXXX"}`

----
### api-dx upload file
```
genetraps.intelliseq.pl:8086/uploadfile
```
> Upload a given file, all params must be passed using '-F' (form) option.  
* **URL**  
  * `/uploadfile`  
* **Method:** `POST`  
*  **URL Params**  
   **Required:**  
   * `file=@[filename]`  
   * `sampleid=[string]`  
   **Not required:**  
   * `newfilename=[string]`
   * `tag=[string]` (can be multiple)  
* **Success Response:**  
  * **Code:** 200  
  * **Content:**    
Sample ls  
`{"id":"file-XXXXXXXXXXXXXXXXXXXXXXXX"}`  

----
### api-dx sample ls
```
genetraps.intelliseq.pl:8086/sample/:sampleid/ls
```
> List out the contents of a sample's rawdata folder by file ids ('byNames' is false by default). If optional parameter 'byNames' is set to true, listing proceeds by file names.  
* **URL**  
* `/sample/:sampleid/ls`  
* **Method:** `GET`  
*  **URL Params**  
 **Required:**  
 * `sampleid=[string]`  
 **Not required:**  
 * `byNames=[boolean]`  
* **Success Response:**  
* **Code:** 200  
* **Content:**    
Sample ls  
`{"file-XXXXXXXXXXXXXXXXXXXXXXXX":{"fileName":"file.fq.gz","tags":["tag"]}}`  
Sample ls (with byNames as true)  
`{"file.fq.gz":{"fileId":"file-XXXXXXXXXXXXXXXXXXXXXXXX","tags":["tag"]}}`

----
### api-dx properties (POST)
```
genetraps.intelliseq.pl:8086/sample/:sampleid/properties
```
> Sets or adds new properties to a sample and returns eventual map of properties (the sample must exist, otherwise returns exception) (returns exception if any property already exists)  
* **URL**  
  * `/sample/:sampleid/properties`  
* **Method:** `POST`  
*  **URL Params**  
   **Required:**  
   * `sampleid=[string]`  
   * `properties=[map]`  
* **Success Response:**  
  * **Code:** 200  
  * **Content:**    
Post properties  
`{"key1":"value1","key2":"value2",...}`  

----
### api-dx properties (GET)
```
genetraps.intelliseq.pl:8086/sample/:sampleid/properties
```
> Returns the map of properties (the sample must exist, otherwise returns exception) (returns exception if t)  
* **URL**  
  * `/sample/:sampleid/properties`  
* **Method:** `GET`  
*  **URL Params**  
   **Required:**  
   * `sampleid=[string]`   
* **Success Response:**  
  * **Code:** 200  
  * **Content:**    
Get properties  
`{"key1":"value1","key2":"value2",...}`  

----
### api-dx properties (PUT)
```
genetraps.intelliseq.pl:8086/sample/:sampleid/properties
```
> Changes values of already existing, specified properties of a sample and returns eventual map of properties (the sample must exist, otherwise returns exception) (returns exception if any of properties doesn't exists)  
* **URL**  
  * `/sample/:sampleid/properties`  
* **Method:** `PUT`  
*  **URL Params**  
   **Required:**  
   * `sampleid=[string]`  
   * `properties=[map]`  
* **Success Response:**  
  * **Code:** 200  
  * **Content:**    
Put properties  
`{"key1":"value1","key2":"value2",...}`    

----
### api-dx properties (DELETE)
```
genetraps.intelliseq.pl:8086/sample/:sampleid/properties
```
> Deletes already existing, specified properties of a sample and returns eventual map of properties (the sample must exist, otherwise returns exception) (returns exception if any of properties doesn't exists)  
* **URL**  
  * `/sample/:sampleid/properties`  
* **Method:** `DELETE`  
*  **URL Params**  
   **Required:**  
   * `sampleid=[string]`  
   * `properties=[map]`  
* **Success Response:**  
  * **Code:** 200  
  * **Content:**    
Delete properties  
`{"key1":"value1","key2":"value2",...}`    
