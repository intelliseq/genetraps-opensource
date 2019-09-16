
# Documentation of using pipelines via genetraps API

**First and foremost, to use any of the following commands you need an authorization token:**
```bash
TOKEN=$(curl -XPOST "web_app:secret@genetraps.intelliseq.pl:8088/oauth/token" -d grant_type=password -d client_id=web_app -d username=$STAGING_USERNAME -d password=$STAGING_PASSWORD | jq -r ".access_token")
```

**You can generate the lowest available sample id (number) using the following command:**
```bash
SAMPLEID=$(curl -H "Authorization: Bearer $TOKEN" "genetraps.intelliseq.pl:8086/sample/create" | jq -r ".response")
```
or
**Using any given number:**
(at risk of using an existing one if unintentionally, although it won't delete the existing files in the sample, only mix in the new ones with them/overwrite some)
```bash
SAMPLEID=
```

**Link to your pipeline - obligatory:**
```bash
WDL=
```
(eg. https://gitlab.com/intelliseq/workflows/raw/master/src/main/wdl/tasks/minimal-task/v1.0/minimal.wdl)

**To upload an input file from your machine to aws (returns special aws path to the file):**
Special aws path can be recognized by '/' at the beginning (is followed by nr of sample, later dirs, ending with name of file)
```bash
curl -H "Authorization: Bearer $TOKEN" "genetraps.intelliseq.pl:8086/sample/$SAMPLEID/file/upload" -F file=@path_to_file
```
To set up another name for a file use `name` flag, to set up a tag - `tag` flag [tags may be multiple]
```bash
curl -H "Authorization: Bearer $TOKEN" "genetraps.intelliseq.pl:8086/sample/$SAMPLEID/file/upload" -F file=@path_to_file -F name=dir/new_filename -F tag=newtag -F tag=anothertag 
```

**To upload a file from url to aws (returns special aws path to the file):**
Special aws path can be recognized by '/' at the beginning (is followed by nr of sample, later optional dirs, ending with name of the file)
```bash
curl -H "Authorization: Bearer $TOKEN" "genetraps.intelliseq.pl:8086/sample/$SAMPLEID/url/upload" -d url='your_url'
```
To give another name, use `name` flag, to set up a tag - `tag` flag [tags may be multiple]  
```bash
curl -H "Authorization: Bearer $TOKEN" "genetraps.intelliseq.pl:8086/sample/$SAMPLEID/url/upload" -d url='your_url' -d name=dir/new_filename -F tag=newtag -F tag=anothertag 
```

**To list out content of a sample**
```bash
curl -H "Authorization: Bearer $TOKEN" "genetraps.intelliseq.pl:8086/sample/$SAMPLEID/ls"
```

**To delete a file in a sample**
```bash
curl -X DELETE -H "Authorization: Bearer $TOKEN" "genetraps.intelliseq.pl:8086/sample/$SAMPLEID/file/delete?path=aws_file_path"
```

**Activating the desired pipeline with desired (requested) outputs**
```bash
JOBID=$(curl -H "Authorization: Bearer $TOKEN" "genetraps.intelliseq.pl:8086/wdl" -H "accept: application/json" \
-d "workflowUrl=$WDL" \
-d "workflowInputs={\"minimal_workflow.minimal.str\":\"testing string\"}" \
-d "labels={\"sampleid\":\"$SAMPLEID\"}" \
-d "req-out={\"proper_name_of_desired_output_file\":\"new_name_for_the_output_file\"}" | jq -r ".id")
```
Flag `req-out`: json with pairs proper_name_of_desired_output_file:new_name_for_the_output_file [can be multiple])  
If you want to leave the original name of output file, then leave the value in key:value pair of req-out empty  
*If you don't specify any `req-out`, all output will be given.*

**Status of the workflow:**
```bash
curl -H "Authorization: Bearer $TOKEN" "genetraps.intelliseq.pl:8086/job/$JOBID/status" | jq -r ".status"
```


## Additional

**Job's properties:**
```bash
curl -H "Authorization: Bearer $TOKEN" "genetraps.intelliseq.pl:8086/job/$JOBID/status"
```

```bash
curl -H "Authorization: Bearer $TOKEN" "genetraps.intelliseq.pl:8086/job/$JOBID/abort"
```

**Job's outputs**
```bash
curl -H "Authorization: Bearer $TOKEN" "genetraps.intelliseq.pl:8086/job/$JOBID/output"
```

**Generate temporary links to download output**
Return json with pairs output_key:link  
With `sub` flag you can limit the resulting links to those which correspond with keys containing the specified substring
```bash
curl -H "Authorization: Bearer $TOKEN" "genetraps.intelliseq.pl:8086/job/$JOBID/output/download/links?sub="
```

**Logs
Generate temporary links to download**  
```bash
curl -H "Authorization: Bearer $TOKEN" "genetraps.intelliseq.pl:8086/job/$JOBID/logs/download/links"
```

**List all jobs associated with a sample**
```bash
curl -H "Authorization: Bearer $TOKEN" "genetraps.intelliseq.pl:8086/sample/$SAMPLEID/jobs"
```


## Example (with minimal.wdl)

```bash
TOKEN=$(curl -XPOST "web_app:secret@genetraps.intelliseq.pl:8088/oauth/token" -d grant_type=password -d client_id=web_app -d username=$STAGING_USERNAME -d password=$STAGING_PASSWORD | jq -r ".access_token")
```

```bash
SAMPLEID=$(curl -H "Authorization: Bearer $TOKEN" "genetraps.intelliseq.pl:8086/sample/create" | jq -r ".response")
```

```bash
WDL=https://gitlab.com/intelliseq/workflows/raw/dev/src/main/wdl/tasks/minimal-task/v1.0/minimal.wdl
```

```bash
JOBID=$(curl -H "Authorization: Bearer $TOKEN" "genetraps.intelliseq.pl:8086/wdl" -H "accept: application/json" \
-d "workflowUrl=$WDL" \
-d "workflowInputs={\"minimal_workflow.minimal.str\":\"testing string\"}" \
-d "labels={\"sampleid\":\"$SAMPLEID\"}" | jq -r ".id")
```

## Example 2 (with minimal v2.0)

```bash
TOKEN=$(curl -XPOST "web_app:secret@genetraps.intelliseq.pl:8088/oauth/token" -d grant_type=password -d client_id=web_app -d username=$STAGING_USERNAME -d password=$STAGING_PASSWORD | jq -r ".access_token")
```

```bash
SAMPLEID=$(curl -H "Authorization: Bearer $TOKEN" "genetraps.intelliseq.pl:8086/sample/create" | jq -r ".response")
```

```bash
WDL=https://gitlab.com/intelliseq/workflows/blob/dev/src/main/wdl/tasks/minimal-task/v2.0/minimal.wdl
```

```bash
FILEPATH=$(curl -H "Authorization: Bearer $TOKEN" "genetraps.intelliseq.pl:8086/sample/$SAMPLEID/file/upload" -F file=@minimal.txt | jq -r ".id")
```

```bash
JOBID=$(curl -H "Authorization: Bearer $TOKEN" "genetraps.intelliseq.pl:8086/wdl" -H "accept: application/json" \
-d "workflowUrl=$WDL" \
-d "workflowInputs={\"minimal_workflow.minimal_.str\":\"$FILEPATH\"}" \
-d "labels={\"sampleid\":\"$SAMPLEID\"}" | jq -r ".id")
```

## Example 3 (with minimal v3.0)

```bash
TOKEN=$(curl -XPOST "web_app:secret@genetraps.intelliseq.pl:8088/oauth/token" -d grant_type=password -d client_id=web_app -d username=$STAGING_USERNAME -d password=$STAGING_PASSWORD | jq -r ".access_token")
```

```bash
SAMPLEID=$(curl -H "Authorization: Bearer $TOKEN" "genetraps.intelliseq.pl:8086/sample/create" | jq -r ".response")
```

```bash
WDL=https://gitlab.com/intelliseq/workflows/blob/dev/src/main/wdl/tasks/minimal-task/v3.0/minimal.wdl
```

```bash
FILEPATH=$(curl -H "Authorization: Bearer $TOKEN" "genetraps.intelliseq.pl:8086/sample/$SAMPLEID/file/upload" -F file=@minimal.txt | jq -r ".id")
```

```bash
JSONINPUT="{\"minimal_workflow.minimal_file.strs\":[\"minimal\"], \"minimal_workflow.minimal_file.files\":[\"$FILEPATH\"]}"
```

```bash
JOBID=$(curl -H "Authorization: Bearer $TOKEN" "genetraps.intelliseq.pl:8086/wdl" -H "accept: application/json" \
-d "workflowUrl=$WDL" \
-d "workflowInputs=$JSONINPUT" \
-d "labels={\"sampleid\":\"$SAMPLEID\"}" | jq -r ".id")
```
