
# Documentation of using pipelines via genetraps API

**First and foremost, to use any of the following commands you need an authorization token:**
```bash
TOKEN=$(curl -XPOST "web_app:secret@genetraps.intelliseq.pl:8088/oauth/token" -d grant_type=password -d client_id=web_app -d username=$STAGING_USERNAME -d password=$STAGING_PASSWORD | jq -r ".access_token")
```

**You can generate the lowest available sample id (number) using the following command:**
```bash
SAMPLEID=$(curl -H "Authorization: Bearer $TOKEN" genetraps.intelliseq.pl:8086/sample/create | jq -r ".response")
```
or
**Using any given number:**
(at risk of using an existing one if unintentionally, although it won't delete the existing files in the sample, only mix in the new ones with them/overwrite some)
```bash
SAMPLEID=
```

**EC2 cromwell dns - obligatory:**
```bash
AWS_ADDRESS=http://ec2-100-25-40-253.compute-1.amazonaws.com/api/workflows/v1
```

**Link to your pipeline - obligatory:**
```bash
WDL=
```
(eg. https://gitlab.com/intelliseq/workflows/raw/dev/src/main/wdl/tasks/minimal-task/v1.0/minimal.wdl)

**To upload an input file from your machine to aws (returns relative path to the file):**
Relative path can be recognized by '/' at the beginning (is followed by nr of sample, later dirs, ending with name of file)
```bash
curl -H "Authorization: Bearer $TOKEN" genetraps.intelliseq.pl:8086/sample/$SAMPLEID/file/upload -F file=@relative_path_to_your_file
```
To set up another name for a file use `new-name` flag, to set up a tag - `tag` [tags can be multiple]
```bash
curl -H "Authorization: Bearer $TOKEN" genetraps.intelliseq.pl:8086/sample/$SAMPLEID/file/upload -F file=@relative_path_to_your_file -F new-name=dir/new_filename -F tag=newtag -F tag=othertag 
```

**Activating the desired pipeline with desired (requested) outputs (req-out [can be multiple]: json with key:value like proper_name_of_desired_output_file:new_name_for_the_output_file)**
```bash
JOBID=$(curl -H "Authorization: Bearer $TOKEN" genetraps.intelliseq.pl:8086/wdl -H "accept: application/json" \
-d "workflowUrl=$WDL" \
-d "workflowInputs={\"minimal_workflow.minimal.str\":\"testing string\"}" \
-d "labels={\"sampleid\":\"$SAMPLEID\"}" \
-d "req-out={\"proper_name_of_desired_output_file\":\"new_name_for_the_output_file\"}" | jq -r ".id")
```
(...-d "req-out={\"minimal_workflow.minimal.out\":\"new_name.out\"}" \...)
if you want to leave the original name of output file, then leave the value in key:value of req-out empty

**Status of the workflow:**
```bash
curl -X GET "$AWS_ADDRESS/query?label=jobId:$JOBID" -H "accept: application/json" | jq -r ".results" | jq -r ".[].status"
```


## Additional

**Job's properties:**
```bash
curl -X GET "$AWS_ADDRESS/query?label=jobId:$JOBID" -H "accept: application/json"
```

**Id of the workflow (one of jobs's properties):**
```bash
WORKFLOWID=$(curl -X GET "$AWS_ADDRESS/query?label=jobId:$JOBID" -H "accept: application/json" | jq -r ".results" | jq -r ".[].id")
```


## Example (with minimal.wdl)

```bash
TOKEN=$(curl -XPOST "web_app:secret@genetraps.intelliseq.pl:8088/oauth/token" -d grant_type=password -d client_id=web_app -d username=$STAGING_USERNAME -d password=$STAGING_PASSWORD | jq -r ".access_token")
```

```bash
SAMPLEID=$(curl -H "Authorization: Bearer $TOKEN" genetraps.intelliseq.pl:8086/sample/create | jq -r ".response")
```

```bash
WDL=https://gitlab.com/intelliseq/workflows/raw/dev/src/main/wdl/tasks/minimal-task/v1.0/minimal.wdl
```

```bash
JOBID=$(curl -H "Authorization: Bearer $TOKEN" genetraps.intelliseq.pl:8086/wdl -H "accept: application/json" \
-d "workflowUrl=$WDL" \
-d "workflowInputs={\"minimal_workflow.minimal.str\":\"testing string\"}" \
-d "labels={\"sampleid\":\"$SAMPLEID\"}" \
-d "req-out={\"minimal_workflow.minimal.out\":\"new_name.out\"}" | jq -r ".id")
```

## Example 2 (with minimal v2.0)

```bash
TOKEN=$(curl -XPOST "web_app:secret@genetraps.intelliseq.pl:8088/oauth/token" -d grant_type=password -d client_id=web_app -d username=$STAGING_USERNAME -d password=$STAGING_PASSWORD | jq -r ".access_token")
```

```bash
SAMPLEID=$(curl -H "Authorization: Bearer $TOKEN" genetraps.intelliseq.pl:8086/sample/create | jq -r ".response")
```

```bash
WDL=https://gitlab.com/intelliseq/workflows/blob/dev/src/main/wdl/tasks/minimal-task/v2.0/minimal.wdl
```

```bash
RELPATH=$(curl -H "Authorization: Bearer $TOKEN" genetraps.intelliseq.pl:8086/sample/$SAMPLEID/file/upload -F file=@minimal.txt | jq -r ".id")
```

```bash
JOBID=$(curl -H "Authorization: Bearer $TOKEN" genetraps.intelliseq.pl:8086/wdl -H "accept: application/json" \
-d "workflowUrl=$WDL" \
-d "workflowInputs={\"minimal_workflow.minimal.str\":\"$RELPATH\"}" \
-d "labels={\"sampleid\":\"$SAMPLEID\"}" \
-d "req-out={\"minimal_file.out\":\"\"}" | jq -r ".id")
```
