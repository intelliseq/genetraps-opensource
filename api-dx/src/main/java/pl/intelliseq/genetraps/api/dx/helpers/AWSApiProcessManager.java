package pl.intelliseq.genetraps.api.dx.helpers;

import com.amazonaws.HttpMethod;
import com.amazonaws.services.dynamodbv2.xspec.S;
import com.amazonaws.services.s3.AmazonS3;
import com.amazonaws.services.s3.AmazonS3ClientBuilder;
import com.amazonaws.services.s3.model.*;
import com.amazonaws.services.s3.transfer.TransferManager;
import com.amazonaws.services.s3.transfer.TransferManagerBuilder;
import com.amazonaws.services.s3.transfer.Upload;
import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.node.ObjectNode;
import com.mashape.unirest.http.HttpResponse;
import com.mashape.unirest.http.Unirest;
import com.mashape.unirest.http.exceptions.UnirestException;
import lombok.extern.log4j.Log4j2;
import org.apache.commons.codec.digest.DigestUtils;
import org.apache.commons.io.IOUtils;
import org.json.JSONArray;
import org.json.JSONObject;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.core.env.Environment;
import org.springframework.dao.DataIntegrityViolationException;
import org.springframework.web.multipart.MultipartFile;
import pl.intelliseq.genetraps.api.dx.exceptions.PropertiesException;

import java.io.BufferedInputStream;
import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.nio.charset.StandardCharsets;
import java.util.*;

@Log4j2
public class AWSApiProcessManager {

    @Autowired
    private Environment env;

    @Autowired
    private AuroraDBManager auroraDBManager;

    @Autowired
    private FilesManager filesManager;

    //TODO: probably cant be final global because diffrent users
    final private AmazonS3 s3Client = AmazonS3ClientBuilder.defaultClient();
//    final private TransferManager s3TransferManager = TransferManagerBuilder.standard().withS3Client(s3Client).build();         // ok? or maybe local

    public void runCreateSample(Integer sampleId) {

        String bucketName = env.getProperty("bucket.default");
        s3Client.putObject(bucketName, String.format("%s/%s/", env.getProperty("samples.folder"), sampleId), "");
    }

    public Integer runCreateSampl(Integer sampleId) {

        String bucketName = env.getProperty("bucket.default");
        s3Client.putObject(bucketName, String.format("%s/%s/", env.getProperty("samples.folder"), sampleId), "");
        return sampleId;
    }

    public String runUrlFetch(String url, Integer sampleId, String newFileName, List<String> tags) throws InterruptedException {

        String bucketName = env.getProperty("bucket.default");
        String fileName = newFileName == null ? url.substring(url.lastIndexOf("/") + 1) : newFileName.charAt(newFileName.length() - 1) == '/' ? String.format("%s%s", newFileName, url.substring(url.lastIndexOf("/") + 1)) : newFileName;
        log.info(url.substring(url.lastIndexOf("/") + 1));
        String fileKey = String.format("%s/%s/%s", env.getProperty("samples.folder"), sampleId, fileName);

        if(s3Client.doesObjectExist(bucketName, fileKey)) {
            throw new InterruptedException("A file with given name (key) already exists");
        }

        long streamContentLength = 0;
        try (BufferedInputStream in = new BufferedInputStream(new URL(url).openStream())) {

            byte[] buffer = new byte[8 * 1024];
            int bytesRead;
            while((bytesRead = in.read(buffer)) != -1) {
                streamContentLength += bytesRead;
            }
            log.info(streamContentLength);
        } catch (IOException e) {
            throw new InterruptedException(e.getMessage());
        }

        try (BufferedInputStream in = new BufferedInputStream(new URL(url).openStream())) {

            TransferManager tm = TransferManagerBuilder.standard()
                    .withS3Client(s3Client)
                    .build();

            Upload upload;
            try {
                ObjectMetadata fileMetadata = new ObjectMetadata();
                fileMetadata.setContentLength(streamContentLength);
                upload = tm.upload(bucketName, fileKey, in, fileMetadata);
            } catch (Exception e) {
                throw new InterruptedException(String.format("File: %s failed to upload #1: %s", fileKey, e.getMessage()));
            }

            try {
                upload.waitForCompletion();
            } catch (InterruptedException e) {
                throw new InterruptedException(String.format("File: %s failed to upload #2: %s", fileKey, e.getMessage()));
            }

        } catch (IOException e) {
            throw new InterruptedException(e.getMessage());
        }

        if(tags != null) {
//            GetObjectTaggingResult getTaggingResult = s3Client.getObjectTagging(new GetObjectTaggingRequest(bucketName, fileKey));
            List<Tag> fileTags = new ArrayList<>();
            for (String tag : tags) {
                fileTags.add(new Tag(tag, ""));
            }
            s3Client.setObjectTagging(new SetObjectTaggingRequest(bucketName, fileKey, new ObjectTagging(fileTags)));
        }

        return String.format("/%s", fileKey);
    }

    public String runFileUpload(MultipartFile mfile, Integer sampleId, String newFileName, List<String> tags) throws InterruptedException {

        String bucketName = env.getProperty("bucket.default");
        String fileName = newFileName == null ? mfile.getOriginalFilename() : newFileName.charAt(newFileName.length() - 1) == '/' ? String.format("%s%s", newFileName, mfile.getOriginalFilename()) : newFileName;
        String fileKey = String.format("%s/%s/%s", env.getProperty("samples.folder"), sampleId, fileName);

        if(s3Client.doesObjectExist(bucketName, fileKey)) {
            throw new InterruptedException("A file with given name (key) already exists");
        }

        TransferManager tm = TransferManagerBuilder.standard()
                .withS3Client(s3Client)
                .build();

        Upload upload;
        try {
            ObjectMetadata fileMetadata = new ObjectMetadata();
            fileMetadata.setContentLength(mfile.getSize());
            upload = tm.upload(bucketName, fileKey, mfile.getInputStream(), fileMetadata);
        } catch (Exception e) {
            throw new InterruptedException(String.format("File: %s failed to upload #1: %s", fileKey, e.getMessage()));
        }

        try {
            upload.waitForCompletion();
        } catch (InterruptedException e) {
            throw new InterruptedException(String.format("File: %s failed to upload #2: %s", fileKey, e.getMessage()));
        }

        if(tags != null) {
//            GetObjectTaggingResult getTaggingResult = s3Client.getObjectTagging(new GetObjectTaggingRequest(bucketName, fileKey));
            List<Tag> fileTags = new ArrayList<>();
            for (String tag : tags) {
                fileTags.add(new Tag(tag, ""));
            }
            s3Client.setObjectTagging(new SetObjectTaggingRequest(bucketName, fileKey, new ObjectTagging(fileTags)));
        }

        return String.format("/%s", fileKey);
    }

    public String runWdl(Integer userId, String workflowUrl, JSONObject workflowInputs, JSONObject labels, JSONObject requestedOutputs) throws InterruptedException {

        // extracts a name of wdl from workflow url, even the ones like "wdl-name.1"
        // gets id of wdl from db unless it's first time, otherwise adds new wdl id to db
        // TODO make wdl ids in db more universal
        Integer wdlId = auroraDBManager.checkWdlId(workflowUrl);
        if(wdlId == null) {
            try {
                auroraDBManager.putWdlToDB(workflowUrl);
                wdlId = auroraDBManager.checkWdlId(workflowUrl);
            } catch (DataIntegrityViolationException e) {
                wdlId = auroraDBManager.checkWdlId(workflowUrl);
            }
        }
        log.info(wdlId);

        // creates correct url links if a string value in workflowInputs begins with "/"
        try {
            for (Iterator<String> it = workflowInputs.keys(); it.hasNext(); ) {
                String key = it.next();
                String value;
                if (workflowInputs.get(key).getClass().equals(JSONArray.class)) {
                    JSONArray values = workflowInputs.getJSONArray(key);
                    for (int i = 0; i < values.length(); i++) {
                        value = values.getString(i);
                        if (value.isEmpty()) continue;
                        if (value.charAt(0) == '/') {
                            values.put(i, getAWSFullPath(value.substring(1)));
                        }
                    }
                    workflowInputs.put(key, values);
                } else {
                    value = workflowInputs.getString(key);
                    if (value.isEmpty()) continue;
                    if (value.charAt(0) == '/') {
                        workflowInputs.put(key, getAWSFullPath(value.substring(1)));
                    }
                }
            }
        } catch (InterruptedException e) {
            throw new InterruptedException(String.format("%s - relative path to the file is not correct or the file is missing", e.getMessage()));
        }

        log.info(workflowInputs.toString());

        try {
            // creates unique 32-char id for a job, using hash
            String responseBody;
            String toHash = String.format("%s%d", auroraDBManager.getUserDetails(userId).getUserName(), System.currentTimeMillis());
            // used to get job status and date
            String jobId = DigestUtils.md5Hex(toHash);
            log.info(jobId);
            labels.put("jobId", jobId);

            HttpResponse<String> response = Unirest.post(env.getProperty("cromwell.server"))
                    .header("accept", "application/json")
                    .field("workflowUrl", workflowUrl)
                    .field("workflowInputs", workflowInputs)
                    .field("workflowType", new ByteArrayInputStream("WDL".getBytes()), "workflowtype")
                    .field("labels", labels)
                    .asString();
            responseBody = response.getBody();
            if (response.getStatus() / 100 != 2)
                throw new InterruptedException(responseBody);

            auroraDBManager.putJobToDB(jobId, userId, wdlId, 0, labels.getInt("sampleid"), requestedOutputs);
            return jobId;

        } catch (UnirestException e) {
            throw new InterruptedException(e.getMessage());
        } catch (Exception e) {
            throw new InterruptedException(e.getMessage());
        }
    }

    private String getAWSFullPath(String fileName) throws InterruptedException {

        try {
            if(s3Client.doesObjectExist(env.getProperty("bucket.default"), fileName))
//                return s3Client.getUrl(env.getProperty("bucket.default"), fileName).toString();
                return String.format("s3://%s/%s", env.getProperty("bucket.default"), fileName);
            else
                throw new InterruptedException(fileName);
        } catch (Exception e) {
            throw new InterruptedException(e.getMessage());
        }
    }

    public JsonNode runGetJobStatus(String jobId) throws InterruptedException {

        JsonNode response = getJobWithLabel("jobId", jobId);
        if(response.get("totalResultsCount").asInt() == 0)
            return new ObjectMapper().createObjectNode().put("id", "No job was found");
        ObjectNode jobsExtractedFromBigJson = (ObjectNode) response.get("results").get(0);
        return jobsExtractedFromBigJson.put("id", jobId);
    }

    public JsonNode runAbortJob(String jobId) throws InterruptedException {

        ObjectNode response = (ObjectNode) getJobWithLabel("jobId", jobId);
        if(response.get("totalResultsCount").asInt() == 0)
            return new ObjectMapper().createObjectNode().put("id", "No job was found");
        String jobCromwellId = response.get("results").get(0).get("id").asText();

        try {
            HttpResponse<com.mashape.unirest.http.JsonNode> responseUnirest = Unirest
                    .post(String.format("%s/%s/abort", env.getProperty("cromwell.server"), jobCromwellId))
                    .header("accept", "application/json")
                    .asJson();
            response = new ObjectMapper().readValue(responseUnirest.getRawBody(), ObjectNode.class);
            if (responseUnirest.getStatus() / 100 != 2)
                throw new InterruptedException(response.toString());
            if(response.get("status").asText().equals("error") || response.get("status").asText().equals("fail"))
                throw new InterruptedException(response.get("message").asText());
        } catch (UnirestException | IOException e) {
            throw new InterruptedException(e.getMessage());
        }

        return response.put("id", jobId);
    }

    public JsonNode runGetJobOutputs(String jobId) throws InterruptedException {

        ObjectNode response = (ObjectNode) getJobWithLabel("jobId", jobId);
        if(response.get("totalResultsCount").asInt() == 0)
            return new ObjectMapper().createObjectNode().put("id", "No job was found");
        String jobCromwellId = response.get("results").get(0).get("id").asText();

        try {
            HttpResponse<com.mashape.unirest.http.JsonNode> responseUnirest = Unirest
                    .get(String.format("%s/%s/outputs", env.getProperty("cromwell.server"), jobCromwellId))
                    .header("accept", "application/json")
                    .asJson();
            response = new ObjectMapper().readValue(responseUnirest.getRawBody(), ObjectNode.class);
            if (responseUnirest.getStatus() / 100 != 2)
                throw new InterruptedException(response.toString());
        } catch (UnirestException | IOException e) {
            throw new InterruptedException(e.getMessage());
        }

        response = (ObjectNode) response.get("outputs");
        for (Iterator<Map.Entry<String, JsonNode>> it = response.fields(); it.hasNext(); ) {
            String outputKey = it.next().getKey();
            String outputLink = response.get(outputKey).asText();
            response.put(outputKey, outputLink.substring(outputLink.lastIndexOf('/') + 1));
        }

        return response;
    }

    public JsonNode runGetJobOutputsDownloadLinks(String jobId, String sub) throws InterruptedException {

        ObjectNode response = (ObjectNode) getJobWithLabel("jobId", jobId);
        if(response.get("totalResultsCount").asInt() == 0)
            return new ObjectMapper().createObjectNode().put("id", "No job was found");
        String jobCromwellId = response.get("results").get(0).get("id").asText();

        try {
            HttpResponse<com.mashape.unirest.http.JsonNode> responseUnirest = Unirest
                    .get(String.format("%s/%s/outputs", env.getProperty("cromwell.server"), jobCromwellId))
                    .header("accept", "application/json")
                    .asJson();
            response = new ObjectMapper().readValue(responseUnirest.getRawBody(), ObjectNode.class);
            if (responseUnirest.getStatus() / 100 != 2)
                throw new InterruptedException(response.toString());
        } catch (UnirestException | IOException e) {
            throw new InterruptedException(e.getMessage());
        }

        String bucketName = env.getProperty("bucket.default");

        java.util.Date expiration = new java.util.Date();
        // an hour
        long expTimeMillis = 1000 * 60 * 60;

        response = (ObjectNode) response.get("outputs");
        LinkedList<String> outputsNotWithSubstring = new LinkedList<>();
        for (Iterator<Map.Entry<String, JsonNode>> it = response.fields(); it.hasNext(); ) {
            String outputKey = it.next().getKey();
            if(!outputKey.contains(sub)) {
                outputsNotWithSubstring.add(outputKey);
                continue;
            }
            String outputLink = response.get(outputKey).asText();
            outputLink = outputLink.substring(outputLink.indexOf(bucketName) + bucketName.length() + 1);
            expiration.setTime(expiration.getTime() + expTimeMillis);
            GeneratePresignedUrlRequest generatePresignedUrlRequest =
                    new GeneratePresignedUrlRequest(bucketName, outputLink)
                            .withMethod(HttpMethod.GET)
                            .withExpiration(expiration);
            response.put(outputKey, s3Client.generatePresignedUrl(generatePresignedUrlRequest).toString().replaceAll("%2F", "/"));
        }
        if(outputsNotWithSubstring.isEmpty())
            return response;
        return response.without(outputsNotWithSubstring);
    }

    public JsonNode runGetJobLogs(String jobId) throws InterruptedException {

        JsonNode response = getJobWithLabel("jobId", jobId);
        if(response.get("totalResultsCount").asInt() == 0)
            return new ObjectMapper().createObjectNode().put("id", "No job was found");
        String jobCromwellId = response.get("results").get(0).get("id").asText();

        try {
            HttpResponse<com.mashape.unirest.http.JsonNode> responseUnirest = Unirest
                    .get(String.format("%s/%s/metadata", env.getProperty("cromwell.server"), jobCromwellId))
                    .header("accept", "application/json")
                    .asJson();
            response = new ObjectMapper().readValue(responseUnirest.getRawBody(), JsonNode.class).get("calls");
            if (responseUnirest.getStatus() / 100 != 2)
                throw new InterruptedException(response.toString());
        } catch (UnirestException | IOException e) {
            throw new InterruptedException(e.getMessage());
        }

        String bucketName = env.getProperty("bucket.default");

        List<String> unnecessaryKeys = Arrays.asList("attempt", "shardIndex");
        ObjectNode failedWorkflowsLogsLinks = new ObjectMapper().createObjectNode();
        // an hour
        long expTimeMillis = 1000 * 60 * 60;

        for (Iterator<Map.Entry<String, JsonNode>> workflowsMetadataIt = response.fields(); workflowsMetadataIt.hasNext(); ) {
            String workflowName = workflowsMetadataIt.next().getKey();
            JsonNode workflowMetadata = response.get(workflowName).get(0);
            if(!workflowMetadata.has("executionStatus") || !workflowMetadata.get("executionStatus").asText().equals("Failed"))
                continue;
            String workflowCromwellId;
            if(workflowMetadata.has("jobId"))
                workflowCromwellId = workflowMetadata.get("jobId").asText();
            else if(workflowMetadata.has("subWorkflowId"))
                workflowCromwellId = workflowMetadata.get("subWorkflowId").asText();
            else
                continue;

            // json with jsons of logs for a bunch of related workflows
            ObjectNode responseLogs;
            try {
                HttpResponse<com.mashape.unirest.http.JsonNode> responseUnirest = Unirest
                        .get(String.format("%s/%s/logs", env.getProperty("cromwell.server"), workflowCromwellId))
                        .header("accept", "application/json")
                        .asJson();
                responseLogs = (ObjectNode) new ObjectMapper().readValue(responseUnirest.getRawBody(), JsonNode.class).get("calls");
                if (responseUnirest.getStatus() / 100 != 2)
                    throw new InterruptedException(responseLogs.toString());
            } catch (UnirestException | IOException e) {
                throw new InterruptedException(e.getMessage());
            }

            for (Iterator<Map.Entry<String, JsonNode>> it = responseLogs.fields(); it.hasNext(); ) {
                // json with logs for one of the workflows
                String workflowLogsName = it.next().getKey();
                ObjectNode workflowLogs = (ObjectNode) responseLogs.get(workflowLogsName).get(0);
                if(workflowLogs.has("stderr")) {
                    String logLink = workflowLogs.get("stderr").asText();
                    logLink = getAWSKeyFromS3BucketKeyFormat(bucketName, logLink);
                    workflowLogs.put("stderr", getAWSPresignedUrl(bucketName, logLink, expTimeMillis));
                }
                if(workflowLogs.has("stdout")) {
                    String logLink = workflowLogs.get("stdout").asText();
                    logLink = getAWSKeyFromS3BucketKeyFormat(bucketName, logLink);
                    workflowLogs.put("stdout", getAWSPresignedUrl(bucketName, logLink, expTimeMillis));
                }
                responseLogs.set(workflowLogsName, workflowLogs.remove(unnecessaryKeys));
//                log.info("--      " + workflowLogs);
            }
//            log.info(responseLogs.toString());
            failedWorkflowsLogsLinks.set(workflowName, responseLogs);
        }
//        log.info("Final Result:     " + failedWorkflowsLogsLinks.toString());
        return failedWorkflowsLogsLinks.put("id", jobId);
    }

    private String getAWSPresignedUrl(String bucketName, String objectLink, long expTimeMillis) {

        java.util.Date expiration = new java.util.Date();
        expiration.setTime(expiration.getTime() + expTimeMillis);
        GeneratePresignedUrlRequest generatePresignedUrlRequest =
                new GeneratePresignedUrlRequest(bucketName, objectLink)
                        .withMethod(HttpMethod.GET)
                        .withExpiration(expiration);
        return s3Client.generatePresignedUrl(generatePresignedUrlRequest).toString().replaceAll("%2F", "/");
    }

    private String getAWSKeyFromS3BucketKeyFormat(String bucketName, String s3bucketkeyformatkey) {

        return s3bucketkeyformatkey.substring(s3bucketkeyformatkey.indexOf(bucketName) + bucketName.length() + 1);
    }

    public JsonNode runGetJobsForSample(String sampleId) throws InterruptedException {

        ObjectNode jobs = new ObjectMapper().createObjectNode();
        JsonNode response;

        List<String> jobsFromDB = auroraDBManager.getJobsWithSampleID(sampleId);
        int jobCount = 0;
        for(String jobId : jobsFromDB) {

            log.info(jobId);
            response = getJobWithLabel("jobId", jobId);
            if(response.get("totalResultsCount").asInt() == 0)
                continue;
            jobCount += response.get("totalResultsCount").asInt();
            jobs.set(jobId, ((ObjectNode) getJobWithLabel("jobId", jobId).get("results").get(0)).without("id"));
        }
        jobs.put("count", jobCount);

        return jobs;
    }

    private JsonNode getJobWithLabel(String label, String value) throws InterruptedException {

        try {
            HttpResponse<com.mashape.unirest.http.JsonNode> responseUnirest = Unirest
                    .get(String.format("%s/query", env.getProperty("cromwell.server")))
                    .header("accept", "application/json")
                    .queryString("label", String.format("%s:%s", label, value))
                    .asJson();
            JsonNode response = new ObjectMapper().readValue(responseUnirest.getRawBody(), ObjectNode.class);
            if (responseUnirest.getStatus() / 100 != 2)
                throw new InterruptedException(response.toString());

            return response;
        } catch (UnirestException | IOException e) {
            throw new InterruptedException(e.getMessage());
        }
    }

    public JSONObject runSampleLs(Integer sampleId, String dir) throws InterruptedException {

        String bucket = env.getProperty("bucket.default");
        try {
            if(!s3Client.doesObjectExist(bucket, String.format("%s/%s/", env.getProperty("samples.folder"), sampleId)))
                throw new InterruptedException(String.format("Sample: %s does not exist", sampleId));
//            else if(!dir.isEmpty() && !s3Client.doesObjectExist(bucket, String.format("samples/%s/%s/", sampleId, dir)))
//                throw new InterruptedException(String.format("Directory: %s does not exist in sample: %s", dir, sampleId));

            List<S3ObjectSummary> objectSummaries = s3Client.listObjectsV2(bucket,
                    dir.isEmpty() ? String.format("%s/%s/", env.getProperty("samples.folder"), sampleId) : String.format("%s/%s/%s/", env.getProperty("samples.folder"), sampleId, dir)).getObjectSummaries();
//            log.info(dir.isEmpty() ? String.format("samples/%s/", sampleId) : String.format("samples/%s/%s/", sampleId, dir));
            JSONObject objects = new JSONObject();
            String objectKey;
            for (S3ObjectSummary object : objectSummaries) {
                objectKey = object.getKey();
                if(objectKey.matches(".*/"))    continue;       // omits id of directories
                objects.put(objectKey.substring(objectKey.lastIndexOf("/") + 1), String.format("/%s", objectKey));
            }

            return objects;
        } catch (Exception e) {
            throw new InterruptedException(e.getMessage());
        }
    }

    public String runDeleteFile(Integer sampleId, String fileRelPath) throws InterruptedException {

        String bucket = env.getProperty("bucket.default");
        if(!fileRelPath.matches(String.format("/samples/%s/.*", sampleId)))
            throw new InterruptedException("There is a problem with your path or sample id");
        try {
            s3Client.deleteObject(bucket, fileRelPath.substring(1));
        } catch (Exception e) {
            throw new InterruptedException("Deletion of file failed: " + e.getMessage());
        }

        return String.format("Deleted: %s", fileRelPath);
    }

    // creates and/or adds new properties to a specified sample folder
    public JsonNode propertiesPost(Integer sampleId, JsonNode propertiesToPost, boolean bePersistent) throws PropertiesException {

        String bucket = env.getProperty("bucket.default");
        StringBuilder builder = new StringBuilder();
        try {
            for (Iterator<Map.Entry<String, JsonNode>> it = propertiesToPost.fields(); it.hasNext(); ) {
                String propertyKey = it.next().getKey();
                if (!s3Client.doesObjectExist(bucket, String.format("%s/%s/%s/%s", env.getProperty("samples.folder"), sampleId, env.getProperty("properties.folder"), propertyKey)))
                    log.info(String.format("%s/%s/%s/%s", env.getProperty("samples.folder"), sampleId, env.getProperty("properties.folder"), propertyKey) + " " + propertiesToPost.get(propertyKey).toString());
//                    s3Client.putObject(bucket, String.format("%s/%s/%s/%s", env.getProperty("samples.folder"), sampleId, env.getProperty("properties.folder"), propertyKey), propertiesToPost.get(propertyKey).toString());
                else if(bePersistent)
                    builder.append(System.lineSeparator()).append(propertyKey);
                else
                    throw new PropertiesException(propertyKey);
            }
            if(builder.length() != 0) {
                throw new PropertiesException(builder.toString());
            }
        } catch (PropertiesException e) {
            if(bePersistent)
                throw new PropertiesException(String.format("Adding new properties finished.%sFollowing properties couldn't be added: %s", System.lineSeparator(), e.getMessage()));
            throw new PropertiesException(String.format("Adding properties interrupted.%sThe following property couldn't be added: %s", System.lineSeparator(), e.getMessage()));
        } catch (Exception e) {
            throw new PropertiesException(e.getMessage());
        }

        return propertiesToPost;
    }

    // returns properties of a specified sample folder
    public JsonNode propertiesGet(Integer sampleId) throws PropertiesException {

        String bucket = env.getProperty("bucket.default");
        try {
            List<S3ObjectSummary> propertySummaries = s3Client.listObjectsV2(bucket, String.format("%s/%s/%s/", env.getProperty("samples.folder"), sampleId, env.getProperty("properties.folder"))).getObjectSummaries();

            if (propertySummaries.isEmpty()) {
                throw new PropertiesException("does not exist");
            }
            // dummy variable to tell if the folder exists -> throws an exception if it doesn't
            ObjectNode properties = new ObjectMapper().createObjectNode();
            String propertyKey;
            for(S3ObjectSummary property : propertySummaries) {
                propertyKey = property.getKey();
                if(propertyKey.matches(".*/"))  continue;
                InputStream propertyData = s3Client.getObject(bucket, propertyKey).getObjectContent();
                properties.put(propertyKey, IOUtils.toString(propertyData, StandardCharsets.UTF_8));
            }
            return properties;
        } catch (Exception e) {
            throw new PropertiesException(String.format("cannot exist, because sample: %d does not exist", sampleId));
        }
    }

    // changes/updates properties of a specified sample folder
    public JsonNode propertiesPut(Integer sampleId, JsonNode propertiesToUpdate, boolean bePersistent) throws PropertiesException {

        String bucket = env.getProperty("bucket.default");
        StringBuilder builder = new StringBuilder();
        try {
            for (Iterator<Map.Entry<String, JsonNode>> it = propertiesToUpdate.fields(); it.hasNext(); ) {
                String propertyKey = it.next().getKey();
                if (s3Client.doesObjectExist(bucket, String.format("%s/%s/%s/%s", env.getProperty("samples.folder"), sampleId, env.getProperty("properties.folder"), propertyKey)))
                    s3Client.putObject(bucket, String.format("%s/%s/%s/%s", env.getProperty("samples.folder"), sampleId, env.getProperty("properties.folder"), propertyKey), propertiesToUpdate.get(propertyKey).toString());
                else if(bePersistent)
                    builder.append(System.lineSeparator()).append(propertyKey);
                else
                    throw new PropertiesException(propertyKey);
            }
            if(builder.length() != 0) {
                throw new PropertiesException(builder.toString());
            }
        } catch (PropertiesException e) {
            if(bePersistent)
                throw new PropertiesException(String.format("Updating properties finished.%sFollowing properties couldn't be updated: %s", System.lineSeparator(), e.getMessage()));
            throw new PropertiesException(String.format("Updating properties interrupted.%sThe following property couldn't be updated: %s", System.lineSeparator(), e.getMessage()));
        } catch (Exception e) {
            throw new PropertiesException(e.getMessage());
        }

        return propertiesToUpdate;
    }

    // deletes properties of a specified sample folder
    public JsonNode propertiesDelete(Integer sampleId, JsonNode propertiesToDelete, boolean bePersistent) throws PropertiesException {

        String bucket = env.getProperty("bucket.default");
        StringBuilder builder = new StringBuilder();
        try {
            for (Iterator<Map.Entry<String, JsonNode>> it = propertiesToDelete.fields(); it.hasNext(); ) {
                String propertyKey = it.next().getKey();
                if (s3Client.doesObjectExist(bucket, String.format("%s/%s/%s/%s", env.getProperty("samples.folder"), sampleId, env.getProperty("properties.folder"), propertyKey)))
                    log.info(String.format("%s/%s/%s/%s", env.getProperty("samples.folder"), sampleId, env.getProperty("properties.folder"), propertyKey));
//                    s3Client.deleteObject(bucket, String.format("%s/%s/%s/%s", env.getProperty("samples.folder"), sampleId, env.getProperty("properties.folder"), propertyKey));
                else if(bePersistent)
                    builder.append(System.lineSeparator()).append(propertyKey);
                else
                    throw new PropertiesException(propertyKey);
            }
            if(builder.length() != 0) {
                throw new PropertiesException(builder.toString());
            }
        } catch (PropertiesException e) {
            if(bePersistent)
                throw new PropertiesException(String.format("Deleting properties finished.%sFollowing properties couldn't be deleted: %s", System.lineSeparator(), e.getMessage()));
            throw new PropertiesException(String.format("Deleting properties interrupted.%sThe following property couldn't be deleted: %s", System.lineSeparator(), e.getMessage()));
        } catch (Exception e) {
            throw new PropertiesException(e.getMessage());
        }

        return propertiesToDelete;
    }
}
