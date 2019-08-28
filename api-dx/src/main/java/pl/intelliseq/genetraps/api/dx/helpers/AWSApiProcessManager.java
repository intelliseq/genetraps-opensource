package pl.intelliseq.genetraps.api.dx.helpers;

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
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

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

        try (BufferedInputStream in = new BufferedInputStream(new URL(url).openStream())) {

            TransferManager tm = TransferManagerBuilder.standard()
                    .withS3Client(s3Client)
                    .build();

            Upload upload;
            try {
                upload = tm.upload(bucketName, fileKey, in, new ObjectMetadata());
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

        return fileKey;
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
            upload = tm.upload(bucketName, fileKey, mfile.getInputStream(), new ObjectMetadata());
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

        return fileKey;
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
                            values.put(i, getAWSUrl(value.substring(1)));
                        }
                    }
                    workflowInputs.put(key, values);
                } else {
                    value = workflowInputs.getString(key);
                    if (value.isEmpty()) continue;
                    if (value.charAt(0) == '/') {
                        workflowInputs.put(key, getAWSUrl(value.substring(1)));
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
            String jobId = DigestUtils.md5Hex(toHash);
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

    public JSONObject runGetJobStatus(String jobId) throws InterruptedException {

        //TODO: change from mashape json to jackson json
        HttpResponse<com.mashape.unirest.http.JsonNode> response;
        try {
            response = Unirest
                    .get(String.format("%s/query", env.getProperty("cromwell.server")))
                    .header("accept", "application/json")
                    .queryString("label", String.format("jobId:%s", jobId))
                    .asJson();
        } catch (UnirestException e) {
            throw new InterruptedException(e.getMessage());
        }
        JSONObject responseBody = response.getBody().getObject().getJSONArray("results").getJSONObject(0);
        if (response.getStatus() / 100 != 2)
            throw new InterruptedException(responseBody.toString());

        return responseBody;
    }

    private String getAWSUrl(String fileName) throws InterruptedException {

        try {
            if(s3Client.doesObjectExist(env.getProperty("bucket.default"), fileName))
                return s3Client.getUrl(env.getProperty("bucket.default"), fileName).toString();
            else
                throw new InterruptedException(fileName);
        } catch (Exception e) {
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
        try {
            s3Client.deleteObject(bucket, fileRelPath);
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
