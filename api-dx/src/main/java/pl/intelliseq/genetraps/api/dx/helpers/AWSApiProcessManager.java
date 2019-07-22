package pl.intelliseq.genetraps.api.dx.helpers;

import com.amazonaws.services.s3.AmazonS3;
import com.amazonaws.services.s3.AmazonS3ClientBuilder;
import com.amazonaws.services.s3.model.*;
import com.amazonaws.services.s3.transfer.TransferManager;
import com.amazonaws.services.s3.transfer.TransferManagerBuilder;
import com.amazonaws.services.s3.transfer.Upload;
import com.mashape.unirest.http.HttpResponse;
import com.mashape.unirest.http.Unirest;
import com.mashape.unirest.http.exceptions.UnirestException;
import lombok.extern.log4j.Log4j2;
import org.apache.commons.codec.digest.DigestUtils;
import org.json.JSONObject;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.core.env.Environment;
import org.springframework.dao.DataIntegrityViolationException;
import org.springframework.web.multipart.MultipartFile;

import java.io.*;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

@Log4j2
public class AWSApiProcessManager {

    @Autowired
    private Environment env;

    @Autowired
    private AuroraDBManager auroraDBManager;

    @Autowired
    private FilesManager filesManager;

    final private AmazonS3 s3Client = AmazonS3ClientBuilder.defaultClient();
    final private TransferManager s3TransferManager = TransferManagerBuilder.standard().withS3Client(s3Client).build();         // ok? or maybe local

    public void runCreateSample(Integer sampleId) {

        String bucketName = env.getProperty("bucket-name");
        s3Client.putObject(bucketName, String.format("samples/%s/", sampleId), "");
    }

    public String runFileUpload(MultipartFile mfile, Integer sampleId, String newFileName, List<String> tags) throws InterruptedException {

        String bucketName = env.getProperty("bucket-name");
        String fileName = newFileName == null ? mfile.getOriginalFilename() : newFileName;
        String fileKey = String.format("samples/%s/%s", sampleId, fileName);

        if(s3Client.doesObjectExist(bucketName, fileKey)) {
            throw new InterruptedException("A file with given name (key) already exists");
        }

        // should rather be final global?
        TransferManager tm = TransferManagerBuilder.standard()
                .withS3Client(s3Client)
                .build();

        File file;
        try {
            file = new File(System.currentTimeMillis()+"fileupload");
            FileOutputStream fos = new FileOutputStream(file);
            fos.write(mfile.getBytes());
            fos.close();
        } catch (Exception e) {
            throw new InterruptedException("File failed to upload #1");
        }
        Upload upload = tm.upload(bucketName, fileKey, file);

        try {
            upload.waitForCompletion();
        } catch (InterruptedException e) {
            throw new InterruptedException("File failed to upload #2");
        }
        file.delete();

        if(tags != null) {
//            GetObjectTaggingResult getTaggingResult = s3Client.getObjectTagging(new GetObjectTaggingRequest(bucketName, fileKey));
            List<Tag> fileTags = new ArrayList<>();
            for (String tag : tags) {
                fileTags.add(new Tag(tag, ""));
            }
            s3Client.setObjectTagging(new SetObjectTaggingRequest(bucketName, fileKey, new ObjectTagging(fileTags)));
        }

        return String.format("/samples/%s/%s", sampleId, fileName);
    }

    public String runWdl(Integer userId, String workflowUrl, JSONObject workflowInputs, JSONObject labels, JSONObject relevantOutput) throws InterruptedException {

        // extracts a name of wdl from workflow url, even the ones like "wdl-name.1"
        String wdlName = workflowUrl;
        // gets id of wdl from db unless it's first time, otherwise adds new wdl id to db
        // TODO make wdl ids in db more universal
        Integer wdlId = auroraDBManager.checkWdlId(wdlName);
        if(wdlId == null) {
            try {
                auroraDBManager.putWdlToDB(wdlName);
                wdlId = auroraDBManager.checkWdlId(wdlName);
            } catch (DataIntegrityViolationException e) {
                wdlId = auroraDBManager.checkWdlId(wdlName);
            }
        }

        // creates correct url links if a string value in workflowInputs begins with "/"
        try {
            for (Iterator<String> it = workflowInputs.keys(); it.hasNext(); ) {
                String key = it.next();
                String value = workflowInputs.getString(key);
                if (value.isEmpty()) continue;
                if (value.charAt(0) == '/')
                    workflowInputs.put(key, getAWSUrl(value));
            }
        } catch (InterruptedException e) {
            throw new InterruptedException(String.format("%s - relative path to the file is not correct or the file is missing", e.getMessage()));
        }

        try {
            // creates unique 32-char id for a job, using hash
            String responseBody;
            String toHash = String.format("%s%d", auroraDBManager.getUserDetails(userId).getUserName(), System.currentTimeMillis());
            String jobId = DigestUtils.md5Hex(toHash);
            labels.put("jobId", jobId);

            HttpResponse<String> response = Unirest.post(env.getProperty("ec2-cromwell-dns"))
                    .header("accept", "application/json")
                    .field("workflowUrl", workflowUrl)
                    .field("workflowInputs", workflowInputs)
                    .field("workflowType", new ByteArrayInputStream("WDL".getBytes()), "workflowtype")
                    .field("labels", labels)
                    .asString();
            responseBody = response.getBody();
            if (response.getStatus() / 100 != 2)
                throw new InterruptedException(responseBody);

            auroraDBManager.putJobToDB(jobId, userId, wdlId, 0, labels.getInt("sampleid"), relevantOutput);
            return jobId;

        } catch (UnirestException e) {
            throw new InterruptedException(e.getMessage());
        } catch (Exception e) {
            throw new InterruptedException(e.getMessage());
        }
    }

    public String getAWSUrl(String fileName) throws InterruptedException {

        try {
            if(s3Client.doesObjectExist(env.getProperty("bucket-name"), fileName.substring(1)))
                return s3Client.getUrl(env.getProperty("bucket-name"), fileName).toString();
            else
                throw new InterruptedException(fileName);
        } catch (Exception e) {
            throw new InterruptedException(e.getMessage());
        }
    }

    public JSONObject runSampleLs(Integer sampleId, String dir) throws InterruptedException {

        String bucket = env.getProperty("bucket-name");
        try {
            if(!s3Client.doesObjectExist(bucket, String.format("samples/%s/", sampleId)))
                throw new InterruptedException(String.format("Sample: %s does not exist", sampleId));
//            else if(!dir.isEmpty() && !s3Client.doesObjectExist(bucket, String.format("samples/%s/%s/", sampleId, dir)))
//                throw new InterruptedException(String.format("Directory: %s does not exist in sample: %s", dir, sampleId));

            List<S3ObjectSummary> objectSummaries = s3Client.listObjectsV2(bucket,
                    dir.isEmpty() ? String.format("samples/%s/", sampleId) : String.format("samples/%s/%s/", sampleId, dir)).getObjectSummaries();
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
}
