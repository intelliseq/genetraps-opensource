package pl.intelliseq.genetraps.api.dx.helpers;

import com.amazonaws.services.s3.AmazonS3;
import com.amazonaws.services.s3.AmazonS3ClientBuilder;
import com.amazonaws.services.s3.model.ObjectTagging;
import com.amazonaws.services.s3.model.S3ObjectSummary;
import com.amazonaws.services.s3.model.SetObjectTaggingRequest;
import com.amazonaws.services.s3.model.Tag;
import com.amazonaws.services.s3.transfer.TransferManager;
import com.amazonaws.services.s3.transfer.TransferManagerBuilder;
import com.amazonaws.services.s3.transfer.Upload;
import com.mashape.unirest.http.HttpResponse;
import com.mashape.unirest.http.Unirest;
import com.mashape.unirest.http.exceptions.UnirestException;
import org.json.JSONObject;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.core.env.Environment;
import org.springframework.web.multipart.MultipartFile;

import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.util.ArrayList;
import java.util.List;

public class AWSApiProcessManager {

    @Autowired
    private Environment env;

    @Autowired
    private AuroraDBManager auroraDBManager;

    final private AmazonS3 s3Client = AmazonS3ClientBuilder.defaultClient();

    public String runCreateSample(Integer sampleId) {

        String bucketName = env.getProperty("bucket-name");
        String sampleFolder = String.format("samples/%s/", sampleId);

        if (!s3Client.doesObjectExist(bucketName, sampleFolder)) {

            s3Client.putObject(bucketName, sampleFolder, "");
//            s3Client.copyObject(bucketName, "samples/0/", bucketName, sampleFolder);

//            return s3Client.getObject(bucketName, sampleId.toString()).getKey();
            return sampleId.toString();
        }
        else {
            // can be changed to returning empty string if a sample of given id already exists (but what with the case of failed creation?)
            return "A sample with given id already exists";
        }
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
//            file = new File(mfile.getOriginalFilename());
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

        return s3Client.getObject(bucketName, fileKey).getKey();
    }

    public String runWdl(Integer userId, String workflowUrl, JSONObject workflowInputs, JSONObject labels) throws InterruptedException {


        // extracts a name of wdl from workflow url, even the ones like "wdl-name.1"
        String wdlName = workflowUrl.replaceFirst(".*/([\\w+\\.\\-]+)\\.wdl.*","$1");
        Integer wdlId = auroraDBManager.checkWdlId(wdlName);
        if(wdlId == null) {
            wdlId = auroraDBManager.getLastWdlId();
            if(wdlId == null)
                wdlId = 1;
            else
                wdlId += 1;
            auroraDBManager.putWdlToDB(wdlName);
        }

        Integer jobId = auroraDBManager.getLastJobId();
        if(jobId == null)
            jobId = 1;
        else
            jobId += 1;
        labels.put("jobid", jobId.toString());

        try {
            HttpResponse<String> response = Unirest.post(env.getProperty("ec2-cromwell-dns"))
                    .header("accept", "application/json")
                    .field("workflowUrl", workflowUrl)
                    .field("workflowInputs", workflowInputs)
                    .field("workflowType", new ByteArrayInputStream("WDL".getBytes()), "workflowtype")
                    .field("labels",labels)
                    .asString();
            String responseBody = response.getBody();
            if(response.getStatus() / 100 != 2)
                throw new InterruptedException(responseBody);
            // extracts the id of job like "7ab0-afd9" and similar
            String cromwellId = responseBody.replaceFirst(".*id\":\"([a-zA-Z0-9\\-]+)\".*", "$1");

            auroraDBManager.putJobToDB(userId, wdlId, cromwellId, labels.getInt("sampleid"));

            return responseBody;

        } catch (UnirestException e) {
            e.printStackTrace();
            throw new InterruptedException(e.getMessage());
        }
    }

    public List<String> runListObjects(String bucketName) {

        List<S3ObjectSummary> objectSummaries;
        if(bucketName == null)
             objectSummaries = s3Client.listObjectsV2(env.getProperty("bucket-name")).getObjectSummaries();
        else
            objectSummaries = s3Client.listObjectsV2(bucketName).getObjectSummaries();
        List<String> objects = new ArrayList<>();
        for (S3ObjectSummary object : objectSummaries) {
            objects.add("" + object.getKey());
        };

        return objects;
    }

    public List<String> runSampleLs(Integer sampleId) {

        List<S3ObjectSummary> objectSummaries = s3Client.listObjectsV2(env.getProperty("bucket-name"), String.format("samples/%s/", sampleId)).getObjectSummaries();
        List<String> objects = new ArrayList<>();
        for (S3ObjectSummary object : objectSummaries) {
            objects.add("" + object.getKey());
        };

        return objects;
    }
}
