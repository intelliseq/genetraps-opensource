package pl.intelliseq.genetraps.api.dx.helpers;

import com.amazonaws.services.s3.AmazonS3;
import com.amazonaws.services.s3.AmazonS3ClientBuilder;
import com.amazonaws.services.s3.model.*;
import com.amazonaws.services.s3.transfer.TransferManager;
import com.amazonaws.services.s3.transfer.TransferManagerBuilder;
import com.amazonaws.services.s3.transfer.Upload;
import lombok.extern.log4j.Log4j2;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.core.env.Environment;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

@Log4j2
public class AWSApiProcessManager {

    @Autowired
    private Environment env;

    final private AmazonS3 s3Client = AmazonS3ClientBuilder.defaultClient();

    // one should be aware that in aws s3 there's no concept of folder-file system, it's more like a tree of keys
    // (but it can be simplified to thinking about it as about folders and files)
    public String runCreateSample(Integer sampleId) {

        String bucketName = env.getProperty("bucket-name");
        String sampleFolder = String.format("samples/%s/", sampleId);

        if (!s3Client.doesObjectExist(bucketName, sampleFolder)) {

            s3Client.copyObject(bucketName, "samples/0/", bucketName, sampleFolder);

//            return s3Client.getObject(bucketName, sampleId.toString()).getKey();
            return sampleId.toString();
        }
        else {
            // can be changed to returning empty string if a sample of given id already exists (but what with the case of failed creation?)
            return "A sample with given id already exists";
        }
    }



    public String runFileUpload(String path, Integer sampleId, String fileName, List<String> tags) throws InterruptedException {

        String bucketName = env.getProperty("bucket-name");
        String fileKey = String.format("samples/%s/%s", sampleId, fileName);

        if(s3Client.doesObjectExist(bucketName, fileKey)) {
            throw new InterruptedException("A file with given name (key) already exists");
        }

        // should rather be final global?
        TransferManager tm = TransferManagerBuilder.standard()
                .withS3Client(s3Client)
                .build();

        Upload upload = tm.upload(bucketName, fileKey, new File(path));

        try {
            upload.waitForCompletion();
//            System.out.println("File uploaded");
        } catch (InterruptedException e) {
            throw new InterruptedException("File failed to upload");
        }

        if(tags != null) {
//            GetObjectTaggingResult getTaggingResult = s3Client.getObjectTagging(new GetObjectTaggingRequest(bucketName, fileKey));
            List<Tag> fileTags = new ArrayList<>();
            for (String tag : tags) {
                fileTags.add(new Tag(tag, ""));
            }
            log.info(fileTags.toString());
            s3Client.setObjectTagging(new SetObjectTaggingRequest(bucketName, fileKey, new ObjectTagging(fileTags)));
        }

        return s3Client.getObject(bucketName, fileKey).getKey();
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
