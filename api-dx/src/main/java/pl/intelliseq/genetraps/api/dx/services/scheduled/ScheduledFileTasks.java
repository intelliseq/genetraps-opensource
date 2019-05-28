package pl.intelliseq.genetraps.api.dx.services.scheduled;

import com.amazonaws.services.s3.AmazonS3;
import com.amazonaws.services.s3.AmazonS3ClientBuilder;
import com.amazonaws.services.s3.model.DeleteObjectsRequest;
import com.amazonaws.services.s3.model.S3ObjectSummary;
import com.mashape.unirest.http.HttpResponse;
import com.mashape.unirest.http.JsonNode;
import com.mashape.unirest.http.Unirest;
import lombok.extern.log4j.Log4j2;
import org.json.JSONObject;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.core.env.Environment;
import org.springframework.jdbc.support.rowset.SqlRowSet;
import pl.intelliseq.genetraps.api.dx.helpers.AuroraDBManager;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

@Log4j2
public class ScheduledFileTasks {

    @Autowired
    private Environment env;

    @Autowired
    private AuroraDBManager auroraDBManager;

    final private AmazonS3 s3Client = AmazonS3ClientBuilder.defaultClient();

    public void deleteAllObjectsWithPrefix(String prefix) {
        List<S3ObjectSummary> objectSummaries = s3Client.listObjectsV2(env.getProperty("bucket-name"), prefix).getObjectSummaries();
        List<String> objects = new ArrayList<>();
        for (S3ObjectSummary object : objectSummaries) {
            objects.add(object.getKey());
        }
        s3Client.deleteObjects(new DeleteObjectsRequest(env.getProperty("bucket-name")).withKeys(objects.toArray(new String[0])));
    }

//    @Scheduled(cron = "0 */5 * * * *")
    public synchronized void checkForNewData() {
        log.info("||| checkForNewData");
        String bucket = env.getProperty("bucket-name");
        if(s3Client.listObjectsV2(bucket, "cromwell-execution/").getKeyCount() > 0) {
            // TODO maybe limit the amount of job entries checked for better performance?
            SqlRowSet jobsToCheck = auroraDBManager.checkJobsStatus();
            List<String> jobsToUpdate = new LinkedList<>();
            String cromwellId, sampleId, outputLink, outputKey;
            JSONObject requestedOutputs;
            Iterator<String> outputKeys;
            while (jobsToCheck.next()) {

                try {
                    HttpResponse<JsonNode> response = Unirest
                            .get(String.format("%s/query?label=jobId:%s", env.getProperty("ec2-cromwell-dns"), jobsToCheck.getString("JobID")))
                            .header("accept", "application/json")
                            .asJson();
                    if (response.getStatus() / 100 != 2)
                        throw new Exception(response.getStatusText());
//                    log.info(response.getBody());

                    // id, name, status
                    JSONObject jobDetails = response.getBody().getObject().getJSONArray("results").getJSONObject(0);
                    if (!jobDetails.getString("status").equals("Succeeded")) continue;
//                    log.info(jobDetails.toString());

                    cromwellId = jobDetails.getString("id");
                    sampleId = jobsToCheck.getString("SampleID");
                    requestedOutputs = new JSONObject(jobsToCheck.getString("Output"));
//                    log.info("cromwellId: " + cromwellId);
//                    log.info("sampleId: " + sampleId);
//                    log.info("output: " + jobsToCheck.getString("Output"));

                    response = Unirest
                            .get(String.format("%s/%s/outputs", env.getProperty("ec2-cromwell-dns"), cromwellId))
                            .header("accept", "application/json")
                            .asJson();
                    if (response.getStatus() / 100 != 2)
                        throw new Exception(response.getStatusText());
//                    log.info(response.getBody());
                    JSONObject availableOutputs = response.getBody().getObject().getJSONObject("outputs");
//                    log.info(availableOutputs.toString());

                    outputKeys = requestedOutputs.keys();
//                    log.info(outputKeys.toString());
                    // TODO make nested try-catch to take care of s3client errors while moving objects
                    while (outputKeys.hasNext()) {
                        outputLink = outputKeys.next();
                        outputKey = availableOutputs.getString(outputLink).replaceFirst(".*/(cromwell-execution/.*)", "$1");
                        s3Client.copyObject(bucket, outputKey, bucket, String.format("samples/%s/%s", sampleId, requestedOutputs.getString(outputLink)));
                        s3Client.deleteObject(bucket, outputKey);
//                        log.info(outputLink + " | " + outputKey);
//                        log.info(bucket + " | " + String.format("samples/%s/%s", sampleId, requestedOutputs.getString(outputLink)));
                    }
                    deleteAllObjectsWithPrefix(String.format("cromwell-execution/%s/%s/", jobDetails.getString("name"), cromwellId));

                    jobsToUpdate.add(jobsToCheck.getString("JobID"));
//                    log.info("");
                } catch (Exception e) {
                    if (!jobsToUpdate.isEmpty())
                        auroraDBManager.updateJobsStatus(jobsToUpdate);
                    log.info(e.getMessage());
                }
            }
            if (!jobsToUpdate.isEmpty())
                auroraDBManager.updateJobsStatus(jobsToUpdate);
        }
    }

}
