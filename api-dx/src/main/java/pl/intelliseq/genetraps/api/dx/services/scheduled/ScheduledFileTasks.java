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
import org.springframework.jdbc.core.JdbcTemplate;
import org.springframework.jdbc.core.PreparedStatementCreatorFactory;
import org.springframework.jdbc.core.RowCallbackHandler;
import org.springframework.scheduling.annotation.Scheduled;
import org.springframework.stereotype.Component;
import pl.intelliseq.genetraps.api.dx.helpers.AuroraDBManager;

import javax.annotation.PostConstruct;
import javax.sql.DataSource;
import java.sql.Types;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

@Log4j2
@Component
public class ScheduledFileTasks {

    @Autowired
    private Environment env;

    @Autowired
    private AuroraDBManager auroraDBManager;

    final private AmazonS3 s3Client = AmazonS3ClientBuilder.defaultClient();

    private JdbcTemplate jdbcTemplate;

    @Autowired
    private DataSource dataSource;

    @PostConstruct
    private void postConstruct() {
        jdbcTemplate = new JdbcTemplate(dataSource);
    }

    public void deleteAllObjectsWithPrefix(String bucket, String prefix) {
        List<S3ObjectSummary> objectSummaries = s3Client.listObjectsV2(env.getProperty("bucket-name"), prefix).getObjectSummaries();
        List<String> objects = new ArrayList<>();
        for (S3ObjectSummary object : objectSummaries) {
            objects.add(object.getKey());
        }
        s3Client.deleteObjects(new DeleteObjectsRequest(bucket).withKeys(objects.toArray(new String[0])));
    }

    // TODO: make cron delete given record and all its workflow output if e.g. number of failed tries reaches 3/5/...
    @Scheduled(cron = "0 */5 * * * *")
    private synchronized void checkForOutput() {
        log.info("||| checkForOutput");
        String bucket = env.getProperty("bucket-name");
        if(s3Client.listObjectsV2(bucket, "cromwell-execution/").getKeyCount() > 0) {

            List<String> jobsToCheck = auroraDBManager.getJobsWithStatusNotYetSucceded();

            for(String jobId : jobsToCheck)
                try {
                    // row lock
                    PreparedStatementCreatorFactory pscf = new PreparedStatementCreatorFactory(
                            "SELECT JobId, SampleId, JobStatus, Output FROM Jobs " +
                                    "WHERE JobID = ? AND JobStatus = 0 FOR UPDATE", Types.CHAR);
                    pscf.setUpdatableResults(true);
//                    pscf.setResultSetType(ResultSet.TYPE_FORWARD_ONLY);

                    RowCallbackHandler rowHandler = resultSet -> {
                        try {
                            // query to cromwell to get details about job with 'jobId' as a label
                            HttpResponse<JsonNode> response = Unirest
                                    .get(String.format("%s/query?label=jobId:%s", env.getProperty("ec2-cromwell-dns"), resultSet.getString("JobID")))
                                    .header("accept", "application/json")
                                    .asJson();
                            if (response.getStatus() / 100 != 2)
                                throw new Exception(response.getStatusText());

                            // extract job details from the response
                            // id, name, status
                            JSONObject jobDetails = response.getBody().getObject().getJSONArray("results").getJSONObject(0);

                            if (!jobDetails.getString("status").equals("Succeeded")) {     // Check whether is successful
                                return;
                            }

                            String cromwellId = jobDetails.getString("id");

                            // get the outputs of the job's workflow
                            response = Unirest
                                    .get(String.format("%s/%s/outputs", env.getProperty("ec2-cromwell-dns"), cromwellId))
                                    .header("accept", "application/json")
                                    .asJson();
                            if (response.getStatus() / 100 != 2)
                                throw new Exception(response.getStatusText());

                            JSONObject availableOutputs = response.getBody().getObject().getJSONObject("outputs");
                            JSONObject requestedOutputs = new JSONObject(resultSet.getString("Output"));

                            // move the requested outputs to sample directory with given sampleId
                            moveObjectsWithJSON(bucket, availableOutputs, requestedOutputs, "samples/" + resultSet.getString("SampleID"));

                            // delete remaining outputs of a workflow
                            deleteAllObjectsWithPrefix(bucket, String.format("cromwell-execution/%s/%s/", jobDetails.getString("name"), cromwellId));

                            log.info(resultSet.getInt("JobStatus"));
                            resultSet.updateInt("JobStatus", 1);

                        } catch (Exception e) {
                            log.info(e.getMessage());
                            try {
                                moveToDebug(bucket, "debug", e.getMessage(), resultSet.getString("JobID"));
                            } catch (Exception e1) {
                                e1.printStackTrace();
                            }
                            resultSet.updateInt("JobStatus", 2);

                        } finally {
                            resultSet.updateInt("SampleID", resultSet.getInt("SampleID"));  // mock update to make sure resultSet is always an updated one
                            resultSet.updateRow();
                        }
                    };

//                    jdbcTemplate.setMaxRows(1);
                    jdbcTemplate.query(pscf.newPreparedStatementCreator(new Object[]{jobId}), rowHandler);

                } catch (Exception e) {
                    log.info(e.getMessage());
                }
        }
    }

    private void moveObjectsWithJSON(String bucket, JSONObject availableObjects, JSONObject requestedObjects, String destinationDirectory) {
        Iterator<String> movedObjects = requestedObjects.keys();
        String outputLink, outputKey, newName;
        // TODO make nested try-catch to take care of s3client errors while moving objects
        while (movedObjects.hasNext()) {
            outputLink = movedObjects.next();
            outputKey = availableObjects.getString(outputLink).replaceFirst(".*/(cromwell-execution/.*)", "$1");
            newName = requestedObjects.getString(outputLink);
            s3Client.copyObject(bucket, outputKey, bucket, String.format("%s/%s", destinationDirectory, newName.length() == 0 ? outputLink : newName));
            s3Client.deleteObject(bucket, outputKey);
        }
    }

    private void moveToDebug(String sourceBucket, String debugDirectory, String errorMessage, String jobId) throws Exception {

        HttpResponse<JsonNode> response = Unirest
                .get(String.format("%s/query?label=jobId:%s", env.getProperty("ec2-cromwell-dns"), jobId))
                .header("accept", "application/json")
                .asJson();
        if (response.getStatus() / 100 != 2)
            throw new Exception(response.getStatusText());
        JSONObject jobDetails = response.getBody().getObject().getJSONArray("results").getJSONObject(0);
        String cromwellId = jobDetails.getString("id");
        response = Unirest
                .get(String.format("%s/%s/outputs", env.getProperty("ec2-cromwell-dns"), cromwellId))
                .header("accept", "application/json")
                .asJson();
        if (response.getStatus() / 100 != 2)
            throw new Exception(response.getStatusText());

        JSONObject objectsToMove = response.getBody().getObject().getJSONObject("outputs");

        Iterator<String> movedObjects = objectsToMove.keys();
        String outputLink, outputKey;

        while (movedObjects.hasNext()) {
            outputLink = movedObjects.next();
            outputKey = objectsToMove.getString(outputLink).replaceFirst(".*/(cromwell-execution/.*)", "$1");
            log.info(outputKey);
            //TODO: not all outputs are listed in Outputs of workflow details, need modification...
            try {
                s3Client.copyObject(sourceBucket, outputKey, sourceBucket, String.format("%s/%s", debugDirectory, outputKey));
                s3Client.deleteObject(sourceBucket, outputKey);
            } catch (Exception e) {
                log.info("err");
            }
        }

        // TODO: put to finally
        s3Client.putObject(sourceBucket, jobId + "\n" + debugDirectory + "/log.err", errorMessage);
    }

}