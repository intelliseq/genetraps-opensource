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
import java.sql.ResultSet;
import java.sql.SQLException;
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

    // TODO: make cron delete given record and all its workflow output if e.g. number of failed tries reaches 3/5/...
    @Scheduled(cron = "0 */5 * * * *")
    private synchronized void checkForOutput() {

        log.info("||| checkForOutput");
        String bucket = env.getProperty("default-bucket");
        String debugDir = env.getProperty("debug-folder");

        List<String> jobsToCheck = auroraDBManager.getJobsWithStatusNotYetSucceded();

        if(!jobsToCheck.isEmpty()) {

//            s3Client.listObjectsV2(bucket, "cromwell-execution/").getKeyCount() > 0

            for(String jobId : jobsToCheck)

                try {
                    // row lock
                    PreparedStatementCreatorFactory pscf = new PreparedStatementCreatorFactory(
                            "SELECT JobId, SampleId, JobStatus, Output FROM Jobs " +
                                    "WHERE JobID = ? AND JobStatus = 0 FOR UPDATE", Types.CHAR);
                    pscf.setUpdatableResults(true);
                    //TODO important
//                    pscf.setResultSetType(ResultSet.TYPE_FORWARD_ONLY);

                    RowCallbackHandler rowHandler = resultSet -> {
                        try {
                            // query to cromwell to get details about job with 'jobId' as a label
                            HttpResponse<JsonNode> response = getCromwellJobWithLabel("jobId", resultSet.getString("JobID"));

                            // extract job details from the response
                            // id, name, status
                            JSONObject jobDetails = getJobDetails(response);

                            // check status of job execution
                            String jobStatus = getJobStatus(jobDetails);
                            if (!jobStatus.equals("Succeeded")) {
                                if(jobStatus.equals("Failed")) {
//                                    moveToDebugWhenJobFailed();
                                    throw new Exception("Job execution - failed");
                                }
                                else
                                    return;
                            }

                            String cromwellJobId = getCromwellJobId(jobDetails);

                            // get the outputs of the job's workflow
                            response = getCromwellJobOutputs(cromwellJobId);

                            JSONObject availableOutputs = getAvailableOuputsForJob(response);
                            JSONObject requestedOutputs = getRequestedOutputsForJob(resultSet);

                            // move the requested outputs to sample directory with given sampleId
                            if(requestedOutputs.length() == 0) {
                                moveAllObjects(bucket, availableOutputs, resultSet.getString("SampleID"));
                            }
                            else {
                                moveObjectsWithJSON(bucket, availableOutputs, requestedOutputs, resultSet.getString("SampleID"));
                            }

                            // delete remaining outputs of a workflow
                            deleteAllObjectsWithPrefix(bucket, String.format("%s/%s/%s/", env.getProperty("cromwell-execution-folder"), jobDetails.getString("name"), cromwellJobId));

                            log.info(resultSet.getInt("JobStatus"));
                            // update job status in db to successful (1)
                            resultSet.updateInt("JobStatus", 1);

                        } catch (Exception e) {
//                            log.info(e.getMessage());
                            try {
                                moveToDebug(bucket, env.getProperty("debug-folder"), String.format("JOB ID: %s\nSAMPLE ID: %s\nERROR MESSAGE:\n%s", jobId, resultSet.getString("SampleID"), e.getMessage()), resultSet.getString("JobID"));
                            } catch (Exception eDebug) {
                                moveToDebugSimple(bucket, jobId, env.getProperty("debug-folder"), String.format("JOB ID: %s\nMOVING TO DEBUG ERROR:\n%sERROR MESSAGE:\n%s", jobId, eDebug.getMessage(), e.getMessage()));
                            }
                            resultSet.updateInt("JobStatus", 2);

                        } finally {
                            resultSet.updateInt("SampleID", resultSet.getInt("SampleID"));  // mock update to make sure resultSet is always an updated one, otherwise the lock won't be removed (I know, so effed)
                            resultSet.updateRow();
                        }
                    };

                    //TODO important
//                    jdbcTemplate.setMaxRows(1);
                    jdbcTemplate.query(pscf.newPreparedStatementCreator(new Object[]{jobId}), rowHandler);

                } catch (Exception e) {
                    log.info(e.getMessage());
                }
        }
    }

    private void moveAllObjects(String bucket, JSONObject availableObjects, String destinationSampleId) {

        Iterator<String> allObjectsIt = availableObjects.keys();
        String object, objectKey;
        while (allObjectsIt.hasNext()) {
            object = allObjectsIt.next();
            objectKey = availableObjects.getString(object).replaceFirst(".*/(cromwell-execution/.*)", "$1");
            s3Client.copyObject(bucket, objectKey, bucket, String.format("%s/%s/%s", env.getProperty("samples-folder"), destinationSampleId, object));
            s3Client.deleteObject(bucket, objectKey);
        }
    }

//    private void moveToDebugWhenJobFailed()

    private HttpResponse<JsonNode> getCromwellJobWithLabel(String label, String key) throws Exception {

        HttpResponse<JsonNode> response = Unirest
                .get(String.format("%s/query?label=%s:%s", env.getProperty("ec2-cromwell-dns"), label, key))
                .header("accept", "application/json")
                .asJson();
        if (response.getStatus() / 100 != 2)
            throw new Exception(response.getStatusText());
        return response;
    }

    private JSONObject getJobDetails(HttpResponse<JsonNode> response) {

        return response.getBody().getObject().getJSONArray("results").getJSONObject(0);
    }

    private String getJobStatus(JSONObject jobDetails) {

        return jobDetails.getString("status");
    }

    private String getCromwellJobId(JSONObject jobDetails) {

        return jobDetails.getString("id");
    }

    private HttpResponse<JsonNode> getCromwellJobOutputs(String cromwellJobId) throws Exception {

        HttpResponse<JsonNode> response = Unirest
                .get(String.format("%s/%s/outputs", env.getProperty("ec2-cromwell-dns"), cromwellJobId))
                .header("accept", "application/json")
                .asJson();
        if (response.getStatus() / 100 != 2)
            throw new Exception(response.getStatusText());
        return response;
    }

    private JSONObject getAvailableOuputsForJob(HttpResponse<JsonNode> response) {

        return response.getBody().getObject().getJSONObject("outputs");
    }

    private JSONObject getRequestedOutputsForJob(ResultSet resultSet) throws SQLException {

        return new JSONObject(resultSet.getString("Output"));
    }

    private void moveObjectsWithJSON(String bucket, JSONObject availableObjects, JSONObject requestedObjects, String destinationSampleId) {

        Iterator<String> reqObjectsIt = requestedObjects.keys();
        String requestedObject, reqObjectKey, newName;
        while (reqObjectsIt.hasNext()) {
            requestedObject = reqObjectsIt.next();
            reqObjectKey = availableObjects.getString(requestedObject).replaceFirst(".*/(cromwell-execution/.*)", "$1");
            newName = requestedObjects.getString(requestedObject);  // newName is a value of pair where requestedObject is key
            s3Client.copyObject(bucket, reqObjectKey, bucket, String.format("%s/%s/%s", env.getProperty("samples-folder"), destinationSampleId, newName.isEmpty() ? requestedObject : newName));
            s3Client.deleteObject(bucket, reqObjectKey);
        }
    }

    public void deleteAllObjectsWithPrefix(String bucket, String prefix) {

        List<S3ObjectSummary> objectSummaries = s3Client.listObjectsV2(env.getProperty("default-bucket"), prefix).getObjectSummaries();
        List<String> objects = new ArrayList<>();
        for (S3ObjectSummary object : objectSummaries) {
            objects.add(object.getKey());
        }
        s3Client.deleteObjects(new DeleteObjectsRequest(bucket).withKeys(objects.toArray(new String[0])));
    }

    private void moveToDebug(String sourceBucket, String debugFolder, String errorMessage, String jobId) throws Exception {

        HttpResponse<JsonNode> response = getCromwellJobWithLabel("jobId", jobId);
        JSONObject jobDetails = getJobDetails(response);
        String cromwellJobId = getCromwellJobId(jobDetails);
        response = getCromwellJobOutputs(cromwellJobId);

        JSONObject objectsToMove = response.getBody().getObject().getJSONObject("outputs");

        Iterator<String> movedObjects = objectsToMove.keys();
        String outputLink, outputKey;

        while (movedObjects.hasNext()) {
            outputLink = movedObjects.next();
            outputKey = objectsToMove.getString(outputLink).replaceFirst(".*/(cromwell-execution/.*)", "$1");
            log.info(outputKey);
            //TODO: not all outputs are listed in Outputs of workflow details, need modification...
            try {
                s3Client.copyObject(sourceBucket, outputKey, sourceBucket, String.format("%s/%s", debugFolder, outputKey));
                s3Client.deleteObject(sourceBucket, outputKey);
            } catch (Exception e) {
                log.info("err");
            }
        }
        // TODO: put to finally
        s3Client.putObject(sourceBucket, String.format("%s/%s/%s", debugFolder, jobId, "log.err"), errorMessage);
    }

    private void moveToDebugSimple(String bucket, String jobId, String debugFolder, String errorMessage) {

        s3Client.putObject(bucket, String.format("%s/%s/%s", debugFolder, jobId, "log.err"), errorMessage);
    }
}
