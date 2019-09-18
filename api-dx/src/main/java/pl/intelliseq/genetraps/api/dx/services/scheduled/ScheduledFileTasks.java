package pl.intelliseq.genetraps.api.dx.services.scheduled;

import com.amazonaws.services.s3.AmazonS3;
import com.amazonaws.services.s3.AmazonS3ClientBuilder;
import com.amazonaws.services.s3.model.DeleteObjectsRequest;
import com.amazonaws.services.s3.model.S3ObjectSummary;
import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.node.ObjectNode;
import com.mashape.unirest.http.HttpResponse;
import com.mashape.unirest.http.Unirest;
import com.mashape.unirest.http.exceptions.UnirestException;
import lombok.extern.log4j.Log4j2;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.core.env.Environment;
import org.springframework.jdbc.core.JdbcTemplate;
import org.springframework.jdbc.core.PreparedStatementCreatorFactory;
import org.springframework.jdbc.core.RowCallbackHandler;
import org.springframework.stereotype.Component;
import pl.intelliseq.genetraps.api.dx.helpers.AuroraDBManager;

import javax.annotation.PostConstruct;
import javax.sql.DataSource;
import java.io.IOException;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Types;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

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

//    @Scheduled(cron = "0 */5 * * * *")
    private synchronized void checkForOutput() {

        log.info("||| checkForOutput");
        String bucket = env.getProperty("bucket.default");
        String debugDir = env.getProperty("debug.folder");

        List<String> jobsToCheck = auroraDBManager.getJobs(true);

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
                            HttpResponse<com.mashape.unirest.http.JsonNode> response = getCromwellJobWithLabel("jobId", resultSet.getString("JobID"));

                            // extract job details from the response
                            // id, name, status
                            JsonNode jobDetails = getJobDetails(response);
                            String cromwellJobId = getCromwellJobId(jobDetails);
                            String workflowName = getWorkflowName(jobDetails);

                            // check status of job execution
                            String jobStatus = getJobStatus(jobDetails);
                            if (!jobStatus.equals("Succeeded")) {
                                if(jobStatus.equals("Failed")) {
//                                    moveToDebugWhenJobFailed();
                                    //TODO: transfer not delete error contents
                                    deleteAllObjectsWithPrefix(bucket, String.format("%s/%s/%s/", env.getProperty("cromwell.execution.folder"), workflowName, cromwellJobId));
                                    throw new Exception("Job execution - failed");
                                }
                                else
                                    log.info(cromwellJobId + ": job still running probably");
                                    return;
                            }

                            // get the outputs of the job's workflow
//                            response = getCromwellJobOutputs(cromwellJobId);

//                            JSONObject availableOutputs = getAvailableOuputsForJob(env.getProperty("bucket.default"), env.getProperty("cromwell.execution.folder"), getWorkflowName(jobDetails), cromwellJobId);
                            JsonNode availableOutputs = getOutputsForJob(cromwellJobId);
                            JsonNode requestedOutputs;
                            try {
                                requestedOutputs = getRequestedOutputsForJob(resultSet);
                            } catch (SQLException | IOException e) {
                                log.info(e.toString());
                                return;
                            }

                            // move the requested outputs to sample directory with given sampleId
                            if(requestedOutputs.size() == 0) {
                                moveAllOutputsToSample(bucket, availableOutputs, resultSet.getString("SampleI"), workflowName, jobId);
//                                moveAllObjects(bucket, availableOutputs, resultSet.getString("SampleID"));
                            }
                            else {
                                moveReqOutputsToSample(bucket, availableOutputs, requestedOutputs, resultSet.getString("SampleID"), workflowName, jobId);
                            }

                            // delete remaining outputs of a workflow
                            deleteAllObjectsWithPrefix(bucket, String.format("%s/%s/%s/", env.getProperty("cromwell.execution.folder"), workflowName, cromwellJobId));

                            log.info(resultSet.getInt("JobStatus"));
                            // update job status in db to successful (1)
                            resultSet.updateInt("JobStatus", 1);

                        } catch (Exception e) {
//                            log.info(e.getMessage());
//                            try {
//                                moveToDebug(bucket, env.getProperty("debug-folder"), String.format("JOB ID: %s\nSAMPLE ID: %s\nERROR MESSAGE:\n%s", jobId, resultSet.getString("SampleID"), e.getMessage()), resultSet.getString("JobID"));
//                            } catch (Exception eDebug) {
//                                moveToDebugSimple(bucket, jobId, env.getProperty("debug-folder"), String.format("JOB ID: %s\nMOVING TO DEBUG ERROR:\n%sERROR MESSAGE:\n%s", jobId, eDebug.getMessage(), e.getMessage()));
//                            }
//                            resultSet.updateInt("JobStatus", 2);

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

    private HttpResponse<com.mashape.unirest.http.JsonNode> getCromwellJobWithLabel(String label, String key) throws Exception {

        HttpResponse<com.mashape.unirest.http.JsonNode> response = Unirest
                .get(String.format("%s/query?label=%s:%s", env.getProperty("cromwell.server"), label, key))
                .header("accept", "application/json")
                .asJson();
        if (response.getStatus() / 100 != 2)
            throw new Exception(response.getStatusText());
        return response;
    }

    private JsonNode getJobDetails(HttpResponse<com.mashape.unirest.http.JsonNode> response) throws IOException {

        return new ObjectMapper().readValue(response.getRawBody(), ObjectNode.class).get("results").get(0);

//        return response.getBody().getObject().getJSONArray("results").getJSONObject(0);
    }

    private String getCromwellJobId(JsonNode jobDetails) {

        return jobDetails.get("id").asText();
    }

    private String getJobStatus(JsonNode jobDetails) {

        return jobDetails.get("status").asText();
    }

    private String getWorkflowName(JsonNode jobDetails) {

        return jobDetails.get("name").asText();
    }

    private JsonNode getOutputsForJob(String cromwellJobId) throws InterruptedException {
        ObjectNode response;
        try {
            HttpResponse<com.mashape.unirest.http.JsonNode> responseUnirest = Unirest
                    .get(String.format("%s/%s/outputs", env.getProperty("cromwell.server"), cromwellJobId))
                    .header("accept", "application/json")
                    .asJson();
            response = new ObjectMapper().readValue(responseUnirest.getRawBody(), ObjectNode.class);
            if (responseUnirest.getStatus() / 100 != 2)
                throw new InterruptedException(response.toString());
        } catch (UnirestException | IOException e) {
            throw new InterruptedException(e.getMessage());
        }
        return response.get("outputs");
    }

    private JsonNode getRequestedOutputsForJob(ResultSet resultSet) throws SQLException, IOException {

        return new ObjectMapper().readValue(resultSet.getString("Output"), JsonNode.class);
    }

    private void moveAllOutputsToSample(String bucket, JsonNode availableOutputs, String destinationSampleId, String workflowName, String jobId) {

        for (Iterator<Map.Entry<String, JsonNode>> it = availableOutputs.fields(); it.hasNext(); ) {
            String outputKey = it.next().getKey();
            String outputLinkKey = availableOutputs.get(outputKey).asText();
            outputLinkKey = outputLinkKey.substring(outputLinkKey.indexOf(bucket) + bucket.length() + 1);
            s3Client.copyObject(bucket, outputLinkKey, bucket, String.format("%s/%s/%s/%s/%s", env.getProperty("samples.folder"), destinationSampleId, workflowName, jobId, outputLinkKey.substring(outputLinkKey.lastIndexOf('/') + 1)));
            s3Client.deleteObject(bucket, outputLinkKey);
        }
    }

    private void moveReqOutputsToSample(String bucket, JsonNode availableOutputs, JsonNode requestedObjects, String destinationSampleId, String workflowName, String jobId) {

        for (Iterator<Map.Entry<String, JsonNode>> it = availableOutputs.fields(); it.hasNext(); ) {
            String outputKey = it.next().getKey();
            if(!requestedObjects.has(outputKey))
                continue;
            String outputLinkKey = availableOutputs.get(outputKey).asText();
            String newName = requestedObjects.get(outputKey).asText();
            outputLinkKey = outputLinkKey.substring(outputLinkKey.indexOf(bucket) + bucket.length() + 1);
            s3Client.copyObject(bucket, outputLinkKey, bucket, String.format("%s/%s/%s/%s/%s", env.getProperty("samples.folder"), destinationSampleId, workflowName, jobId, newName.isEmpty() ?  outputLinkKey.substring(outputLinkKey.lastIndexOf('/') + 1) : newName));
            s3Client.deleteObject(bucket, outputLinkKey);
        }
    }

    private void deleteAllObjectsWithPrefix(String bucket, String prefix) {

        List<S3ObjectSummary> objectSummaries = s3Client.listObjectsV2(env.getProperty("bucket.default"), prefix).getObjectSummaries();
        List<String> objects = new ArrayList<>();
        for (S3ObjectSummary object : objectSummaries) {
            objects.add(object.getKey());
        }
        s3Client.deleteObjects(new DeleteObjectsRequest(bucket).withKeys(objects.toArray(new String[0])));
    }

//    private void moveToDebug(String sourceBucket, String debugFolder, String errorMessage, String jobId) throws Exception {
//
//        HttpResponse<JsonNode> response = getCromwellJobWithLabel("jobId", jobId);
//        JSONObject jobDetails = getJobDetails(response);
//        String cromwellJobId = getCromwellJobId(jobDetails);
//        response = getCromwellJobOutputs(cromwellJobId);
//
//        JSONObject objectsToMove = response.getBody().getObject().getJSONObject("outputs");
//
//        Iterator<String> movedObjects = objectsToMove.keys();
//        String outputLink, outputKey;
//
//        while (movedObjects.hasNext()) {
//            outputLink = movedObjects.next();
//            outputKey = objectsToMove.getString(outputLink).replaceFirst(".*/(cromwell-execution/.*)", "$1");
//            log.info(outputKey);
//            //TODO: not all outputs are listed in Outputs of workflow details, need modification...
//            try {
//                s3Client.copyObject(sourceBucket, outputKey, sourceBucket, String.format("%s/%s", debugFolder, outputKey));
//                s3Client.deleteObject(sourceBucket, outputKey);
//            } catch (Exception e) {
//                log.info("err");
//            }
//        }
//        // TODO: put to finally
//        s3Client.putObject(sourceBucket, String.format("%s/%s/%s", debugFolder, jobId, "log.err"), errorMessage);
//    }

    private void moveToDebugSimple(String bucket, String jobId, String debugFolder, String errorMessage) {

        s3Client.putObject(bucket, String.format("%s/%s/%s", debugFolder, jobId, "log.err"), errorMessage);
    }


//    private void moveAllObjects(String bucket, JSONObject availableObjects, String destinationSampleId) {
//
//        Iterator<String> allObjectsIt = availableObjects.keys();
//        String object, objectKey;
//        while (allObjectsIt.hasNext()) {
//            object = allObjectsIt.next();
//            objectKey = availableObjects.getString(object).replaceFirst(".*/(cromwell-execution/.*)", "$1");
//            s3Client.copyObject(bucket, objectKey, bucket, String.format("%s/%s/%s", env.getProperty("samples.folder"), destinationSampleId, object));
//            s3Client.deleteObject(bucket, objectKey);
//        }
//    }

//    private HttpResponse<JsonNode> getCromwellJobOutputs(String cromwellJobId) throws Exception {
//
//        HttpResponse<com.mashape.unirest.http.JsonNode> response = Unirest
//                .get(String.format("%s/%s/outputs", env.getProperty("cromwell.server"), cromwellJobId))
//                .header("accept", "application/json")
//                .asJson();
//        if (response.getStatus() / 100 != 2)
//            throw new Exception(response.getStatusText());
//        return response;
//    }

//    private JSONObject getAvailableOuputsForJob(String bucket, String cromwellExecutionFolder, String workflowName, String jobId) {
//
//        List<S3ObjectSummary> objectSummaries = s3Client.listObjectsV2(bucket, String.format("%s/%s/%s/", cromwellExecutionFolder, workflowName, jobId)).getObjectSummaries();
////            log.info(dir.isEmpty() ? String.format("samples/%s/", sampleId) : String.format("samples/%s/%s/", sampleId, dir));
//        JSONObject objects = new JSONObject();
//        String objectKey;
//        for (S3ObjectSummary object : objectSummaries) {
//            objectKey = object.getKey();
//            if(objectKey.matches(".*/"))    continue;       // omits id of directories
//            objects.put(objectKey.substring(objectKey.lastIndexOf("/") + 1), String.format("/%s", objectKey));
//        }
//
//        return objects;
//    }
}
