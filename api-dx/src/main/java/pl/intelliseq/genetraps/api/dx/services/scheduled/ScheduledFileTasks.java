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

    public void deleteAllObjectsWithPrefix(String prefix) {
        List<S3ObjectSummary> objectSummaries = s3Client.listObjectsV2(env.getProperty("bucket-name"), prefix).getObjectSummaries();
        List<String> objects = new ArrayList<>();
        for (S3ObjectSummary object : objectSummaries) {
            objects.add(object.getKey());
        }
        s3Client.deleteObjects(new DeleteObjectsRequest(env.getProperty("bucket-name")).withKeys(objects.toArray(new String[0])));
    }

    @Scheduled(cron = "0 */5 * * * *")
    public synchronized void checkForOutput() {
        log.info("||| checkForNewData");
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

                    RowCallbackHandler rch = resultSet -> {
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
                            String sampleId = resultSet.getString("SampleID");
                            JSONObject requestedOutputs = new JSONObject(resultSet.getString("Output"));

                            // get the outputs of the job's workflow
                            response = Unirest
                                    .get(String.format("%s/%s/outputs", env.getProperty("ec2-cromwell-dns"), cromwellId))
                                    .header("accept", "application/json")
                                    .asJson();
                            if (response.getStatus() / 100 != 2)
                                throw new Exception(response.getStatusText());
                            JSONObject availableOutputs = response.getBody().getObject().getJSONObject("outputs");

                            Iterator<String> outputKeys = requestedOutputs.keys();
                            String outputLink, outputKey;

                            // move the requested outputs to sample directory of given sampleId
                            // TODO make nested try-catch to take care of s3client errors while moving objects
                            while (outputKeys.hasNext()) {
                                outputLink = outputKeys.next();
                                outputKey = availableOutputs.getString(outputLink).replaceFirst(".*/(cromwell-execution/.*)", "$1");
                                s3Client.copyObject(bucket, outputKey, bucket, String.format("samples/%s/%s", sampleId, requestedOutputs.getString(outputLink)));
                                s3Client.deleteObject(bucket, outputKey);
                            }
                            // delete remaining outputs of a workflow
                            deleteAllObjectsWithPrefix(String.format("cromwell-execution/%s/%s/", jobDetails.getString("name"), cromwellId));

                            log.info(resultSet.getInt("JobStatus"));
                            resultSet.updateInt("JobStatus", 1);

                        } catch (Exception e) {
                            log.info(e.getMessage());
                        } finally {
                            resultSet.updateInt("SampleID", resultSet.getInt("SampleID"));  // mock update to make sure resultSet is always an updated one
                            resultSet.updateRow();
                        }
                    };

//                    jdbcTemplate.setMaxRows(1);
                    jdbcTemplate.query(pscf.newPreparedStatementCreator(new Object[]{jobId}), rch);

                } catch (Exception e) {
                    log.info(e.getMessage());
                }
        }

    }

}