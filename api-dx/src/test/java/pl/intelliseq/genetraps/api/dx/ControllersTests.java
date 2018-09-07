package pl.intelliseq.genetraps.api.dx;

import com.dnanexus.DXJob;
import com.dnanexus.JobState;
import com.fasterxml.jackson.databind.JsonNode;
import lombok.extern.log4j.Log4j2;
import org.hamcrest.Matchers;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.test.context.SpringBootTest;
import org.springframework.boot.test.context.SpringBootTest.WebEnvironment;
import org.springframework.boot.test.web.client.TestRestTemplate;
import org.springframework.http.HttpEntity;
import org.springframework.http.HttpHeaders;
import org.springframework.http.MediaType;
import org.springframework.test.context.ActiveProfiles;
import org.springframework.test.context.junit4.SpringRunner;
import org.springframework.util.LinkedMultiValueMap;
import org.springframework.util.MultiValueMap;
import pl.intelliseq.genetraps.api.dx.helpers.DxApiProcessManager;

import static org.junit.Assert.assertThat;


@RunWith(SpringRunner.class)
@SpringBootTest(webEnvironment = WebEnvironment.RANDOM_PORT)
@Log4j2
@ActiveProfiles("test")
public class ControllersTests {

    @Autowired
    private TestRestTemplate restTemplate;

    @Autowired
    private DxApiProcessManager processManager;

    private String sampleLeft = "http://resources.intelliseq.pl/kamilant/test-data/fastq/capn3.1.fq.gz";
    private String sampleRight = "http://resources.intelliseq.pl/kamilant/test-data/fastq/capn3.2.fq.gz";
    private String interval = "chr15:42377802-42397802";

    private Integer mkDir() {
        return restTemplate.getForObject("/mkdir", JsonNode.class).get("response").asInt();
    }

    private DXJob.Describe waitUntilJobIsDone(String jobId) {
        try {
            JobState jobState;
            DXJob.Describe describe;
            do {
                describe = DXJob.getInstance(jobId).describe();
                jobState = describe.getState();
                Thread.sleep(5000);
            } while (jobState != JobState.DONE);
            return describe;
        } catch (InterruptedException e) {
            throw new RuntimeException(e);
        }
    }

    private DXJob.Describe upload(String sampleUrl, Integer sampleid, String... tags) {
        HttpHeaders uploadHeaders = new HttpHeaders();
        uploadHeaders.setContentType(MediaType.APPLICATION_FORM_URLENCODED);

        MultiValueMap<String, String> uploadValueMap = new LinkedMultiValueMap<String, String>();
        uploadValueMap.add("url", sampleUrl);
        uploadValueMap.add("sampleid", sampleid.toString());
        for (String tag : tags) {
            uploadValueMap.add("tag", tag);
        }

        HttpEntity<MultiValueMap<String, String>> uploadEntity = new HttpEntity<MultiValueMap<String, String>>(uploadValueMap, uploadHeaders);

        JsonNode response = this.restTemplate.postForObject("/upload", uploadEntity, JsonNode.class);
        assertThat(response.get("id").textValue(), Matchers.matchesPattern("job-\\w*"));

        return waitUntilJobIsDone(response.get("id").textValue());
    }

    private DXJob.Describe fastqc(String fileId) {
        HttpHeaders uploadHeaders = new HttpHeaders();
        uploadHeaders.setContentType(MediaType.APPLICATION_FORM_URLENCODED);

        MultiValueMap<String, String> uploadValueMap = new LinkedMultiValueMap<String, String>();
        uploadValueMap.add("fileId", fileId);

        HttpEntity<MultiValueMap<String, String>> uploadEntity = new HttpEntity<MultiValueMap<String, String>>(uploadValueMap, uploadHeaders);

        JsonNode response = this.restTemplate.postForObject("/fastqc", uploadEntity, JsonNode.class);
        assertThat(response.get("id").textValue(), Matchers.matchesPattern("job-\\w*"));

        return waitUntilJobIsDone(response.get("id").textValue());
    }

    private DXJob.Describe bwa(Integer sampleid) {
        HttpHeaders uploadHeaders = new HttpHeaders();
        uploadHeaders.setContentType(MediaType.APPLICATION_FORM_URLENCODED);

        MultiValueMap<String, String> uploadValueMap = new LinkedMultiValueMap<>();
        uploadValueMap.add("sampleid", sampleid.toString());

        HttpEntity<MultiValueMap<String, String>> uploadEntity = new HttpEntity<>(uploadValueMap, uploadHeaders);

        JsonNode response = this.restTemplate.postForObject("/bwa", uploadEntity, JsonNode.class);
        assertThat(response.get("id").textValue(), Matchers.matchesPattern("job-\\w*"));

        return waitUntilJobIsDone(response.get("id").textValue());
    }

    private DXJob.Describe gatkhc(Integer sampleid, String interval) {
        HttpHeaders uploadHeaders = new HttpHeaders();
        uploadHeaders.setContentType(MediaType.APPLICATION_FORM_URLENCODED);

        MultiValueMap<String, String> uploadValueMap = new LinkedMultiValueMap<>();
        uploadValueMap.add("sampleid", sampleid.toString());
        uploadValueMap.add("interval", interval);

        HttpEntity<MultiValueMap<String, String>> uploadEntity = new HttpEntity<>(uploadValueMap, uploadHeaders);

        JsonNode response = this.restTemplate.postForObject("/gatkhc", uploadEntity, JsonNode.class);
        assertThat(response.get("id").textValue(), Matchers.matchesPattern("job-\\w*"));

        return waitUntilJobIsDone(response.get("id").textValue());
    }

    @Test
    public void upload() {
        //Integer sampleid = mkDir();
        //DXJob.Describe upload1 = upload(sampleLeft, sampleid, "left");
        //DXJob.Describe upload2 = upload(sampleRight, sampleid, "right");
        //String file1Id = upload1.getOutput(JsonNode.class).get("file").get("$dnanexus_link").asText();

        //log.info(restTemplate.getForObject("/describe/" + file1Id, JsonNode.class));
        //log.info(fastqc(file1Id));
        //log.info(bwa(sampleid));
        //log.info(gatkhc(sampleid, interval));

    }




}
