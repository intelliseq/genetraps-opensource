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
import org.springframework.http.*;
import org.springframework.test.context.ActiveProfiles;
import org.springframework.test.context.junit4.SpringRunner;
import org.springframework.util.LinkedMultiValueMap;
import org.springframework.util.MultiValueMap;
import pl.intelliseq.genetraps.api.dx.helpers.DxApiProcessManager;

import static org.hamcrest.CoreMatchers.is;
import static org.hamcrest.text.IsEqualIgnoringCase.equalToIgnoringCase;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertThat;
import static org.junit.Assert.assertTrue;
import static pl.intelliseq.genetraps.api.dx.TestUser.PSYDUCK;


@RunWith(SpringRunner.class)
@SpringBootTest(webEnvironment = WebEnvironment.RANDOM_PORT)
@Log4j2
//@ActiveProfiles("test")
public class ControllersTests {

    @Autowired
    private TestRestTemplate restTemplate;

    @Autowired
    private DxApiProcessManager processManager;

    private String sampleLeft = "http://resources.intelliseq.pl/kamilant/test-data/fastq/capn3.1.fq.gz";
    private String sampleRight = "http://resources.intelliseq.pl/kamilant/test-data/fastq/capn3.2.fq.gz";
    private String interval = "chr15:42377802-42397802";

    public ResponseEntity<JsonNode> getForResponseEnity(TestUser user, String url) {

        HttpHeaders headers = new HttpHeaders();
        headers.set("Authorization", "Bearer "+user.getAccessToken());

        HttpEntity<String> entity = new HttpEntity<>("parameters", headers);

        ResponseEntity<JsonNode> response = restTemplate.exchange(url, HttpMethod.GET, entity, JsonNode.class);

        return response;
    }

    private Integer mkDir() {
        return getForResponseEnity(PSYDUCK, "/mkdir").getBody().get("response").asInt();
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

    private DXJob.Describe upload(TestUser user, String sampleUrl, Integer sampleid, String... tags) {
        HttpHeaders uploadHeaders = new HttpHeaders();
        uploadHeaders.setContentType(MediaType.APPLICATION_FORM_URLENCODED);
        uploadHeaders.set("Authorization", "Bearer "+user.getAccessToken());

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

    private DXJob.Describe fastqc(TestUser user, String fileId) {
        HttpHeaders uploadHeaders = new HttpHeaders();
        uploadHeaders.setContentType(MediaType.APPLICATION_FORM_URLENCODED);
        uploadHeaders.set("Authorization", "Bearer "+user.getAccessToken());

        MultiValueMap<String, String> uploadValueMap = new LinkedMultiValueMap<String, String>();
        uploadValueMap.add("fileId", fileId);

        HttpEntity<MultiValueMap<String, String>> uploadEntity = new HttpEntity<MultiValueMap<String, String>>(uploadValueMap, uploadHeaders);

        JsonNode response = this.restTemplate.postForObject("/fastqc", uploadEntity, JsonNode.class);
        assertThat(response.get("id").textValue(), Matchers.matchesPattern("job-\\w*"));

        return waitUntilJobIsDone(response.get("id").textValue());
    }

    private DXJob.Describe bwa(TestUser user, Integer sampleid) {
        HttpHeaders uploadHeaders = new HttpHeaders();
        uploadHeaders.setContentType(MediaType.APPLICATION_FORM_URLENCODED);
        uploadHeaders.set("Authorization", "Bearer "+user.getAccessToken());

        MultiValueMap<String, String> uploadValueMap = new LinkedMultiValueMap<>();
        uploadValueMap.add("sampleid", sampleid.toString());

        HttpEntity<MultiValueMap<String, String>> uploadEntity = new HttpEntity<>(uploadValueMap, uploadHeaders);

        JsonNode response = this.restTemplate.postForObject("/bwa", uploadEntity, JsonNode.class);
        assertThat(response.get("id").textValue(), Matchers.matchesPattern("job-\\w*"));

        return waitUntilJobIsDone(response.get("id").textValue());
    }

    private DXJob.Describe gatkhc(TestUser user, Integer sampleid, String interval) {
        HttpHeaders uploadHeaders = new HttpHeaders();
        uploadHeaders.setContentType(MediaType.APPLICATION_FORM_URLENCODED);
        uploadHeaders.set("Authorization", "Bearer "+user.getAccessToken());

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
        Integer sampleid = mkDir();
        DXJob.Describe upload1 = upload(PSYDUCK, sampleLeft, sampleid, "left");
        DXJob.Describe upload2 = upload(PSYDUCK, sampleRight, sampleid, "right");

        ResponseEntity<JsonNode> ls = getForResponseEnity(PSYDUCK, String.format("/sample/%s/ls", sampleid));
        assertTrue(ls.getStatusCode().is2xxSuccessful());
        ResponseEntity<JsonNode> lsByNames = getForResponseEnity(PSYDUCK, String.format("/sample/%s/ls?byNames=true", sampleid));
        assertTrue(lsByNames.getStatusCode().is2xxSuccessful());
        log.info(ls.getBody());
        log.info(lsByNames.getBody());

        String file1Id = upload1.getOutput(JsonNode.class).get("file").get("$dnanexus_link").asText();

        ResponseEntity<JsonNode> describe = getForResponseEnity(PSYDUCK, "/describe/" + file1Id);
        assertTrue(describe.getStatusCode().is2xxSuccessful());
        log.info(describe.getBody());

        log.info(fastqc(PSYDUCK, file1Id));
        log.info(bwa(PSYDUCK, sampleid));
        log.info(gatkhc(PSYDUCK, sampleid, interval));

    }




}
