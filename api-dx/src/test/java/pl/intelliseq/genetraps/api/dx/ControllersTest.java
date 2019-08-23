package pl.intelliseq.genetraps.api.dx;

import com.dnanexus.DXJob;
import com.dnanexus.JobState;
import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import lombok.extern.log4j.Log4j2;
import org.hamcrest.Matchers;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.test.autoconfigure.web.servlet.AutoConfigureMockMvc;
import org.springframework.boot.test.context.SpringBootTest;
import org.springframework.boot.test.context.SpringBootTest.WebEnvironment;
import org.springframework.boot.test.web.client.TestRestTemplate;
import org.springframework.http.*;
import org.springframework.mock.web.MockMultipartFile;
import org.springframework.test.context.testng.AbstractTestNGSpringContextTests;
import org.springframework.test.web.servlet.MockMvc;
import org.springframework.test.web.servlet.request.MockMvcRequestBuilders;
import org.springframework.test.web.servlet.result.MockMvcResultMatchers;
import org.springframework.util.LinkedMultiValueMap;
import org.springframework.util.MultiValueMap;
import org.testng.annotations.Test;
import pl.intelliseq.genetraps.api.dx.helpers.DxApiProcessManager;

import static org.hamcrest.CoreMatchers.is;
import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.text.IsEqualIgnoringCase.equalToIgnoringCase;
import static org.testng.Assert.assertEquals;
import static org.testng.AssertJUnit.assertTrue;
import static pl.intelliseq.genetraps.api.dx.UserTest.ADMIN;
import static pl.intelliseq.genetraps.api.dx.UserTest.PSYDUCK;


//@RunWith(SpringRunner.class)
@SpringBootTest(webEnvironment = WebEnvironment.RANDOM_PORT)
@AutoConfigureMockMvc
@Log4j2
//@ActiveProfiles("test")
public class ControllersTest extends AbstractTestNGSpringContextTests {

    @Autowired
    private TestRestTemplate restTemplate;

    @Autowired
    private DxApiProcessManager processManager;

    @Autowired
    private MockMvc mockMvc;

    private String sampleLeft = "http://resources.intelliseq.pl/kamilant/test-data/fastq/capn3.1.fq.gz";
    private String sampleRight = "http://resources.intelliseq.pl/kamilant/test-data/fastq/capn3.2.fq.gz";
    private String interval = "chr15:42377802-42397802";

    public ResponseEntity<JsonNode> getForResponseEnity(UserTest user, String url) {

        HttpHeaders headers = new HttpHeaders();
        headers.set("Authorization", "Bearer " + user.getAccessToken());

        HttpEntity<String> entity = new HttpEntity<>("parameters", headers);

        ResponseEntity<JsonNode> response = restTemplate.exchange(url, HttpMethod.GET, entity, JsonNode.class);

        return response;
    }

    private Integer mkDir() {
        ResponseEntity<JsonNode> out = getForResponseEnity(PSYDUCK, "/sample/new");
        log.debug(out);
        return out.getBody().get("response").asInt();
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

    private DXJob.Describe upload(UserTest user, String sampleUrl, Integer sampleid, String... tags) {
        HttpHeaders uploadHeaders = new HttpHeaders();
        uploadHeaders.setContentType(MediaType.APPLICATION_FORM_URLENCODED);
        uploadHeaders.set("Authorization", "Bearer " + user.getAccessToken());

        MultiValueMap<String, String> uploadValueMap = new LinkedMultiValueMap<String, String>();
        uploadValueMap.add("url", sampleUrl);
        for (String tag : tags) {
            uploadValueMap.add("tag", tag);
        }

        HttpEntity<MultiValueMap<String, String>> uploadEntity = new HttpEntity<MultiValueMap<String, String>>(uploadValueMap, uploadHeaders);

        JsonNode response = this.restTemplate.postForObject(String.format("/sample/%s/urlupload", sampleid.toString()), uploadEntity, JsonNode.class);
        log.debug(response);
        assertThat(response.get("id").textValue(), Matchers.matchesPattern("job-\\w*"));

        return waitUntilJobIsDone(response.get("id").textValue());
    }

//    private DXJob.Describe fastqc(UserTest user, String fileId) {
//        HttpHeaders uploadHeaders = new HttpHeaders();
//        uploadHeaders.setContentType(MediaType.APPLICATION_FORM_URLENCODED);
//        uploadHeaders.set("Authorization", "Bearer " + user.getAccessToken());
//
//        MultiValueMap<String, String> uploadValueMap = new LinkedMultiValueMap<String, String>();
//        uploadValueMap.add("fileId", fileId);
//
//        HttpEntity<MultiValueMap<String, String>> uploadEntity = new HttpEntity<MultiValueMap<String, String>>(uploadValueMap, uploadHeaders);
//
//        JsonNode response = this.restTemplate.postForObject("/fastqc", uploadEntity, JsonNode.class);
//        assertThat(response.get("id").textValue(), Matchers.matchesPattern("job-\\w*"));
//
//        return waitUntilJobIsDone(response.get("id").textValue());
//    }

//    private DXJob.Describe bwa(UserTest user, Integer sampleid) {
//        HttpHeaders uploadHeaders = new HttpHeaders();
//        uploadHeaders.setContentType(MediaType.APPLICATION_FORM_URLENCODED);
//        uploadHeaders.set("Authorization", "Bearer " + user.getAccessToken());
//
//        MultiValueMap<String, String> uploadValueMap = new LinkedMultiValueMap<>();
//        uploadValueMap.add("sampleId", sampleid.toString());
//
//        HttpEntity<MultiValueMap<String, String>> uploadEntity = new HttpEntity<>(uploadValueMap, uploadHeaders);
//
//        JsonNode response = this.restTemplate.postForObject("/bwa", uploadEntity, JsonNode.class);
//        log.debug("R " + response);
//        assertThat(response.get("id").textValue(), Matchers.matchesPattern("job-\\w*"));
//
//        return waitUntilJobIsDone(response.get("id").textValue());
//    }

//    private DXJob.Describe gatkhc(UserTest user, Integer sampleid, String interval) {
//        HttpHeaders uploadHeaders = new HttpHeaders();
//        uploadHeaders.setContentType(MediaType.APPLICATION_FORM_URLENCODED);
//        uploadHeaders.set("Authorization", "Bearer " + user.getAccessToken());
//
//        MultiValueMap<String, String> uploadValueMap = new LinkedMultiValueMap<>();
//        uploadValueMap.add("sampleId", sampleid.toString());
//        uploadValueMap.add("interval", interval);
//
//        HttpEntity<MultiValueMap<String, String>> uploadEntity = new HttpEntity<>(uploadValueMap, uploadHeaders);
//
//        JsonNode response = this.restTemplate.postForObject("/gatkhc", uploadEntity, JsonNode.class);
//        assertThat(response.get("id").textValue(), Matchers.matchesPattern("job-\\w*"));
//
//        return waitUntilJobIsDone(response.get("id").textValue());
//    }

//    @Test
    public void propertiesAndFileUploadTest() {
        Integer sampleId = mkDir();
        try {
            MockMultipartFile multipartFile = new MockMultipartFile("file", "multipart", "text/plain", "multipartTest - good".getBytes());
            String response = new ObjectMapper().readTree(
                    mockMvc.perform(MockMvcRequestBuilders.multipart(String.format("/sample/%s/fileupload", sampleId))
                            .file(multipartFile)
                            .header("Authorization", "Bearer " + PSYDUCK.getAccessToken())
                            .param("newfilename", "multipartTest")
                            .param("tag", "tag1")
                            .param("tag", "tag2"))
                            .andExpect(MockMvcResultMatchers.status().isAccepted())
                            .andReturn().getResponse().getContentAsString()
            ).get("id").textValue();
            assertThat(response, Matchers.matchesPattern("file-\\w*"));
            log.info(response);

            String jsonString = "{\"first\":\"ok\", \"second\":\"notok\"}";
            log.info(jsonString);
            response = new ObjectMapper().readTree(
                    mockMvc.perform(MockMvcRequestBuilders.post(String.format("/sample/%s/properties", sampleId))
                            .contentType(MediaType.APPLICATION_JSON)
                            .content(jsonString).header("Authorization", "Bearer " + PSYDUCK.getAccessToken()))
                            .andExpect(MockMvcResultMatchers.status().isCreated())
                            .andReturn().getResponse().getContentAsString()
            ).toString();
            assertThat(response, Matchers.matchesPattern("\\{\"first\":\"ok\",\"second\":\"notok\"\\}"));
            log.info("POST(new): " + response);
            response = new ObjectMapper().readTree(
                    mockMvc.perform(MockMvcRequestBuilders.get(String.format("/sample/%s/properties", sampleId))
                            .header("Authorization", "Bearer " + PSYDUCK.getAccessToken()))
                            .andReturn().getResponse().getContentAsString()
            ).toString();
            assertThat(response, Matchers.matchesPattern("\\{\"first\":\"ok\",\"second\":\"notok\"\\}"));
            log.info("GET: " + response);

            String jsonStringPost = "{\"third\":\"good\"}";
            response = new ObjectMapper().readTree(
                    mockMvc.perform(MockMvcRequestBuilders.post(String.format("/sample/%s/properties", sampleId))
                            .contentType(MediaType.APPLICATION_JSON)
                            .content(jsonStringPost)
                            .header("Authorization", "Bearer " + PSYDUCK.getAccessToken()))
                            .andExpect(MockMvcResultMatchers.status().isCreated())
                            .andReturn().getResponse().getContentAsString()
            ).toString();
            assertThat(response, Matchers.matchesPattern("\\{\"first\":\"ok\",\"second\":\"notok\",\"third\":\"good\"\\}"));
            log.info("POST(add): " + response);

            String jsonStringPut = "{\"second\":\"ok\"}";
            response = new ObjectMapper().readTree(
                    mockMvc.perform(MockMvcRequestBuilders.put(String.format("/sample/%s/properties", sampleId))
                            .contentType(MediaType.APPLICATION_JSON)
                            .content(jsonStringPut)
                            .header("Authorization", "Bearer " + PSYDUCK.getAccessToken()))
                            .andReturn().getResponse().getContentAsString()
            ).toString();
            assertThat(response, Matchers.matchesPattern("\\{\"first\":\"ok\",\"second\":\"ok\",\"third\":\"good\"\\}"));
            log.info("PUT: " + response);

            String jsonStringDelete = "{\"second\":\"ok\"}";
            response = new ObjectMapper().readTree(
                    mockMvc.perform(MockMvcRequestBuilders.delete(String.format("/sample/%s/properties", sampleId))
                            .contentType(MediaType.APPLICATION_JSON)
                            .content(jsonStringDelete)
                            .header("Authorization", "Bearer " + PSYDUCK.getAccessToken()))
                            .andReturn().getResponse().getContentAsString()
            ).toString();
            assertThat(response, Matchers.matchesPattern("\\{\"first\":\"ok\",\"third\":\"good\"\\}"));
            log.info("DELETE: " + response);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    @Test
    public void priviligesTest() {
        ResponseEntity<JsonNode> psyduck = getForResponseEnity(PSYDUCK, "/user/privileges");
        log.debug(psyduck);
        assertTrue(psyduck.getStatusCode().is2xxSuccessful());

        ResponseEntity<JsonNode> admin = getForResponseEnity(ADMIN, "/user/privileges");
        assertTrue(admin.getStatusCode().is2xxSuccessful());
        admin.getBody().elements().forEachRemaining(n -> assertThat(n.asText(), is(equalToIgnoringCase(Roles.ADMIN.toString()))));


    }

//BUG mkDir can faint - no sync!!!! BUG

//    @Test
    public void simpleUpload() {
        Integer sampleid = mkDir();
        DXJob.Describe upload1 = upload(PSYDUCK, sampleLeft, sampleid, "left");

        String file1Id = upload1.getOutput(JsonNode.class).get("file").get("$dnanexus_link").asText();

        ResponseEntity<JsonNode> describe = getForResponseEnity(PSYDUCK, String.format("/sample/%s/describe/", file1Id));
        log.debug("D: " + describe);
    }

//    @Test
    public void upload() {
        Integer sampleId = mkDir();
        DXJob.Describe upload1 = upload(PSYDUCK, sampleLeft, sampleId, "left");
        DXJob.Describe upload2 = upload(PSYDUCK, sampleRight, sampleId, "right");

        ResponseEntity<JsonNode> ls = getForResponseEnity(PSYDUCK, String.format("/sample/%s/ls", sampleId));
        assertTrue(ls.getStatusCode().is2xxSuccessful());
        ResponseEntity<JsonNode> lsByNames = getForResponseEnity(PSYDUCK, String.format("/sample/%s/ls?byNames=true", sampleId));
        assertTrue(lsByNames.getStatusCode().is2xxSuccessful());
        log.info(ls.getBody());
        log.info(lsByNames.getBody());

        String file1Id = upload1.getOutput(JsonNode.class).get("file").get("$dnanexus_link").asText();

        ResponseEntity<JsonNode> describe = getForResponseEnity(PSYDUCK, String.format("/sample/%s/describe/", file1Id));
        log.debug("D: " + describe);
        assertTrue(describe.getStatusCode().is2xxSuccessful());
        log.info(describe.getBody());

//        log.info(fastqc(PSYDUCK, file1Id));
//        log.info(bwa(PSYDUCK, sampleId));
//        log.info(gatkhc(PSYDUCK, sampleId, interval));

    }

    @Test
    public void usersGroupsTest() {
        ResponseEntity<JsonNode> admin = getForResponseEnity(ADMIN, "/groups");
        assertTrue(admin.getStatusCode().is2xxSuccessful());
        JsonNode node = admin.getBody();
        assertEquals(node.size(), 2);
    }


}