package pl.intelliseq.genetraps.api.dx;

import org.apache.log4j.Logger;
import org.hamcrest.Matchers;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;
import org.json.simple.parser.ParseException;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.test.context.SpringBootTest;
import org.springframework.boot.test.context.SpringBootTest.WebEnvironment;
import org.springframework.boot.test.web.client.TestRestTemplate;
import org.springframework.http.HttpEntity;
import org.springframework.http.HttpHeaders;
import org.springframework.http.MediaType;
import org.springframework.test.context.junit4.SpringRunner;
import org.springframework.util.LinkedMultiValueMap;
import org.springframework.util.MultiValueMap;
import pl.intelliseq.genetraps.api.dx.enums.JobStates;

import java.util.HashMap;
import java.util.Map;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertThat;


@RunWith(SpringRunner.class)
@SpringBootTest(webEnvironment = WebEnvironment.RANDOM_PORT)
public class UploadTests {

	Logger log = Logger.getLogger(UploadTests.class);

	@Autowired
	private TestRestTemplate restTemplate;

    public String sampleUrl = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeGencodeV11/supplemental/gencode.v11.tRNAs.gtf.gz";

    @Test
    public void upload(){
        JSONParser jsonParser = new JSONParser();


        HttpHeaders uploadHeaders = new HttpHeaders();
        uploadHeaders.setContentType(MediaType.APPLICATION_FORM_URLENCODED);

        MultiValueMap<String, String> uploadValueMap = new LinkedMultiValueMap<String, String>();
        uploadValueMap.add("url", sampleUrl);
        uploadValueMap.add("sampleNumber", "1");

        HttpEntity<MultiValueMap<String, String>> uploadEntity = new HttpEntity<MultiValueMap<String, String>>(uploadValueMap, uploadHeaders);

        String response = this.restTemplate.postForObject("/upload", uploadEntity, String.class);
        System.out.println(response);
//        String jobId = response.substring(7,response.length()-2);

        String jobId = null;
        try {
            jobId = (String) ((JSONObject) jsonParser.parse(response)).get("id");
        } catch (ParseException e) {
            throw new RuntimeException(e);
        }

        assertThat(jobId, Matchers.matchesPattern("job-\\w*"));

        /**
         * Check state
         */

        try{
            JobStates jobState;
            do {
                response = this.restTemplate.getForObject("/state/{jobId}", String.class, jobId);
                String state = (String) ((JSONObject) jsonParser.parse(response)).get("state");
                jobState = JobStates.getEnum(state);
                Thread.sleep(5000);
            }while (jobState != JobStates.DONE) ;
        }catch (InterruptedException | ParseException e){
            throw new RuntimeException(e);
        }

    }


	
}
