package pl.intelliseq.genetraps.api.dx;

import org.apache.log4j.Logger;
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

import java.util.HashMap;
import java.util.Map;

import static org.junit.Assert.assertEquals;


@RunWith(SpringRunner.class)
@SpringBootTest(webEnvironment = WebEnvironment.RANDOM_PORT)
public class UploadTests {

	Logger log = Logger.getLogger(UploadTests.class);

	@Autowired
	private TestRestTemplate restTemplate;

    public String sampleUrl = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeGencodeV11/supplemental/gencode.v11.tRNAs.gtf.gz";

    @Test
    public void upload() {
        HttpHeaders uploadHeaders = new HttpHeaders();
        uploadHeaders.setContentType(MediaType.APPLICATION_FORM_URLENCODED);

        MultiValueMap<String, String> uploadValueMap = new LinkedMultiValueMap<String, String>();
        uploadValueMap.add("url", sampleUrl);

        HttpEntity<MultiValueMap<String, String>> uploadEntity = new HttpEntity<MultiValueMap<String, String>>(uploadValueMap, uploadHeaders);

        String response = this.restTemplate.postForObject("/upload", uploadEntity, String.class);
        System.out.println(response);
        System.out.println(response.substring(7,response.length()-2));

        /**
         * Check state
         */

        response = this.restTemplate.getForObject("/state/{jobId}", String.class, response.substring(7,response.length()-2));
        System.out.println(response);
    }


	
}
