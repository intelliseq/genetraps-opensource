package pl.intelliseq.genetraps.api.dx;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.IOException;
import org.apache.log4j.Logger;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.test.context.SpringBootTest;
import org.springframework.test.context.junit4.SpringRunner;
import org.springframework.util.LinkedMultiValueMap;
import org.springframework.util.MultiValueMap;

import pl.intelliseq.genetraps.api.dx.parser.DxJob;
import pl.intelliseq.genetraps.api.dx.parser.DxRunner;
import org.springframework.boot.test.context.SpringBootTest.WebEnvironment;
import org.springframework.boot.test.web.client.TestRestTemplate;

@RunWith(SpringRunner.class)
@SpringBootTest(webEnvironment = WebEnvironment.RANDOM_PORT)
public class ApiDxApplicationTests {

	Logger log = Logger.getLogger(ApiDxApplicationTests.class);

	@Autowired
	private TestRestTemplate restTemplate;
	
	@Test
	public void contextLoads() {
	}
	
	@Test
	public void testVep() {
		
		//MultiValueMap<String, Object> map = new LinkedMultiValueMap<String, Object>();
		//map.add("variant", "rs669");
		
		String body = this.restTemplate.getForObject("/touch", String.class);
		//String body = this.restTemplate.postForObject("/touch", map, String.class);
		System.out.println(body);
		assertThat(body.contains("Job ID"));
	}
	
	//@Test
	public void dxTest() throws IOException, InterruptedException {
    	String result = DxRunner.runCommand("printf \"Y\\n\" | dx run touch -iname=test");
    	
    	//log.info("Result: " + result);
    	//log.info("Result: " + this.getJobId(result + "\ntadam"));
    	
    	String jobId = this.getJobId(result);
    	
    	//log.info(DxJob.getDxJobById(jobId)); 
    	
    	for (int i = 1; i <= 30; i++) {
    		log.info("i: " + i);
    		Thread.sleep(500);
    		log.info(DxJob.getDxJobById(jobId));
    	}
	}

	private String getJobId(String result) {
		try {
			int indexOfJobId = result.indexOf("Job ID") + 8;
			int endIndexOfJobId = result.substring(indexOfJobId).indexOf("\n");
			if (endIndexOfJobId == -1) {
				return result.substring(indexOfJobId);
			}
			return result.substring(indexOfJobId, indexOfJobId + endIndexOfJobId);
		} catch (Exception e) {
			return null;
		}
	}
	
}
