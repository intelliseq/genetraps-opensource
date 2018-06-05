package pl.intelliseq.explorare;

import org.junit.Test;
import org.junit.runner.RunWith;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.test.context.SpringBootTest;
import org.springframework.boot.test.web.client.TestRestTemplate;
import org.springframework.http.HttpEntity;
import org.springframework.http.HttpHeaders;
import org.springframework.http.MediaType;
import org.springframework.test.context.junit4.SpringJUnit4ClassRunner;



import static org.junit.Assert.*;

/**
 * @author Marcin Piechota
 */
@RunWith(SpringJUnit4ClassRunner.class)
@SpringBootTest(webEnvironment = SpringBootTest.WebEnvironment.RANDOM_PORT)
public class GetGenesControllerTest {
	
	@Autowired
	private TestRestTemplate restTemplate;

	@Test
	public void getDiseasesTest() {
		String url = "/get-diseases";
		String query = "[{\"id\":\"HP:0000006\",\"name\":\"Autosomal dominant inheritance\"},{\"id\":\"HP:0012188\",\"name\":\"Hyperemesis gravidarum\"},{\"id\":\"HP:0000836\",\"name\":\"Hyperthyroidism\"}]";
		String requestJson = "{\"hpoTerms\":" + query + "}";

		HttpHeaders headers = new HttpHeaders();
		headers.setContentType(MediaType.APPLICATION_JSON);

		HttpEntity<String> entity = new HttpEntity<String>(requestJson,headers);
		String body = this.restTemplate.postForObject(url, entity, String.class);
		assertTrue(body.contains("HYPERTHYROIDISM, FAMILIAL GESTATIONAL"));
		
		
		query = "[{\"id\":\"HP:0000836\",\"name\":\"Hyperthyroidism\"}]";
		requestJson = "{\"hpoTerms\":" + query + "}";
		entity = new HttpEntity<String>(requestJson,headers);
		body = this.restTemplate.postForObject(url, entity, String.class);
		System.out.println(body);
	}
	
}