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
public class GetPhenotypeTest {
	
	@Autowired
	private TestRestTemplate restTemplate;

	@Test
	public void diseaseAutocomplete() {
		String url = "/get-phenotype-by-id?phenotypeid=HP:0000836";
		String body = this.restTemplate.getForObject(url, String.class);
		System.out.println(body);
	}
	
}