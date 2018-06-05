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
public class RestControllerTest {
	
	@Autowired
	private TestRestTemplate restTemplate;

	@Test
	public void helloTest() {
		String body = this.restTemplate.getForObject("/hello", String.class);
		assertEquals(body.toString(), "{\"response\":\"explorare service is alive\"}");
	}
	
	@Test
	public void parseTextTest() {
		String url = "/parse-text";
		String query = "The patient is 15 y.o boy with no prior somatic problems. The first prodromal symptoms occurred in June, 2015. The patient had transient no provoked pain in lower extremities and in right lower quadrant of the abdomen. It resolved spontaneously. Then in August, 2015 a severe pain occurred in the back and again in the right lower quadrant of the abdomen. This time the patient was transferred to ER. ECHO showed dissection of the aorta type DeBackeyIII/Stanford B. The valves and Ascending aorta appeared normal. CT of the abdomen showed the dissection started 3cm beneath of the origin of the Left subclavian artery and ended 27 cm into the right common iliac artery; infarct of the lower pole of the right kidney; double supply of the right kidney –lower pole from the false channel of the aorta and the upper pole from the true channel of the aorta; left iliac common artery is supplied from the true channel. As so far, the patient is treated symptomatically. In addition, he has been treated for the past 10 years for obsessive-compulsive and ADHD type of behavior. Medications: perazine, fluoxetine, haroperidol (on demand). Physical examination beside macrocephaly, mild scoliosis, flat foot, and abdominal tenderness, was normal.";
		String requestJson = "{\"query\":\"" + query + "\"}";
		HttpHeaders headers = new HttpHeaders();
		headers.setContentType(MediaType.APPLICATION_JSON);

		HttpEntity<String> entity = new HttpEntity<String>(requestJson,headers);
		String body = this.restTemplate.postForObject(url, entity, String.class);
		System.out.println(body);
	}

	
}