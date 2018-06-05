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
public class ParseTextControllerTest {
	
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
		String query = "Warburg et al. (1993) used the designation Micro syndrome for an autosomal recessive syndrome comprising microcephaly, microcornea, congenital cataract, mental retardation, optic atrophy, and hypogenitalism. They described an affected brother and sister and their male cousin. The sibs were offspring of a consanguineous Pakistani marriage; the parents of the cousin denied consanguinity. Agenesis of the corpus callosum, prominent root of the nose, large anteverted ears, facial hypertrichosis, small pupils with posterior synechiae, hypotonia, mild to moderate spastic palsy with hip dislocations, and hormonal dysfunction, presumably of hypothalamic origin, were other features. The children were almost blind, whether or not the cataracts had been operated on. The electroretinographic responses indicated dysfunction of both retinal rods and cones, and the visual evoked potentials confirmed optic nerve atrophy. The children were late walkers and were incontinent of urine and stools. In the differential diagnosis, Warburg et al. (1993) considered COFS syndrome (214150), CAMAK/CAMFAK syndromes (212540), Martsolf syndrome (212720), lethal Rutledge syndrome (270400), and lethal Neu-Laxova syndrome (256520).  ";
		String requestJson = "{\"query\":\"" + query + "\"}";
		HttpHeaders headers = new HttpHeaders();
		headers.setContentType(MediaType.APPLICATION_JSON);

		HttpEntity<String> entity = new HttpEntity<String>(requestJson,headers);
		String body = this.restTemplate.postForObject(url, entity, String.class);
		System.out.println(body);
	}

	@Test
	public void diseaseAutocomplete() {
		String url = "/disease-autocomplete?firstLetters=mic";
		//String firsLetters = "mic";
		//String requestJson = "{\"query\":\"" + query + "\"}";
		//HttpHeaders headers = new HttpHeaders();
		//headers.setContentType(MediaType.APPLICATION_JSON);

		//HttpEntity<String> entity = new HttpEntity<String>(requestJson,headers);
		//String body = this.restTemplate.postForObject(url, entity, String.class);
		String body = this.restTemplate.getForObject(url, String.class);
		System.out.println(body);
	}
	
}