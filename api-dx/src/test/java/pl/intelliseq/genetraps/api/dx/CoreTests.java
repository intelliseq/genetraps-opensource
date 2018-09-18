package pl.intelliseq.genetraps.api.dx;

import com.fasterxml.jackson.databind.JsonNode;
import lombok.extern.log4j.Log4j2;

import static org.hamcrest.CoreMatchers.is;
import static org.hamcrest.text.IsEqualIgnoringCase.equalToIgnoringCase;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.test.context.SpringBootTest;
import org.springframework.boot.test.context.SpringBootTest.WebEnvironment;
import org.springframework.boot.test.web.client.TestRestTemplate;
import org.springframework.http.HttpEntity;
import org.springframework.http.HttpHeaders;
import org.springframework.http.HttpMethod;
import org.springframework.http.ResponseEntity;
import org.springframework.test.context.junit4.SpringRunner;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertThat;
import static org.junit.Assert.assertTrue;
import static pl.intelliseq.genetraps.api.dx.TestUser.PSYDUCK;


@RunWith(SpringRunner.class)
@SpringBootTest(webEnvironment = WebEnvironment.RANDOM_PORT)
@Log4j2
//@ActiveProfiles("test")
public class CoreTests {

    @Autowired
    private TestRestTemplate restTemplate;

    @Test
    public void empty() {
        assert true;
    }

    @Test
    public void helloTest() {
        String response = this.restTemplate.getForObject("/hello", String.class);

        assertEquals(response, "{\"status\": \"up\"}");
    }

    @Test
    public void user() {

        HttpHeaders headers = new HttpHeaders();
        headers.set("Authorization", "Bearer "+PSYDUCK.getAccessToken());

        HttpEntity<String> entity = new HttpEntity<>("parameters", headers);

        ResponseEntity<JsonNode> response = restTemplate.exchange("/user", HttpMethod.GET, entity, JsonNode.class);

        assertTrue(response.getStatusCode().is2xxSuccessful());
        assertThat(response.getBody().get("Username").asText(), is(equalToIgnoringCase(PSYDUCK.toString())));


    }

    @Test
    public void logging(){
        log.fatal("FATAL");
        log.error("ERROR");
        log.warn("WARN");
        log.info("INFO");
        log.debug("DEBUG");
        log.trace("TRACE");
    }

}