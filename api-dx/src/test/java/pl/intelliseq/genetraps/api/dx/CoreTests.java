package pl.intelliseq.genetraps.api.dx;

import com.fasterxml.jackson.databind.JsonNode;
import lombok.extern.log4j.Log4j2;

import static org.hamcrest.CoreMatchers.is;

import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.test.context.SpringBootTest;
import org.springframework.boot.test.context.SpringBootTest.WebEnvironment;
import org.springframework.boot.test.web.client.TestRestTemplate;
import org.springframework.core.io.Resource;
import org.springframework.http.HttpEntity;
import org.springframework.http.HttpHeaders;
import org.springframework.http.HttpMethod;
import org.springframework.http.ResponseEntity;
import org.springframework.test.context.junit4.SpringRunner;
import org.springframework.test.context.testng.AbstractTestNGSpringContextTests;
import org.testng.annotations.*;
import pl.intelliseq.genetraps.api.dx.helpers.AuroraDBManager;

import java.io.IOException;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertThat;
import static org.junit.Assert.assertTrue;
import static pl.intelliseq.genetraps.api.dx.TestUser.DEVIL;
import static pl.intelliseq.genetraps.api.dx.TestUser.PSYDUCK;


//@RunWith(SpringRunner.class)
@SpringBootTest(webEnvironment = WebEnvironment.RANDOM_PORT)
@Log4j2
//@ActiveProfiles("test")
public class CoreTests extends AbstractTestNGSpringContextTests {

    @Autowired
    private TestRestTemplate restTemplate;

    @Autowired
    private AuroraDBManager auroraDBManager;

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
    public void badge(){
        Resource resource = restTemplate.getForObject("/status", Resource.class);
        log.info(resource.exists());
    }

    @Test
    public void user() {

        HttpHeaders headers = new HttpHeaders();
        headers.set("Authorization", "Bearer "+PSYDUCK.getAccessToken());

        HttpEntity<String> entity = new HttpEntity<>("parameters", headers);

        ResponseEntity<JsonNode> response = restTemplate.exchange("/user", HttpMethod.GET, entity, JsonNode.class);

        log.debug(response);

        assertTrue(response.getStatusCode().is2xxSuccessful());
        assertEquals(response.getBody().get("id").asInt(), PSYDUCK.getId().intValue());


    }

    @Test
    public void simpleUserTest(){
        User psyduck = auroraDBManager.getUserDetails(PSYDUCK.getId());
        log.info(psyduck);
        log.info(auroraDBManager.getUserPrivileges(psyduck.getId()));
        log.info(auroraDBManager.getUserPrivileges(DEVIL.getId()));
        log.info(auroraDBManager.getUserPrivilegesToSample(1,1));
    }

//    @Test
//    public void rootPrivileges(){
//        log.info(auroraDBManager.getRootPrivileges());
//    }

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