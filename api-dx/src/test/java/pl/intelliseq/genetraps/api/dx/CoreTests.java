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

import static org.junit.Assert.assertEquals;


@RunWith(SpringRunner.class)
@SpringBootTest(webEnvironment = WebEnvironment.RANDOM_PORT)
public class CoreTests {

    Logger log = Logger.getLogger(CoreTests.class);

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

}
