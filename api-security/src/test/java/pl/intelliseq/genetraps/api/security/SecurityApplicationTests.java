package pl.intelliseq.genetraps.api.security;

import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.mashape.unirest.http.HttpResponse;
import com.mashape.unirest.http.Unirest;
import com.mashape.unirest.http.exceptions.UnirestException;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.test.context.SpringBootTest;
import org.springframework.core.env.Environment;
import org.springframework.test.context.junit4.SpringRunner;

import java.io.IOException;


@RunWith(SpringRunner.class)
@SpringBootTest(webEnvironment=SpringBootTest.WebEnvironment.RANDOM_PORT)
public class SecurityApplicationTests {

    @Autowired
    Environment env;

	@Test
	public void contextLoads() {
	}

	@Test
	public void test() throws UnirestException, IOException {
		HttpResponse<String> response = Unirest.post(String.format("http://localhost:%s/oauth/token", env.getProperty("local.server.port")))
				.header("content-type", "application/x-www-form-urlencoded")
				.header("authorization", "Basic d2ViX2FwcDpzZWNyZXQ=")
				.body("grant_type=password&username=admin&password=welcome1")
				.asString();

		ObjectMapper mapper = new ObjectMapper();
		JsonNode actualObj = mapper.readTree(response.getBody());

		assert actualObj.get("access_token").asText().length() > 0;
	}

}
